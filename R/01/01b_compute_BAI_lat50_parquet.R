###############################################################################
# 01b_compute_BAI_lat50_parquet.R 
#
# Goal:
#   - Read 01a filtered site metadata (lat >= 50N, meta OK)
#   - Parse local RWL (prefer Tucson .rwl, fallback NOAA -noaa.rwl)
#   - Assign TreeID, compute BAI per year per tree, add log transforms + Age_1901
#   - Save ONLY Parquet outputs
#   - Print ONLY high-level progress + QC summaries (no per-file dplR chatter)
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, fs, dplyr, stringr, readr, purrr, tibble,
  tidyr, furrr, future, parallel,
  dplR, arrow
)

# ----------------------------- CONFIG ----------------------------------------

CFG <- list(
  PARALLEL = TRUE,
  MAX_WORKERS = 4,
  
  CHUNK_SIZE_SITES = 75L,
  
  SCALE_THRESHOLD_MEDIAN = 20,  # median rw > 20 -> treat as 0.01 mm
  EPS_LOG = 1e-4,
  YEAR_REF_AGE = 1901L,
  
  # console
  PRINT_CHUNK_SUMMARY = TRUE
)

num_cores <- max(1, parallel::detectCores() - 1)
workers <- if (isTRUE(CFG$PARALLEL)) min(num_cores, CFG$MAX_WORKERS) else 1
future::plan(future::multisession, workers = workers)

# ----------------------------- PATHS -----------------------------------------

step_00a <- "00a_download_ITRDB_RWL_global"
dir_00a  <- here::here("output", step_00a)
dir_rwl  <- fs::path(dir_00a, "data", "rwl")

step_01a <- "01a_analyse_ITRDB"
dir_01a  <- here::here("output", step_01a)
in_lat50 <- fs::path(dir_01a, "tables", "01a_ITRDB_site_metadata_lat50.csv")

step_01b <- "01b_analyse_ITRDB"
out_dir   <- here::here("output", step_01b)
dir_data  <- fs::path(out_dir, "data")
dir_parts <- fs::path(dir_data, "bai_parts")
dir_tab   <- fs::path(out_dir, "tables")
dir_dbg   <- fs::path(out_dir, "debug")
fs::dir_create(c(dir_data, dir_parts, dir_tab, dir_dbg))

stopifnot(fs::dir_exists(dir_rwl))
stopifnot(fs::file_exists(in_lat50))

cat("=== 01b Compute BAI (Parquet) ===\n")
cat("Input lat50: ", in_lat50, "\n", sep = "")
cat("RWL dir:     ", dir_rwl, "\n", sep = "")
cat("Workers:     ", workers, "\n", sep = "")
cat("Chunk size:  ", CFG$CHUNK_SIZE_SITES, " sites\n\n", sep = "")

# ----------------------------- LOAD SITES ------------------------------------

sites <- readr::read_csv(in_lat50, show_col_types = FALSE) %>% distinct()

req_cols <- c("site_id", "region", "continent", "lat", "lon", "species_name", "rwl_tucson_file", "rwl_noaa_file")
missing_cols <- setdiff(req_cols, names(sites))
if (length(missing_cols) > 0) stop("01a lat50 is missing columns: ", paste(missing_cols, collapse = ", "))

cat("Sites (lat>=50, meta OK): ", nrow(sites), "\n", sep = "")
cat("Site columns: ", paste(names(sites), collapse = ", "), "\n\n", sep = "")

# ----------------------------- HELPERS ---------------------------------------

quiet_read_rwl <- function(path) {
  x <- NULL
  capture.output({
    x <- suppressMessages(suppressWarnings(dplR::read.rwl(path)))
  })
  x
}

quiet_read_tucson <- function(path) {
  y <- NULL
  capture.output({
    y <- suppressMessages(suppressWarnings(dplR::read.tucson(path)))
  })
  y
}

safe_read_rwl <- function(path) {
  if (is.na(path) || !fs::file_exists(path)) return(list(ok = FALSE, data = NULL, note = "Missing file"))
  tryCatch({
    x <- quiet_read_rwl(path)
    list(ok = TRUE, data = x, note = NA_character_)
  }, error = function(e) {
    # optional fallback
    y <- tryCatch(quiet_read_tucson(path), error = function(e2) NULL)
    if (!is.null(y)) return(list(ok = TRUE, data = y, note = "Fallback read.tucson"))
    list(ok = FALSE, data = NULL, note = paste0("Read error: ", e$message))
  })
}

choose_and_read_site <- function(row) {
  tucson_path <- if (!is.na(row$rwl_tucson_file) && row$rwl_tucson_file != "") {
    fs::path(dir_rwl, row$region, row$rwl_tucson_file)
  } else NA_character_
  
  noaa_path <- if (!is.na(row$rwl_noaa_file) && row$rwl_noaa_file != "") {
    fs::path(dir_rwl, row$region, row$rwl_noaa_file)
  } else NA_character_
  
  if (!is.na(tucson_path) && fs::file_exists(tucson_path)) {
    r1 <- safe_read_rwl(tucson_path)
    if (isTRUE(r1$ok)) {
      return(list(ok = TRUE, data = r1$data, chosen_type = "tucson_rwl", chosen_path = as.character(tucson_path), note = r1$note))
    }
  }
  
  if (!is.na(noaa_path) && fs::file_exists(noaa_path)) {
    r2 <- safe_read_rwl(noaa_path)
    if (isTRUE(r2$ok)) {
      return(list(ok = TRUE, data = r2$data, chosen_type = "noaa_rwl", chosen_path = as.character(noaa_path), note = r2$note))
    }
  }
  
  note <- paste(
    "tucson:", ifelse(is.na(tucson_path), "NA", as.character(tucson_path)),
    "| noaa:", ifelse(is.na(noaa_path), "NA", as.character(noaa_path)),
    sep = " "
  )
  list(ok = FALSE, data = NULL, chosen_type = NA_character_, chosen_path = NA_character_, note = note)
}

rwl_to_long <- function(site_id, rwl_df) {
  yrs <- suppressWarnings(as.integer(rownames(rwl_df)))
  if (all(is.na(yrs))) {
    return(tibble(site_id = site_id, series = character(), year = integer(), rw_raw = numeric()))
  }
  
  tibble(year = yrs) %>%
    bind_cols(tibble::as_tibble(rwl_df, .name_repair = "minimal")) %>%
    tidyr::pivot_longer(
      cols = -year,
      names_to = "series",
      values_to = "rw_raw",
      values_drop_na = FALSE
    ) %>%
    mutate(
      site_id = site_id,
      series = trimws(series)
    ) %>%
    select(site_id, series, year, rw_raw)
}

infer_scale_mm <- function(rw_raw_vec) {
  x <- rw_raw_vec
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  if (length(x) == 0) return(1)
  med <- median(x)
  if (is.na(med)) return(1)
  if (med > CFG$SCALE_THRESHOLD_MEDIAN) return(0.01)
  1
}

compute_bai_tree <- function(df_tree) {
  df_tree <- df_tree %>% arrange(year)
  
  rw <- df_tree$rw_mm
  rw_clean <- ifelse(is.na(rw) | !is.finite(rw) | rw <= 0, NA_real_, rw)
  
  rw_fill <- ifelse(is.na(rw_clean), 0, rw_clean)
  radius <- cumsum(rw_fill)
  
  bai <- pi * (radius^2 - dplyr::lag(radius)^2)
  bai <- ifelse(is.na(rw_clean), NA_real_, bai)
  bai <- ifelse(is.finite(bai) & bai > 0, bai, NA_real_)
  
  first_year <- suppressWarnings(min(df_tree$year[!is.na(rw_clean)], na.rm = TRUE))
  if (!is.finite(first_year)) first_year <- NA_integer_
  
  age_1901 <- if (!is.na(first_year) && first_year <= CFG$YEAR_REF_AGE) {
    as.integer(CFG$YEAR_REF_AGE - first_year + 1L)
  } else {
    NA_integer_
  }
  
  df_tree %>%
    mutate(
      rw_mm = rw_clean,
      radius_mm = radius,
      bai = bai,
      log_bai = ifelse(!is.na(bai) & bai > 0, log(bai), NA_real_),
      log_bai_eps = log(dplyr::coalesce(bai, 0) + CFG$EPS_LOG),
      age_1901 = age_1901
    )
}

process_one_site <- function(row) {
  sid <- row$site_id
  
  rr <- choose_and_read_site(row)
  
  if (!isTRUE(rr$ok)) {
    qc_site <- tibble(
      site_id = sid,
      region = row$region,
      continent = row$continent,
      chosen_type = NA_character_,
      chosen_path = NA_character_,
      read_ok = FALSE,
      note = rr$note,
      n_trees = 0L,
      n_rows = 0L,
      scale_used = NA_real_
    )
    return(list(bai = NULL, qc_site = qc_site, qc_tree = NULL))
  }
  
  long <- rwl_to_long(sid, rr$data)
  
  scale <- infer_scale_mm(long$rw_raw)
  long <- long %>%
    mutate(
      rw_raw = suppressWarnings(as.numeric(rw_raw)),
      rw_mm  = rw_raw * scale
    ) %>%
    left_join(
      tibble(
        site_id = row$site_id,
        region = row$region,
        continent = row$continent,
        lat = row$lat,
        lon = row$lon,
        species_name = row$species_name
      ),
      by = "site_id"
    ) %>%
    mutate(
      tree_id = paste0(site_id, "__", series)
    )
  
  bai_long <- long %>%
    group_by(site_id, tree_id, series, region, continent, lat, lon, species_name) %>%
    group_modify(~ compute_bai_tree(.x)) %>%
    ungroup() %>%
    mutate(
      chosen_type = rr$chosen_type,
      chosen_path = rr$chosen_path,
      scale_used = scale
    ) %>%
    select(
      site_id, tree_id, series, year,
      continent, region, lat, lon, species_name,
      chosen_type, scale_used,
      rw_mm, bai, log_bai, log_bai_eps, age_1901
    )
  
  qc_tree <- bai_long %>%
    group_by(site_id, tree_id, continent, region, chosen_type, scale_used) %>%
    summarise(
      start_year = suppressWarnings(min(year, na.rm = TRUE)),
      end_year   = suppressWarnings(max(year, na.rm = TRUE)),
      n_years_total = n(),
      n_rw_missing  = sum(is.na(rw_mm)),
      n_bai_missing = sum(is.na(bai)),
      n_bai_pos     = sum(!is.na(bai)),
      age_1901 = suppressWarnings(max(age_1901, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(age_1901 = ifelse(is.finite(age_1901), as.integer(age_1901), NA_integer_))
  
  qc_site <- tibble(
    site_id = sid,
    region = row$region,
    continent = row$continent,
    chosen_type = rr$chosen_type,
    chosen_path = rr$chosen_path,
    read_ok = TRUE,
    note = rr$note,
    n_trees = dplyr::n_distinct(bai_long$tree_id),
    n_rows = nrow(bai_long),
    scale_used = scale
  )
  
  list(bai = bai_long, qc_site = qc_site, qc_tree = qc_tree)
}

write_part_parquet <- function(df, part_idx) {
  out_path <- fs::path(dir_parts, sprintf("01b_bai_long_part_%04d.parquet", part_idx))
  arrow::write_parquet(df, out_path, compression = "zstd")
  out_path
}

# ----------------------------- RUN (CHUNKED) ---------------------------------

site_rows <- sites %>% arrange(site_id)
n_sites <- nrow(site_rows)
chunks <- split(site_rows, ceiling(seq_len(n_sites) / CFG$CHUNK_SIZE_SITES))

cat("Chunks: ", length(chunks), "\n", sep = "")
cat("Output BAI columns (Parquet):\n")
cat(paste(
  c("site_id","tree_id","series","year","continent","region","lat","lon","species_name",
    "chosen_type","scale_used","rw_mm","bai","log_bai","log_bai_eps","age_1901"),
  collapse = ", "
), "\n\n", sep = "")

qc_sites_all <- list()
qc_trees_all <- list()
part_files <- character()

for (i in seq_along(chunks)) {
  chunk_df <- chunks[[i]]
  
  res_list <- furrr::future_pmap(
    .l = chunk_df,
    .f = function(...) process_one_site(tibble(...)),
    .options = furrr::furrr_options(seed = TRUE)
  )
  
  bai_chunk <- dplyr::bind_rows(purrr::map(res_list, "bai"))
  qc_sites_chunk <- dplyr::bind_rows(purrr::map(res_list, "qc_site"))
  qc_trees_chunk <- dplyr::bind_rows(purrr::map(res_list, "qc_tree"))
  
  qc_sites_all[[i]] <- qc_sites_chunk
  qc_trees_all[[i]] <- qc_trees_chunk
  
  if (nrow(bai_chunk) > 0) {
    p <- write_part_parquet(bai_chunk, i)
    part_files <- c(part_files, as.character(p))
  }
  
  if (isTRUE(CFG$PRINT_CHUNK_SUMMARY)) {
    ok_n <- sum(qc_sites_chunk$read_ok, na.rm = TRUE)
    fail_n <- sum(!qc_sites_chunk$read_ok, na.rm = TRUE)
    trees_n <- sum(qc_sites_chunk$n_trees, na.rm = TRUE)
    rows_n <- sum(qc_sites_chunk$n_rows, na.rm = TRUE)
    cat(sprintf("Chunk %d/%d | sites=%d | ok=%d fail=%d | trees=%d | rows=%d\n",
                i, length(chunks), nrow(chunk_df), ok_n, fail_n, trees_n, rows_n))
  }
}

qc_sites <- dplyr::bind_rows(qc_sites_all)
qc_trees <- dplyr::bind_rows(qc_trees_all)

failures <- qc_sites %>%
  filter(!read_ok) %>%
  select(site_id, continent, region, chosen_type, chosen_path, note)

# ----------------------------- SAVE QC (PARQUET) ------------------------------

out_qc_sites <- fs::path(dir_tab, "01b_qc_sites.parquet")
out_qc_trees <- fs::path(dir_tab, "01b_qc_trees.parquet")
out_fail     <- fs::path(dir_tab, "01b_read_failures.parquet")

arrow::write_parquet(qc_sites, out_qc_sites, compression = "zstd")
arrow::write_parquet(qc_trees, out_qc_trees, compression = "zstd")
arrow::write_parquet(failures, out_fail, compression = "zstd")

# ----------------------------- FINAL PRINTS -----------------------------------

cat("\n=== 01b DONE ===\n")
cat("Sites processed: ", nrow(qc_sites), "\n", sep = "")
cat("Read OK:         ", sum(qc_sites$read_ok, na.rm = TRUE), "\n", sep = "")
cat("Read FAIL:       ", sum(!qc_sites$read_ok, na.rm = TRUE), "\n", sep = "")
cat("Trees total:     ", sum(qc_sites$n_trees, na.rm = TRUE), "\n", sep = "")
cat("Rows total:      ", sum(qc_sites$n_rows, na.rm = TRUE), "\n", sep = "")
cat("Parquet parts:   ", length(part_files), "\n", sep = "")

cat("\nChosen type:\n")
print(table(qc_sites$chosen_type, useNA = "ifany"))

cat("\nScale used (mm factor):\n")
print(table(qc_sites$scale_used, useNA = "ifany"))

cat("\nSaved Parquet:\n")
cat("  BAI parts dir: ", dir_parts, "\n", sep = "")
cat("  qc_sites:      ", out_qc_sites, "\n", sep = "")
cat("  qc_trees:      ", out_qc_trees, "\n", sep = "")
cat("  failures:      ", out_fail, "\n", sep = "")

cat("\nParquet QC columns:\n")
cat("qc_sites: ", paste(names(qc_sites), collapse = ", "), "\n", sep = "")
cat("qc_trees: ", paste(names(qc_trees), collapse = ", "), "\n", sep = "")
cat("failures: ", paste(names(failures), collapse = ", "), "\n", sep = "")

if (nrow(failures) > 0) {
  cat("\nTop failures (first 10):\n")
  print(failures %>% head(10))
}
