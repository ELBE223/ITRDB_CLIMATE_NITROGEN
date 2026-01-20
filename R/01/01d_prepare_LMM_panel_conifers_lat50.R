###############################################################################
# 01d_prepare_LMM_panel_conifers_lat50.R
#
# Goal:
#   Final panel dataset for LMM (LAT>=50 from 01a):
#   - Keep conifers via 01c genus list (global > 15 sites)
#   - Join site-level metadata (elevation_m + optional genus fallback)
#   - Convert BAI to cm^2 and remove outliers (BAI_cm2 > 100)
#   - Output parquet with:
#     year, BAI_cm2, BAI, log_BAI, log_BAI_eps, age, continent, site_id, tree_id, tree_uid,
#     lat, lon, elevation_m, genus
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, dplyr, readr, stringr, fs, arrow, tibble)

# ----------------------------- CONFIG ----------------------------------------

CFG <- list(
  PATH_META_LAT50 = here::here("output", "01a_analyse_ITRDB", "tables", "01a_ITRDB_site_metadata_lat50.csv"),
  PATH_KEEP_GENUS = here::here("output", "01c_qc_descriptives_lat50_conifers", "01c_conifers_global_genus_filtered.csv"),
  
  # Optional hint (file or dir). Leave NULL for auto-detect.
  BAI_SOURCE_HINT = NULL,
  
  OUT_DIR  = here::here("output", "01d_prepare_LMM_panel"),
  OUT_FILE = "01d_LMM_panel_conifers_lat50_genus15.parquet",
  
  YEAR_MIN = 1901L,
  YEAR_MAX = 2010L,
  
  EPS = 1e-4,
  
  # Unit conversion + outlier filter
  BAI_IN_MM2 = TRUE,          # if TRUE: convert to cm^2 via /100
  BAI_CM2_MAX = 100           # drop rows with BAI_cm2 > 100
)

# ----------------------------- HELPERS ---------------------------------------

has_required_fields <- function(field_names) {
  core_ok <- all(c("site_id", "tree_id", "year") %in% field_names)
  bai_ok  <- any(c("bai_mm2", "BAI", "bai") %in% field_names)
  core_ok && bai_ok
}

open_bai_dataset <- function(x) {
  ds <- arrow::open_dataset(x, format = "parquet")
  f  <- names(ds$schema)
  if (!has_required_fields(f)) return(NULL)
  ds
}

find_bai_dataset <- function(base_out, hint = NULL) {
  if (!is.null(hint)) {
    if (fs::dir_exists(hint)) {
      ds <- open_bai_dataset(hint)
      if (!is.null(ds)) return(list(ds = ds, source = hint))
    }
    if (fs::file_exists(hint) && stringr::str_detect(hint, "\\.parquet$")) {
      ds <- open_bai_dataset(hint)
      if (!is.null(ds)) return(list(ds = ds, source = hint))
    }
  }
  
  all_parquet <- fs::dir_ls(base_out, recurse = TRUE, type = "file", glob = "*.parquet")
  if (length(all_parquet) == 0) stop("No parquet files found under output/. Did 01b write parquet?")
  
  cand_01b <- all_parquet[str_detect(str_to_lower(all_parquet), "01b")]
  cand_any <- all_parquet
  
  rank_paths <- function(paths) {
    p <- str_to_lower(paths)
    score <- 0L +
      ifelse(str_detect(p, "bai"), 4L, 0L) +
      ifelse(str_detect(p, "parts"), 2L, 0L) +
      ifelse(str_detect(p, "parquet"), 1L, 0L)
    paths[order(score, decreasing = TRUE)]
  }
  
  cand_01b <- rank_paths(cand_01b)
  cand_any <- rank_paths(cand_any)
  
  try_paths <- function(paths) {
    if (length(paths) == 0) return(NULL)
    
    parent_dirs <- unique(dirname(paths))
    parent_dirs <- rank_paths(parent_dirs)
    
    for (d in parent_dirs) {
      files_in_dir <- paths[dirname(paths) == d]
      if (length(files_in_dir) >= 2) {
        ds <- open_bai_dataset(d)
        if (!is.null(ds)) return(list(ds = ds, source = d))
      }
    }
    
    for (f in head(paths, 200)) {
      ds <- open_bai_dataset(f)
      if (!is.null(ds)) return(list(ds = ds, source = f))
    }
    
    NULL
  }
  
  res <- try_paths(cand_01b)
  if (!is.null(res)) return(res)
  
  res <- try_paths(cand_any)
  if (!is.null(res)) return(res)
  
  stop(
    "Could not find a parquet dataset with columns {site_id, tree_id, year, bai_*}. ",
    "Found parquet files (first 20):\n",
    paste(head(all_parquet, 20), collapse = "\n")
  )
}

# ----------------------------- PATHS -----------------------------------------

stopifnot(fs::file_exists(CFG$PATH_META_LAT50))
stopifnot(fs::file_exists(CFG$PATH_KEEP_GENUS))
fs::dir_create(CFG$OUT_DIR)
path_out <- fs::path(CFG$OUT_DIR, CFG$OUT_FILE)

cat("=== 01d Prepare LMM Panel ===\n")
cat("Meta lat50: ", CFG$PATH_META_LAT50, "\n", sep = "")
cat("Keep genus: ", CFG$PATH_KEEP_GENUS, "\n", sep = "")
cat("Output:     ", path_out, "\n\n", sep = "")

# ----------------------------- LOAD KEEP GENERA (01c) -------------------------

keep_tbl <- readr::read_csv(CFG$PATH_KEEP_GENUS, show_col_types = FALSE)
if (!("Genus" %in% names(keep_tbl))) stop("Keep-genus file missing column 'Genus': ", CFG$PATH_KEEP_GENUS)

keep_genera <- keep_tbl %>%
  mutate(Genus = str_to_title(as.character(Genus))) %>%
  pull(Genus) %>%
  unique() %>%
  na.omit()

cat("Genera kept (01c, global >15): ", length(keep_genera), "\n", sep = "")
cat("Kept: ", paste(keep_genera, collapse = ", "), "\n\n", sep = "")

# ----------------------------- LOAD META (01a LAT>=50) ------------------------

meta <- readr::read_csv(CFG$PATH_META_LAT50, show_col_types = FALSE)

meta_req <- c("site_id", "continent", "lat", "lon")
missing_meta <- setdiff(meta_req, names(meta))
if (length(missing_meta) > 0) stop("Meta file missing columns: ", paste(missing_meta, collapse = ", "))

meta2 <- meta %>%
  mutate(
    site_id = as.character(site_id),
    continent_meta = as.character(continent),
    lat_meta = as.numeric(lat),
    lon_meta = as.numeric(lon),
    elevation_m = if ("elevation_m" %in% names(.)) as.numeric(.data$elevation_m) else NA_real_,
    genus_meta = if ("genus" %in% names(.)) str_to_title(as.character(.data$genus)) else NA_character_
  ) %>%
  select(site_id, continent_meta, lat_meta, lon_meta, elevation_m, genus_meta)

# ----------------------------- LOAD BAI DATASET (01b) -------------------------

base_out <- here::here("output")
bai_res <- find_bai_dataset(base_out, hint = CFG$BAI_SOURCE_HINT)

cat("BAI source found: ", bai_res$source, "\n\n", sep = "")

ds <- bai_res$ds
schema_names <- names(ds$schema)

bai_col <- dplyr::case_when(
  "bai_mm2" %in% schema_names ~ "bai_mm2",
  "BAI" %in% schema_names ~ "BAI",
  "bai" %in% schema_names ~ "bai",
  TRUE ~ NA_character_
)
if (is.na(bai_col)) stop("Cannot find BAI column (expected bai_mm2/BAI/bai).")

age_col <- dplyr::case_when( # Alter/Index bezogen auf 1901
  "age" %in% schema_names ~ "age",
  "age_1901" %in% schema_names ~ "age_1901", 
  TRUE ~ NA_character_
)

select_cols <- c("site_id", "tree_id", "year", bai_col, "continent", "lat", "lon", "species_name")
if (!is.na(age_col)) select_cols <- c(select_cols, age_col)

bai <- ds %>%
  select(all_of(select_cols)) %>%
  filter(year >= CFG$YEAR_MIN, year <= CFG$YEAR_MAX) %>%
  collect()

cat("Rows loaded (year-filtered): ", nrow(bai), "\n", sep = "")
cat("Columns loaded: ", paste(names(bai), collapse = ", "), "\n\n", sep = "")

# ----------------------------- JOIN + BUILD FINAL ----------------------------

df <- bai %>%
  mutate(
    site_id = as.character(site_id),
    tree_id = as.character(tree_id),
    year = as.integer(year),
    BAI_raw = as.numeric(.data[[bai_col]]),
    continent = as.character(continent),
    lat = as.numeric(lat),
    lon = as.numeric(lon)
  ) %>%
  left_join(meta2, by = "site_id") %>%
  mutate(
    continent = coalesce(continent, continent_meta),
    lat = coalesce(lat, lat_meta),
    lon = coalesce(lon, lon_meta),
    
    genus = coalesce(
      genus_meta,
      str_to_title(str_extract(as.character(species_name), "^[A-Za-z]+"))
    )
  ) %>%
  filter(genus %in% keep_genera)

# Age: prefer existing (age or age_1901), else compute within tree
if (!is.na(age_col)) {
  df <- df %>% mutate(age = as.integer(.data[[age_col]]))
} else {
  df <- df %>%
    group_by(site_id, tree_id) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(age = row_number()) %>%
    ungroup()
}

# Convert BAI to cm^2 if needed + remove outliers
df <- df %>%
  mutate(
    BAI_cm2 = if (isTRUE(CFG$BAI_IN_MM2)) BAI_raw / 100 else BAI_raw
  ) %>%
  filter(!is.na(BAI_cm2), BAI_cm2 > 0, BAI_cm2 <= CFG$BAI_CM2_MAX)

# Logs + IDs (based on cm^2)
df <- df %>%
  mutate(
    log_BAI = log(BAI_cm2),
    log_BAI_eps = log(BAI_cm2 + CFG$EPS),
    tree_uid = paste(site_id, tree_id, sep = "__"),
    BAI = BAI_cm2  # optional alias for convenience
  )

final <- df %>%
  transmute(
    year,
    BAI_cm2,
    BAI,
    log_BAI,
    log_BAI_eps,
    age,
    continent,
    site_id,
    tree_id,
    tree_uid,
    lat,
    lon,
    elevation_m,
    genus
  ) %>%
  filter(
    !is.na(year),
    !is.na(site_id),
    !is.na(tree_id),
    !is.na(continent),
    !is.na(lat), !is.na(lon),
    !is.na(genus)
  )

# ----------------------------- QC PRINTS -------------------------------------

cat("=== QC (final) ===\n")
cat("Rows:       ", nrow(final), "\n", sep = "")
cat("Sites:      ", n_distinct(final$site_id), "\n", sep = "")
cat("Trees:      ", n_distinct(final$tree_uid), "\n", sep = "")
cat("Year range: ", min(final$year, na.rm = TRUE), " - ", max(final$year, na.rm = TRUE), "\n", sep = "")

cat("\nCounts by continent:\n")
print(table(final$continent, useNA = "ifany"))

cat("\nTop genera (final, top 15):\n")
print(final %>% count(genus, name = "N") %>% arrange(desc(N)) %>% head(15))

cat("\nMissing elevation_m: ", sum(is.na(final$elevation_m)), "\n", sep = "")

cat("\nBAI_cm2 summary:\n")
print(summary(final$BAI_cm2))

cat("\nPreview final (first 5 rows):\n")
print(head(final, 5))

# ----------------------------- SAVE PARQUET ----------------------------------

arrow::write_parquet(final, path_out, compression = "zstd")

cat("\n=== 01d DONE ===\n")
cat("Saved parquet: ", path_out, "\n", sep = "")
cat("Columns: ", paste(names(final), collapse = ", "), "\n", sep = "")
