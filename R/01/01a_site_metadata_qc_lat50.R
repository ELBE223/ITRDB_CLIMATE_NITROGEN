###############################################################################
# 01a_site_metadata_qc_lat50.R
#
# Goal:
#   Build clean ITRDB site metadata table from LOCAL NOAA metadata TXT files:
#     output/00a_download_ITRDB_RWL_global/data/meta_txt/<region>/*-rwl-noaa.txt
#   Assign SiteID + Continent, run QC flags, and export:
#     - FULL metadata table (for 01b)
#     - FILTERED table (lat >= 50N + meta OK + valid lat/lon)
#
# Output:
#   output/01a_analyse_ITRDB/tables/
#     01a_ITRDB_site_metadata.csv
#     01a_ITRDB_site_metadata_lat50.csv
#     01a_meta_failures.csv
#     01a_missing_latlon.csv
#     01a_qc_summary.csv
#
# Updates:
#   1) Russia split by longitude: rus* lon < 60 => Europe else Asia (no downloads)
#   2) Hemisphere enforcement by FIELD meaning:
#        - Northernmost lat always +, Southernmost always -
#        - Easternmost lon always +, Westernmost always -
#      This fixes cases where NOAA header values miss N/S/E/W letters (e.g. Alaska).
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, fs, dplyr, stringr, readr, purrr, tibble, furrr, future, parallel
)

# ----------------------------- CONFIG ----------------------------------------

CFG <- list(
  MODE = "test_then_all",      # "test_then_all" | "all"
  N_TEST = 100L,
  AUTO_RUN_ALL_IF_PASS = TRUE,
  
  PARALLEL = FALSE,
  MAX_WORKERS = 6,
  
  LAT_MIN_KEEP = 50,
  
  MIN_META_OK_RATE = 0.95,
  MIN_LATLON_RATE  = 0.90,
  MIN_SPECIES_RATE = 0.85
)

num_cores <- max(1, parallel::detectCores() - 1)
workers <- if (isTRUE(CFG$PARALLEL)) min(num_cores, CFG$MAX_WORKERS) else 1

if (isTRUE(CFG$PARALLEL)) {
  future::plan(future::multisession, workers = workers)
} else {
  future::plan(future::sequential)
}

# ----------------------------- PATHS -----------------------------------------

step_00a <- "00a_download_ITRDB_RWL_global"
in_dir   <- here::here("output", step_00a)
in_idx   <- fs::path(in_dir, "tables", "00a_ITRDB_download_status.csv")

dir_rwl      <- fs::path(in_dir, "data", "rwl")
dir_meta_txt <- fs::path(in_dir, "data", "meta_txt")

step_01a <- "01a_analyse_ITRDB"
out_dir  <- here::here("output", step_01a)
dir_tab  <- fs::path(out_dir, "tables")
dir_dbg  <- fs::path(out_dir, "debug")
fs::dir_create(c(dir_tab, dir_dbg))

stopifnot(fs::file_exists(in_idx))
stopifnot(fs::dir_exists(dir_rwl))
stopifnot(fs::dir_exists(dir_meta_txt))

cat("=== 01a Site Metadata + QC (OFFLINE) ===\n")
cat("Mode:     ", CFG$MODE, "\n", sep = "")
cat("Index:    ", in_idx, "\n", sep = "")
cat("RWL dir:  ", dir_rwl, "\n", sep = "")
cat("META dir: ", dir_meta_txt, "\n", sep = "")
cat("Workers:  ", workers, "\n\n", sep = "")

# ----------------------------- INDEX -> SITES --------------------------------

files_df <- readr::read_csv(in_idx, show_col_types = FALSE)

required_cols <- c("region", "filename", "file_type", "status")
missing_cols <- setdiff(required_cols, names(files_df))
if (length(missing_cols) > 0) stop("00a index missing columns: ", paste(missing_cols, collapse = ", "))

files_df <- files_df %>%
  mutate(
    site_base = case_when(
      file_type == "noaa_rwl"      ~ str_to_lower(str_remove(filename, "(?i)-noaa\\.rwl$")),
      file_type == "tucson_rwl"    ~ str_to_lower(str_remove(filename, "(?i)\\.rwl$")),
      file_type == "noaa_meta_txt" ~ str_to_lower(str_remove(filename, "(?i)-rwl-noaa\\.txt$")),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(site_base)) %>%
  filter(str_detect(site_base, "^[a-z]{2,4}[0-9]{3}$"))

files_ok <- files_df %>% filter(status %in% c("downloaded", "skipped"))

site_files_df <- files_ok %>%
  group_by(region, site_base) %>%
  summarise(
    rwl_noaa_file   = filename[file_type == "noaa_rwl"][1],
    rwl_tucson_file = filename[file_type == "tucson_rwl"][1],
    meta_txt_file   = filename[file_type == "noaa_meta_txt"][1],
    .groups = "drop"
  ) %>%
  arrange(region, site_base) %>%
  mutate(
    site_id  = paste(region, site_base, sep = "__"),
    site_uid = row_number(),
    meta_txt_path = ifelse(
      is.na(meta_txt_file),
      NA_character_,
      as.character(fs::path(dir_meta_txt, region, meta_txt_file))
    ),
    meta_txt_exists = !is.na(meta_txt_path) & fs::file_exists(meta_txt_path)
  )

cat("Unique sites found (status ok): ", nrow(site_files_df), "\n", sep = "")
cat("Sites with local META TXT:      ", sum(site_files_df$meta_txt_exists, na.rm = TRUE), "\n\n", sep = "")

if (sum(site_files_df$meta_txt_exists, na.rm = TRUE) == 0) {
  stop("No local meta txt found. Run updated 00a first (download *-rwl-noaa.txt).")
}

# ----------------------------- HELPERS ---------------------------------------

escape_regex <- function(x) {
  if (length(x) == 0 || is.na(x)) return(x)
  specials <- c(".", "|", "(", ")", "[", "]", "{", "}", "^", "$", "*", "+", "?")
  for (ch in specials) x <- gsub(ch, paste0("\\\\", ch), x, fixed = TRUE)
  x
}

normalize_noaa_lines <- function(lines) {
  if (is.null(lines) || length(lines) == 0) return(character())
  txt <- paste(lines, collapse = "\n")
  txt <- gsub("\r", "", txt, fixed = TRUE)
  txt <- gsub("[[:space:]]+#", "\n#", txt)
  unlist(strsplit(txt, "\n", fixed = TRUE))
}

trim_na <- function(x) {
  x <- trimws(x)
  ifelse(x == "" | is.na(x), NA_character_, x)
}

null_to_empty <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (length(x) > 1) x <- x[1]
  if (is.na(x)) return("")
  as.character(x)
}

parse_numeric <- function(x) {
  s <- null_to_empty(x)
  suppressWarnings(as.numeric(str_replace_all(s, ",", ".")))
}

parse_latlon_value <- function(x) {
  s <- null_to_empty(x)
  if (s == "") return(NA_real_)
  num <- parse_numeric(str_extract(s, "[-+]?[0-9]*\\.?[0-9]+"))
  if (is.na(num)) return(NA_real_)
  if (str_detect(s, "[SsWw]")) return(-abs(num))
  num
}

get_field_multi2 <- function(h, labels) {
  for (lab in labels) {
    pat <- paste0("^\\s*#\\s*", escape_regex(lab), "\\s*:")
    m <- grep(pat, h, value = TRUE, ignore.case = TRUE, perl = TRUE)
    if (length(m)) {
      val <- sub("^\\s*#\\s*[^:]+\\s*:\\s*", "", m[1], perl = TRUE)
      return(list(value = trim_na(val), label = lab))
    }
  }
  list(value = NA_character_, label = NA_character_)
}

infer_continent <- function(region) {
  if (str_detect(region, "^northamerica")) return("NorthAmerica")
  if (str_detect(region, "^southamerica")) return("SouthAmerica")
  if (str_detect(region, "^europe")) return("Europe")
  if (str_detect(region, "^asia")) return("Asia")
  if (str_detect(region, "^africa")) return("Africa")
  if (str_detect(region, "^australia")) return("Australia")
  "Unknown"
}

lon_wrap_mean <- function(a, b) {
  to360 <- function(x) ifelse(is.na(x), NA_real_, ifelse(x < 0, x + 360, x))
  aa <- to360(a); bb <- to360(b)
  m <- mean(c(aa, bb), na.rm = TRUE)
  if (is.na(m)) return(NA_real_)
  if (m > 180) m - 360 else m
}

LAB_LAT_N <- c("Northernmost_Latitude", "Northernmost Latitude", "North Latitude",
               "North_Latitude", "N_Latitude", "Latitude_North")
LAB_LAT_S <- c("Southernmost_Latitude", "Southernmost Latitude", "South Latitude",
               "South_Latitude", "S_Latitude", "Latitude_South")
LAB_LON_E <- c("Easternmost_Longitude", "Easternmost Longitude", "East Longitude",
               "East_Longitude", "E_Longitude", "Longitude_East")
LAB_LON_W <- c("Westernmost_Longitude", "Westernmost Longitude", "West Longitude",
               "West_Longitude", "W_Longitude", "Longitude_West")

LAB_LAT_SINGLE <- c("Latitude", "Lat", "Site_Latitude", "Site Latitude")
LAB_LON_SINGLE <- c("Longitude", "Lon", "Long", "Site_Longitude", "Site Longitude")

LAB_ELEV <- c("Elevation_m", "Elevation (m)", "Elevation", "Elev_m", "Elevation_meters", "Elev")
LAB_SITE_NAME <- c("Site_Name", "Site name")
LAB_SPECIES_NAME <- c("Species_Name", "Species name", "Species", "Tree Species",
                      "Tree_Species", "Taxon", "SpeciesName")
LAB_SPECIES_CODE <- c("Tree_Species_Code", "Tree species code", "Species_Code",
                      "Species code", "SpeciesCode")
LAB_COLLECTION_NAME <- c("Collection_Name", "Collection name")
LAB_DBH <- c("DBH_cm", "DBH", "Diameter_cm", "Diameter")
LAB_POE <- c("Pith_Offset_mm", "Pith_Offset", "Pith offset", "POE")

extract_genus <- function(species_name) {
  s <- null_to_empty(species_name)
  g <- str_extract(s, "^[A-Za-z]+")
  g <- str_to_title(g)
  ifelse(is.na(g) | g == "", NA_character_, g)
}

infer_leaf_type <- function(genus) {
  conifer_genera <- c(
    "Abies","Cedrus","Chamaecyparis","Cupressus","Juniperus","Larix",
    "Picea","Pinus","Pseudotsuga","Sequoia","Taxus","Thuja","Tsuga"
  )
  if (is.na(genus) || genus == "") return(NA_character_)
  if (genus %in% conifer_genera) return("conifer")
  "unknown"
}

# ----------------------------- PARSER ----------------------------------------

process_site_metadata_local <- function(region, site_base, meta_txt_path) {
  
  base_row <- tibble(
    region = region,
    continent = infer_continent(region),
    site_base = site_base,
    site_id = paste(region, site_base, sep = "__"),
    meta_txt_path = meta_txt_path,
    meta_status = NA_character_,
    meta_note = NA_character_,
    
    site_code = toupper(site_base),
    site_name = NA_character_,
    
    lat = NA_real_,
    lon = NA_real_,
    lat_source = NA_character_,
    lon_source = NA_character_,
    
    elevation_m = NA_real_,
    first_year = NA_integer_,
    last_year = NA_integer_,
    
    species_name = NA_character_,
    species_code = NA_character_,
    genus = NA_character_,
    leaf_type = NA_character_,
    
    dbh_cm = NA_real_,
    poe_mm = NA_real_
  )
  
  if (is.na(meta_txt_path) || !fs::file_exists(meta_txt_path)) {
    return(base_row %>% mutate(meta_status = "FAIL", meta_note = "Local meta txt missing"))
  }
  
  tryCatch({
    
    raw_lines <- tryCatch(readLines(meta_txt_path, warn = FALSE, encoding = "UTF-8"),
                          error = function(e) NULL)
    
    if (is.null(raw_lines) || length(raw_lines) == 0) {
      return(base_row %>% mutate(meta_status = "FAIL", meta_note = "Cannot read local meta txt"))
    }
    
    lines <- normalize_noaa_lines(raw_lines)
    h <- lines[grepl("^\\s*#", lines)]
    if (length(h) == 0) {
      return(base_row %>% mutate(meta_status = "FAIL", meta_note = "No header lines found"))
    }
    
    # ---- LAT: bounds + enforce hemisphere by FIELD meaning ----
    north <- get_field_multi2(h, LAB_LAT_N)
    south <- get_field_multi2(h, LAB_LAT_S)
    lat_single <- get_field_multi2(h, LAB_LAT_SINGLE)
    
    north_v <- parse_latlon_value(north$value)
    south_v <- parse_latlon_value(south$value)
    single_lat_v <- parse_latlon_value(lat_single$value)
    
    # REMOVED: forced abs() lines for lat
    
    lat_bounds <- if (all(is.na(c(north_v, south_v)))) NA_real_ else mean(c(north_v, south_v), na.rm = TRUE)
    lat_val <- dplyr::coalesce(lat_bounds, single_lat_v)
    
    lat_src <- if (!is.na(lat_bounds)) {
      paste0("bounds:", paste(na.omit(c(north$label, south$label)), collapse = "+"))
    } else if (!is.na(single_lat_v)) {
      paste0("single:", lat_single$label)
    } else {
      NA_character_
    }
    
    # ---- LON: bounds + enforce hemisphere by FIELD meaning ----
    east <- get_field_multi2(h, LAB_LON_E)
    west <- get_field_multi2(h, LAB_LON_W)
    lon_single <- get_field_multi2(h, LAB_LON_SINGLE)
    
    east_v <- parse_latlon_value(east$value)
    west_v <- parse_latlon_value(west$value)
    single_lon_v <- parse_latlon_value(lon_single$value)
    
    # REMOVED: forced abs() lines for lon
    
    lon_bounds <- NA_real_
    if (!all(is.na(c(east_v, west_v)))) {
      if (!is.na(east_v) && !is.na(west_v) && abs(east_v - west_v) > 180) {
        lon_bounds <- lon_wrap_mean(east_v, west_v)
      } else {
        lon_bounds <- mean(c(east_v, west_v), na.rm = TRUE)
      }
    }
    
    lon_val <- dplyr::coalesce(lon_bounds, single_lon_v)
    
    # Fix missing 'W' for Americas: NOAA headers sometimes store 148 (meaning 148W)
    if (stringr::str_detect(region, "^(northamerica|southamerica)") &&
        !is.na(lon_val) && lon_val > 0) {
      lon_val <- -abs(lon_val)
    }
    
    lon_src <- if (!is.na(lon_bounds)) {
      paste0("bounds:", paste(na.omit(c(east$label, west$label)), collapse = "+"))
    } else if (!is.na(single_lon_v)) {
      paste0("single:", lon_single$label)
    } else {
      NA_character_
    }
    
    elev_val <- parse_numeric(get_field_multi2(h, LAB_ELEV)$value)
    if (!is.na(elev_val) && elev_val <= -900) elev_val <- NA_real_
    
    fyear_val <- suppressWarnings(as.integer(get_field_multi2(h, c("First_Year", "First year"))$value))
    lyear_val <- suppressWarnings(as.integer(get_field_multi2(h, c("Last_Year", "Last year"))$value))
    
    site_name_val    <- get_field_multi2(h, LAB_SITE_NAME)$value
    species_name_val <- get_field_multi2(h, LAB_SPECIES_NAME)$value
    species_code_val <- get_field_multi2(h, LAB_SPECIES_CODE)$value
    
    dbh_raw <- get_field_multi2(h, LAB_DBH)$value
    poe_raw <- get_field_multi2(h, LAB_POE)$value
    
    dbh_cm_val <- parse_numeric(str_extract(null_to_empty(dbh_raw), "[-+]?[0-9]*\\.?[0-9]+"))
    poe_mm_val <- parse_numeric(str_extract(null_to_empty(poe_raw), "[-+]?[0-9]*\\.?[0-9]+"))
    
    site_code_extracted <- toupper(get_field_multi2(h, LAB_COLLECTION_NAME)$value)
    if (is.na(site_code_extracted) || site_code_extracted == "") site_code_extracted <- toupper(site_base)
    
    genus_val <- extract_genus(species_name_val)
    leaf_type_val <- infer_leaf_type(genus_val)
    
    base_row %>%
      mutate(
        meta_status  = "OK",
        meta_note    = NA_character_,
        site_code    = site_code_extracted,
        site_name    = site_name_val,
        lat          = lat_val,
        lon          = lon_val,
        lat_source   = lat_src,
        lon_source   = lon_src,
        elevation_m  = elev_val,
        first_year   = fyear_val,
        last_year    = lyear_val,
        species_name = species_name_val,
        species_code = species_code_val,
        genus        = genus_val,
        leaf_type    = leaf_type_val,
        dbh_cm       = dbh_cm_val,
        poe_mm       = poe_mm_val
      )
    
  }, error = function(e) {
    base_row %>% mutate(meta_status = "FAIL", meta_note = paste0("ERROR: ", e$message))
  })
}

run_batch <- function(batch_df, tag) {
  cat("\n--- Parsing LOCAL meta txt:", tag, "---\n")
  
  input <- batch_df %>% select(region, site_base, meta_txt_path)
  
  if (isTRUE(CFG$PARALLEL)) {
    res <- furrr::future_pmap_dfr(
      .l = input,
      .f = function(region, site_base, meta_txt_path) {
        process_site_metadata_local(region, site_base, meta_txt_path)
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    res <- purrr::pmap_dfr(
      .l = input,
      .f = function(region, site_base, meta_txt_path) {
        process_site_metadata_local(region, site_base, meta_txt_path)
      }
    )
  }
  
  out <- res %>%
    left_join(
      site_files_df %>% select(region, site_base, site_uid, rwl_noaa_file, rwl_tucson_file),
      by = c("region", "site_base")
    ) %>%
    mutate(
      lat_in_range = !is.na(lat) & lat >= -90 & lat <= 90,
      lon_in_range = !is.na(lon) & lon >= -180 & lon <= 180,
      has_latlon   = lat_in_range & lon_in_range,
      has_species  = !is.na(species_name),
      keep_lat50   = has_latlon & lat >= CFG$LAT_MIN_KEEP
    )
  
  # Fix transcontinental Russia (no downloads): split by lon ~ Ural (60E)
  out <- out %>%
    mutate(
      continent_raw = continent,
      continent = case_when(
        has_latlon & str_detect(site_base, "^rus") & lon < 60  ~ "Europe",
        has_latlon & str_detect(site_base, "^rus") & lon >= 60 ~ "Asia",
        TRUE ~ continent
      )
    )
  
  qc <- list(
    meta_ok_rate = mean(out$meta_status == "OK", na.rm = TRUE),
    latlon_rate  = mean(out$has_latlon, na.rm = TRUE),
    species_rate = mean(out$has_species, na.rm = TRUE)
  )
  
  cat("\n--- QC:", tag, "---\n")
  cat("Rows:           ", nrow(out), "\n", sep = "")
  cat(sprintf("Meta OK rate:   %.3f\n", qc$meta_ok_rate))
  cat(sprintf("Lat/Lon rate:   %.3f\n", qc$latlon_rate))
  cat(sprintf("Species rate:   %.3f\n", qc$species_rate))
  
  list(df = out, qc = qc)
}

qc_passes <- function(qc) {
  isTRUE(qc$meta_ok_rate >= CFG$MIN_META_OK_RATE &&
           qc$latlon_rate  >= CFG$MIN_LATLON_RATE  &&
           qc$species_rate >= CFG$MIN_SPECIES_RATE)
}

# ----------------------------- RUN -------------------------------------------

if (identical(CFG$MODE, "test_then_all")) {
  
  test_df <- site_files_df %>% slice_sample(n = min(CFG$N_TEST, nrow(site_files_df)))
  test_run <- run_batch(test_df, tag = paste0("TEST_", nrow(test_df)))
  
  test_csv <- fs::path(dir_dbg, paste0("01a_ITRDB_site_metadata_TEST_", nrow(test_df), ".csv"))
  readr::write_csv(test_run$df, test_csv)
  cat("Saved test CSV: ", test_csv, "\n", sep = "")
  
  if (!qc_passes(test_run$qc)) {
    cat("\n!!! QC FAILED. Not running ALL.\n")
    fail_examples <- test_run$df %>%
      filter(meta_status != "OK" | !has_latlon | is.na(species_name)) %>%
      select(site_id, meta_txt_path, meta_status, meta_note, lat, lon, species_name, lat_source, lon_source) %>%
      slice_head(n = 30)
    out_fail <- fs::path(dir_dbg, "01a_TEST_fail_examples.csv")
    readr::write_csv(fail_examples, out_fail)
    cat("Saved fail examples: ", out_fail, "\n", sep = "")
    quit(save = "no", status = 0)
  }
  
  if (!isTRUE(CFG$AUTO_RUN_ALL_IF_PASS)) {
    cat("\nQC PASSED but AUTO_RUN_ALL_IF_PASS=FALSE -> stop after test.\n")
    quit(save = "no", status = 0)
  }
}

all_run <- run_batch(site_files_df, tag = "ALL")
meta_df <- all_run$df

# ----------------------------- EXPORT ----------------------------------------

out_full  <- fs::path(dir_tab, "01a_ITRDB_site_metadata.csv")
out_lat50 <- fs::path(dir_tab, "01a_ITRDB_site_metadata_lat50.csv")

meta_fail <- meta_df %>%
  filter(meta_status != "OK") %>%
  select(site_id, site_uid, continent, continent_raw, region, site_base, meta_txt_path, meta_status, meta_note)

missing_latlon <- meta_df %>%
  filter(meta_status == "OK" & !has_latlon) %>%
  select(site_id, site_uid, continent, continent_raw, region, site_base, meta_txt_path, lat, lon, lat_source, lon_source)

qc_summary <- tibble(
  n_sites      = nrow(meta_df),
  n_meta_ok    = sum(meta_df$meta_status == "OK", na.rm = TRUE),
  n_has_latlon = sum(meta_df$has_latlon, na.rm = TRUE),
  n_keep_lat50 = sum(meta_df$keep_lat50, na.rm = TRUE),
  meta_ok_rate = mean(meta_df$meta_status == "OK", na.rm = TRUE),
  latlon_rate  = mean(meta_df$has_latlon, na.rm = TRUE),
  species_rate = mean(meta_df$has_species, na.rm = TRUE)
)

meta_lat50 <- meta_df %>%
  filter(meta_status == "OK", has_latlon, lat >= CFG$LAT_MIN_KEEP)

readr::write_csv(meta_df, out_full)
readr::write_csv(meta_lat50, out_lat50)
readr::write_csv(meta_fail, fs::path(dir_tab, "01a_meta_failures.csv"))
readr::write_csv(missing_latlon, fs::path(dir_tab, "01a_missing_latlon.csv"))
readr::write_csv(qc_summary, fs::path(dir_tab, "01a_qc_summary.csv"))

# ----------------------------- FINAL PRINTS ----------------------------------

cat("\n=== 01a DONE (OFFLINE) ===\n")
cat("Saved FULL:    ", out_full, "\n", sep = "")
cat("Saved LAT>=50: ", out_lat50, "\n", sep = "")

cat("\nQC summary:\n")
print(qc_summary)

cat("\nCounts by continent (ALL):\n")
print(table(meta_df$continent, useNA = "ifany"))

cat("\nCounts by continent_raw (ALL):\n")
print(table(meta_df$continent_raw, useNA = "ifany"))

cat("\nCounts by continent (LAT>=50 + OK + valid lat/lon):\n")
print(table(meta_lat50$continent, useNA = "ifany"))

cat("\nMeta status table:\n")
print(table(meta_df$meta_status, useNA = "ifany"))

cat("\nLat source (top 10):\n")
print(head(sort(table(meta_df$lat_source, useNA = "ifany"), decreasing = TRUE), 10))

cat("\nLon source (top 10):\n")
print(head(sort(table(meta_df$lon_source, useNA = "ifany"), decreasing = TRUE), 10))

cat("\nOutput columns (FULL):\n")
cat(paste(names(meta_df), collapse = ", "), "\n")

cat("\nPreview FULL (first 5 rows):\n")
print(meta_df %>% select(site_id, site_uid, continent, continent_raw, lat, lon, species_name, meta_status, meta_note, keep_lat50) %>% head(5))

cat("\nPreview LAT>=50 (first 5 rows):\n")
print(meta_lat50 %>% select(site_id, site_uid, continent, continent_raw, lat, lon, species_name) %>% head(5))

cat("\nMeta failures saved:     ", fs::path(dir_tab, "01a_meta_failures.csv"), "\n", sep = "")
cat("Missing lat/lon saved:  ", fs::path(dir_tab, "01a_missing_latlon.csv"), "\n", sep = "")
cat("QC summary saved:       ", fs::path(dir_tab, "01a_qc_summary.csv"), "\n", sep = "")
