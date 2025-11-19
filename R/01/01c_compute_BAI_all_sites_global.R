###############################################################################
# 01c_compute_BAI_all_sites_global.R
#
# Description:
#   Compute Basal Area Increment (BAI) for all sites with fallback to NOAA RWL
#   if Tucson RWL reading fails.
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  furrr,
  dplR,
  tidyr,
  stringr
)

num_cores <- parallel::detectCores() - 1
plan(multisession, workers = num_cores)

# --- 2) Directories & Paths --------------------------------------------------

step_01a <- "01a_download_ITRDB_RWL_global"
step_01b <- "01b_build_BAI_base_global"
current_name <- "01c_compute_BAI_all_sites_global"

input_rwl_base <- here::here("output", step_01a, "data", "rwl")
input_meta_csv <- here::here("output", step_01b, "tables", "01b_ITRDB_site_metadata.csv")

base_out_dir <- here::here("output", current_name)
dir_data     <- path(base_out_dir, "data")
dir_bai      <- path(dir_data, "bai")
dir_tables   <- path(base_out_dir, "tables")
dir_figures  <- path(base_out_dir, "figures")

dir_create(c(dir_data, dir_bai, dir_tables, dir_figures))

cat("=== 01c Compute BAI All Sites Global ===\n\n")
cat("RWL files from: ", input_rwl_base, "\n")
cat("Metadata from:  ", input_meta_csv, "\n")
cat("Output:         ", base_out_dir, "\n\n")

# --- 3) Load & Prepare Data --------------------------------------------------

if (!file_exists(input_meta_csv)) {
  stop("Metadata file not found. Run 01b first.")
}

meta_df <- read_csv(input_meta_csv, show_col_types = FALSE)

process_df <- meta_df %>%
  filter(has_rwl_tucson_local == TRUE | has_rwl_noaa_local == TRUE) %>%
  mutate(
    full_rwl_tucson_path = if_else(
      has_rwl_tucson_local,
      as.character(path(input_rwl_base, region, rwl_tucson_file)),
      NA_character_
    ),
    full_rwl_noaa_path = if_else(
      has_rwl_noaa_local,
      as.character(path(input_rwl_base, region, rwl_noaa_file)),
      NA_character_
    )
  )

if (nrow(process_df) == 0) stop("No local RWL files found in metadata.")

cat(sprintf("Identified %d sites to process.\n\n", nrow(process_df)))

unique_regions <- unique(process_df$region)
for (reg in unique_regions) {
  dir_create(path(dir_bai, reg))
}

# --- 4) BAI Computation Function with Fallback -------------------------------

compute_single_site_bai <- function(site_code, region, rwl_tucson_path, rwl_noaa_path, 
                                    site_name, species_code) {
  
  out_file <- path(dir_bai, region, paste0(site_code, "_BAI.csv"))
  
  make_row <- function(status, method = NA_character_, n_series = NA_integer_,
                       first_year = NA_integer_, last_year = NA_integer_,
                       bai_mean_cm2 = NA_real_, bai_max_cm2 = NA_real_,
                       bai_median_cm2 = NA_real_, error_detail = NA_character_) {
    tibble(
      region = region, site_code = site_code, site_name = site_name,
      species_code = species_code, status = status, method = method,
      n_series = n_series, first_year = first_year, last_year = last_year,
      bai_mean_cm2 = bai_mean_cm2, bai_max_cm2 = bai_max_cm2,
      bai_median_cm2 = bai_median_cm2, error_detail = error_detail
    )
  }
  
  # Try Tucson first, then NOAA
  paths_to_try <- c(rwl_tucson_path, rwl_noaa_path)
  paths_to_try <- paths_to_try[!is.na(paths_to_try)]
  
  if (length(paths_to_try) == 0) {
    return(make_row("Missing RWL file", error_detail = "No RWL files available"))
  }
  
  rwl <- NULL
  used_path <- NA_character_
  
  for (p in paths_to_try) {
    if (!file_exists(p)) next
    
    rwl <- tryCatch(
      suppressWarnings(dplR::read.rwl(p)),
      error = function(e) NULL
    )
    
    if (!is.null(rwl) && nrow(rwl) > 0 && ncol(rwl) > 0) {
      used_path <- p
      break
    }
  }
  
  if (is.null(rwl)) {
    return(make_row("Read Error", error_detail = "All RWL files failed to read"))
  }
  
  file_type <- if_else(str_detect(used_path, "noaa"), "NOAA", "Tucson")
  
  # Clean data
  rwl[rwl <= 0] <- NA
  rwl[rwl > 50] <- NA
  
  valid_counts <- colSums(!is.na(rwl))
  rwl <- rwl[, valid_counts >= 10, drop = FALSE]
  
  if (ncol(rwl) == 0) {
    return(make_row("Insufficient Data", method = file_type,
                    error_detail = "No series with >=10 valid observations"))
  }
  
  # Compute BAI
  bai_mm2 <- tryCatch(
    suppressWarnings(dplR::bai.in(rwl)),
    error = function(e) NULL
  )
  
  if (is.null(bai_mm2) || all(is.na(as.matrix(bai_mm2)))) {
    return(make_row("BAI Calc Failed", method = file_type, n_series = ncol(rwl),
                    error_detail = "BAI computation failed"))
  }
  
  # Convert to cm² and filter
  bai_cm2 <- bai_mm2 / 100
  bai_cm2[bai_cm2 > 100] <- NA
  bai_cm2[bai_cm2 <= 0] <- NA
  
  years <- suppressWarnings(as.integer(rownames(bai_cm2)))
  if (any(is.na(years))) {
    return(make_row("Invalid Years", method = file_type, n_series = ncol(bai_cm2),
                    error_detail = "Year parsing failed"))
  }
  
  # Save
  bai_out <- bai_cm2 %>%
    mutate(year = years) %>%
    select(year, everything())
  
  write_csv(bai_out, out_file)
  
  # Summary stats
  vals <- as.vector(as.matrix(bai_cm2))
  vals <- vals[!is.na(vals) & vals > 0]
  
  if (length(vals) == 0) {
    return(make_row("No Valid Data", method = file_type, n_series = ncol(bai_cm2),
                    first_year = min(years), last_year = max(years),
                    error_detail = "All BAI values invalid after filtering"))
  }
  
  make_row(
    status = "OK", method = file_type, n_series = ncol(bai_cm2),
    first_year = min(years), last_year = max(years),
    bai_mean_cm2 = mean(vals), bai_max_cm2 = max(vals),
    bai_median_cm2 = median(vals)
  )
}

# --- 5) Execute Parallel Processing ------------------------------------------

cat("=== STARTING BAI COMPUTATION ===\n")

bai_stats_df <- process_df %>%
  select(site_code, region, full_rwl_tucson_path, full_rwl_noaa_path, 
         site_name, species_code) %>%
  future_pmap_dfr(
    .f = function(site_code, region, full_rwl_tucson_path, full_rwl_noaa_path,
                  site_name, species_code) {
      compute_single_site_bai(site_code, region, full_rwl_tucson_path, 
                              full_rwl_noaa_path, site_name, species_code)
    },
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

# --- 6) Save Summary & Report ------------------------------------------------

cat("\n=== RESULTS SUMMARY ===\n")

successes <- bai_stats_df %>% filter(status == "OK")
failures  <- bai_stats_df %>% filter(status != "OK")

cat("Successfully computed: ", nrow(successes), " sites\n")
cat("Failed/Skipped:        ", nrow(failures), " sites\n")
cat("Success rate:          ", 
    sprintf("%.1f%%", 100 * nrow(successes) / nrow(bai_stats_df)), "\n\n")

if (nrow(failures) > 0) {
  cat("Failure breakdown:\n")
  print(table(failures$status))
  cat("\n")
}

if (nrow(successes) > 0) {
  cat("File types used:\n")
  print(table(successes$method))
  cat("\n")
}

summary_path <- path(dir_tables, "01c_BAI_summary_stats.csv")
write_csv(bai_stats_df, summary_path)
write_csv(successes, path(dir_tables, "01c_BAI_successes.csv"))
write_csv(failures, path(dir_tables, "01c_BAI_failures.csv"))

cat("=== OUTPUT FILES ===\n")
cat("  Summary stats: ", summary_path, "\n")
cat("  Successes:     ", path(dir_tables, "01c_BAI_successes.csv"), "\n")
cat("  Failures:      ", path(dir_tables, "01c_BAI_failures.csv"), "\n")
cat("  BAI data:      ", dir_bai, "\n")

cat("\n=== 01c Compute BAI All Sites COMPLETE ===\n")