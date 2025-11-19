###############################################################################
# 01e_extract_BAI_boreal_global.R
#
# Description:
#   Extract BAI data for boreal forest sites (>= 50°N) that passed QC.
#   Creates a combined dataset with site metadata and BAI time series.
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  stringr,
  tidyr,
  purrr,
  furrr
)

num_cores <- parallel::detectCores() - 1
plan(multisession, workers = num_cores)

# --- 2) Directories & Paths --------------------------------------------------

step_01c <- "01c_compute_BAI_all_sites_global"
step_01d <- "01d_QC_BAI_global"

qc_results_path <- here::here("output", step_01d, "tables", "01d_QC_results.csv")
bai_root <- here::here("output", step_01c, "data", "bai")

base_out_dir <- here::here("output", "01e_extract_BAI_boreal_global")
dir_tables   <- path(base_out_dir, "tables")
dir_data     <- path(base_out_dir, "data")

dir_create(c(dir_tables, dir_data))

cat("=== 01e Extract BAI Boreal Global ===\n\n")
cat("QC results: ", qc_results_path, "\n")
cat("BAI dir:    ", bai_root, "\n")
cat("Output:     ", base_out_dir, "\n\n")

# --- 3) Load QC Results & Filter ---------------------------------------------

if (!file_exists(qc_results_path)) {
  stop("QC results not found. Run 01d_QC_BAI_global first.")
}

qc_data <- read_csv(qc_results_path, show_col_types = FALSE)
cat("Total sites in QC results: ", nrow(qc_data), "\n")

passed_sites <- qc_data %>%
  filter(qc_status == "PASS") %>%
  mutate(
    path_bai = file.path(bai_root, region, paste0(site_code, "_BAI.csv"))
  )

cat("Sites that passed QC:      ", nrow(passed_sites), "\n\n")

if (nrow(passed_sites) == 0) {
  stop("No sites passed QC. Check 01d results.")
}

# --- 4) Extract BAI Data -----------------------------------------------------

cat("=== EXTRACTING BAI DATA ===\n")

extract_bai_with_metadata <- function(site_code, path_bai, site_name, lat, lon, 
                                      elevation_m, genus, scientific_name, species_code, 
                                      region) {
  
  if (!file_exists(path_bai)) {
    return(NULL)
  }
  
  bai_data <- tryCatch(
    read_csv(path_bai, show_col_types = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(bai_data) || nrow(bai_data) == 0) {
    return(NULL)
  }
  
  bai_long <- bai_data %>%
    pivot_longer(
      cols = -year,
      names_to = "series_id",
      values_to = "bai_cm2"
    ) %>%
    filter(!is.na(bai_cm2)) %>%
    mutate(
      site_code       = site_code,
      site_name       = site_name,
      lat             = lat,
      lon             = lon,
      elevation_m     = elevation_m,
      genus           = genus,
      scientific_name = scientific_name,
      species_code    = species_code,
      region          = region
    ) %>%
    select(
      site_code, site_name, series_id, year, bai_cm2,
      lat, lon, elevation_m, genus, scientific_name, species_code,
      region
    )
  
  return(bai_long)
}

bai_combined <- passed_sites %>%
  select(
    site_code, path_bai, site_name, lat, lon, elevation_m,
    genus, scientific_name, species_code, region
  ) %>%
  future_pmap_dfr(
    extract_bai_with_metadata,
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )

cat("\nExtracted BAI records:\n")
cat("  Total observations: ", nrow(bai_combined), "\n")
cat("  Unique sites:       ", n_distinct(bai_combined$site_code), "\n")
cat("  Unique series:      ", n_distinct(bai_combined$series_id), "\n")
cat("  Year range:         ", min(bai_combined$year, na.rm = TRUE), " - ",
    max(bai_combined$year, na.rm = TRUE), "\n\n")

# --- 5) Save Combined Dataset ------------------------------------------------

combined_path <- path(dir_data, "01e_BAI_boreal_combined.csv")
write_csv(bai_combined, combined_path)

cat("Combined BAI dataset saved to:\n")
cat("  ", combined_path, "\n\n")

# --- 6) Summary Statistics ---------------------------------------------------

cat("=== SUMMARY STATISTICS ===\n")

site_summary <- bai_combined %>%
  group_by(site_code, site_name, lat, lon, elevation_m, genus, 
           scientific_name, species_code, region) %>%
  summarise(
    n_series      = n_distinct(series_id),
    n_years       = n_distinct(year),
    year_min      = min(year, na.rm = TRUE),
    year_max      = max(year, na.rm = TRUE),
    mean_bai      = mean(bai_cm2, na.rm = TRUE),
    median_bai    = median(bai_cm2, na.rm = TRUE),
    sd_bai        = sd(bai_cm2, na.rm = TRUE),
    total_obs     = n(),
    .groups = "drop"
  )

write_csv(site_summary, path(dir_tables, "01e_site_summary.csv"))

cat("\nSite-level summary:\n")
print(site_summary %>% head(10), width = 120)

genus_summary <- bai_combined %>%
  group_by(genus) %>%
  summarise(
    n_sites       = n_distinct(site_code),
    n_series      = n_distinct(series_id),
    n_observations = n(),
    mean_bai      = mean(bai_cm2, na.rm = TRUE),
    median_bai    = median(bai_cm2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sites))

write_csv(genus_summary, path(dir_tables, "01e_genus_summary.csv"))

cat("\nGenus-level summary:\n")
print(genus_summary, n = 15)

region_summary <- bai_combined %>%
  group_by(region) %>%
  summarise(
    n_sites       = n_distinct(site_code),
    n_series      = n_distinct(series_id),
    n_observations = n(),
    mean_bai      = mean(bai_cm2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sites))

write_csv(region_summary, path(dir_tables, "01e_region_summary.csv"))

cat("\nRegion-level summary:\n")
print(region_summary, n = 20)

# Temporal coverage
temporal_summary <- bai_combined %>%
  group_by(year) %>%
  summarise(
    n_sites       = n_distinct(site_code),
    n_series      = n_distinct(series_id),
    n_observations = n(),
    mean_bai      = mean(bai_cm2, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(temporal_summary, path(dir_tables, "01e_temporal_coverage.csv"))

cat("\nTemporal coverage saved.\n")

# --- 7) Data Quality Checks --------------------------------------------------

cat("\n=== DATA QUALITY CHECKS ===\n")

na_summary <- bai_combined %>%
  summarise(
    total_rows = n(),
    na_bai = sum(is.na(bai_cm2)),
    na_year = sum(is.na(year)),
    na_site = sum(is.na(site_code)),
    pct_na_bai = round(100 * mean(is.na(bai_cm2)), 2)
  )

print(na_summary)

range_summary <- bai_combined %>%
  summarise(
    min_bai = min(bai_cm2, na.rm = TRUE),
    max_bai = max(bai_cm2, na.rm = TRUE),
    q01_bai = quantile(bai_cm2, 0.01, na.rm = TRUE),
    q99_bai = quantile(bai_cm2, 0.99, na.rm = TRUE)
  )

cat("\nBAI range:\n")
print(range_summary)

outliers <- bai_combined %>%
  filter(bai_cm2 > 50 | bai_cm2 < 0.01) %>%
  select(site_code, series_id, year, bai_cm2, genus)

if (nrow(outliers) > 0) {
  cat(sprintf("\nWarning: %d potential outliers (BAI > 50 or < 0.01 cm²)\n", 
              nrow(outliers)))
  write_csv(outliers, path(dir_tables, "01e_potential_outliers.csv"))
} else {
  cat("\nNo extreme outliers detected.\n")
}

# --- 8) Output Summary -------------------------------------------------------

cat("\n=== OUTPUT FILES ===\n")
cat("  Combined BAI data:    ", combined_path, "\n")
cat("  Site summary:         ", path(dir_tables, "01e_site_summary.csv"), "\n")
cat("  Genus summary:        ", path(dir_tables, "01e_genus_summary.csv"), "\n")
cat("  Region summary:       ", path(dir_tables, "01e_region_summary.csv"), "\n")
cat("  Temporal coverage:    ", path(dir_tables, "01e_temporal_coverage.csv"), "\n")
if (nrow(outliers) > 0) {
  cat("  Potential outliers:   ", path(dir_tables, "01e_potential_outliers.csv"), "\n")
}

cat("\n=== 01e Extract BAI Boreal Global COMPLETE ===\n")
cat(sprintf("Successfully extracted BAI data for %d sites\n", 
            n_distinct(bai_combined$site_code)))
cat(sprintf("Total records: %d observations across %d series\n",
            nrow(bai_combined), n_distinct(bai_combined$series_id)))