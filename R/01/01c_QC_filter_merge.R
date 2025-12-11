###############################################################################
# 01c_QC_filter_merge.R
#
# Description:
#   1. Merge BAI statistics (01b) with Site Metadata (01a).
#   2. Apply strict filters:
#      - Latitude >= 50 (Northern Hemisphere High Latitudes)
#      - Series count >= 10 (Robustness)
#      - Time span >= 50 years (Dhyani et al. 2025 requirement)
#   3. Load BAI time series ONLY for the valid sites.
#   4. Save cleaned long-format dataset.
#
# Input:
#   - output/01a_analyse_ITRDB/tables/01a_ITRDB_site_metadata.csv (Lat/Lon)
#   - output/01b_compute_BAI/tables/01b_BAI_summary_stats.csv (Years/Counts)
#   - output/01b_compute_BAI/data/bai/ (Actual Data)
#
# Output:
#   - output/01c_QC_filter_merge/data/01c_BAI_long_filtered.rds
#   - output/01c_QC_filter_merge/tables/01c_QC_pass_list.csv
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  purrr,
  tidyr,
  stringr
)

# --- 2) Directories & Paths --------------------------------------------------

step_01a <- "01a_analyse_ITRDB"
step_01b <- "01b_compute_BAI"
current  <- "01c_QC_filter_merge"

# Inputs
path_meta  <- here::here("output", step_01a, "tables", "01a_ITRDB_site_metadata.csv")
path_stats <- here::here("output", step_01b, "tables", "01b_BAI_summary_stats.csv")
dir_bai    <- here::here("output", step_01b, "data", "bai")

# Outputs
base_out_dir <- here::here("output", current)
dir_data     <- path(base_out_dir, "data")
dir_tables   <- path(base_out_dir, "tables")

dir_create(c(dir_data, dir_tables))

cat("=== 01c QC, Filter & Merge ===\n\n")

# --- 3) Load Metadata & Stats ------------------------------------------------

if (!file_exists(path_meta) || !file_exists(path_stats)) {
  stop("Input files missing. Check 01a and 01b tables.")
}

cat("Loading metadata and stats...\n")
df_meta  <- read_csv(path_meta, show_col_types = FALSE)
df_stats <- read_csv(path_stats, show_col_types = FALSE)

# Join them (Stats has the success info, Meta has Lat/Lon)
# We filter df_stats to only successful BAI calculations first
df_merged <- df_stats %>%
  filter(status == "OK") %>%
  left_join(
    df_meta %>% select(site_code, region, lat, lon, species_name, elevation_m),
    by = c("site_code", "region")
  ) %>%
  mutate(
    time_span = last_year - first_year + 1
  )

cat(sprintf("Total sites with successful BAI: %d\n", nrow(df_merged)))

# --- 4) Apply Filters --------------------------------------------------------

# Define Thresholds
THRESH_LAT    <- 50   # >= 50 Degrees North
THRESH_SERIES <- 10   # >= 10 Trees
THRESH_YEARS  <- 50   # >= 50 Years (Paper Dhyani et al. 2025)

cat("\n--- APPLYING FILTERS ---\n")

# 1. Latitude Filter
pass_lat <- df_merged %>% filter(lat >= THRESH_LAT)
cat(sprintf("1. Latitude >= %d°N:     %d sites remaining (Dropped: %d)\n", 
            THRESH_LAT, nrow(pass_lat), nrow(df_merged) - nrow(pass_lat)))

# 2. Series Count Filter
pass_series <- pass_lat %>% filter(n_series >= THRESH_SERIES)
cat(sprintf("2. Series (Trees) >= %d:   %d sites remaining (Dropped: %d)\n", 
            THRESH_SERIES, nrow(pass_series), nrow(pass_lat) - nrow(pass_series)))

# 3. Time Span Filter
final_sites <- pass_series %>% filter(time_span >= THRESH_YEARS)
cat(sprintf("3. Time Span >= %d yrs:   %d sites remaining (Dropped: %d)\n", 
            THRESH_YEARS, nrow(final_sites), nrow(pass_series) - nrow(final_sites)))

cat(sprintf("\nFinal Selection: %d sites ready for analysis.\n", nrow(final_sites)))

if (nrow(final_sites) == 0) {
  stop("No sites left after filtering! Check your thresholds.")
}

# Save the list of passing sites
write_csv(final_sites, path(dir_tables, "01c_QC_pass_list.csv"))

# --- 5) Load Actual BAI Data for Valid Sites ---------------------------------

cat("\n--- Loading BAI Data for selected sites ---\n")

# Helper to construct path and load
load_bai_for_site <- function(site_code, region) {
  # Path construction: output/01b.../data/bai/{region}/{site_code}_BAI.csv
  fpath <- path(dir_bai, region, paste0(site_code, "_BAI.csv"))
  
  if (!file_exists(fpath)) return(NULL)
  
  # Read and pivot to long format
  read_csv(fpath, show_col_types = FALSE) %>%
    pivot_longer(
      cols = -year,
      names_to = "tree_id",
      values_to = "bai"
    ) %>%
    filter(!is.na(bai)) %>%
    mutate(site_code = site_code) %>%
    select(site_code, year, tree_id, bai)
}

# Apply over the filtered list
# Using map2_dfr to iterate over site_code and region
bai_long_clean <- purrr::map2_dfr(
  .x = final_sites$site_code,
  .y = final_sites$region,
  .f = load_bai_for_site
)

# Join Metadata back into the long format (optional, but useful for 02a stats)
# Doing a lighter join (just basics)
bai_long_clean <- bai_long_clean %>%
  left_join(
    final_sites %>% select(site_code, region, species_code, lat, lon),
    by = "site_code"
  )

cat("Data loaded. Total rows (observations): ", nrow(bai_long_clean), "\n")

# --- 6) Save Final Dataset ---------------------------------------------------

out_rds <- path(dir_data, "01c_BAI_long_filtered.rds")
saveRDS(bai_long_clean, out_rds)

cat("\n=== OUTPUT SAVED ===\n")
cat("Clean Data (RDS): ", out_rds, "\n")
cat("Site List (CSV):  ", path(dir_tables, "01c_QC_pass_list.csv"), "\n")
cat("\nProceed to 01d (Combine with Climate Data).\n")