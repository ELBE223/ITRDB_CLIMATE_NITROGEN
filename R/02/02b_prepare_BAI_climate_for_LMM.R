###############################################################################
# 02b_prepare_BAI_climate_for_LMM.R
#
# Description:
#   Combine boreal BAI data with climate and nitrogen deposition data.
#   Extracts site-level annual values for:
#   - CRU TS v4.09 (temperature, precipitation)
#   - ISIMIP3a N-deposition (ndep-nhx, ndep-noy)
#   Restricted to common time window 1901-2010.
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  tidyr,
  tibble,
  arrow,
  terra,
  stringr
)

# --- 2) Time Window & Paths --------------------------------------------------

year_min <- 1901L
year_max <- 2010L

step_01e <- "01e_extract_BAI_boreal_global"
step_02a <- "02a_download_climate_and_ndep"
current_name <- "02b_prepare_BAI_climate_for_LMM"

bai_data_path <- here::here("output", step_01e, "data", "01e_BAI_boreal_combined.csv")
bai_site_path <- here::here("output", step_01e, "tables", "01e_site_summary.csv")

clim_base <- here::here("output", step_02a, "data", "climate")
cru_dir   <- path(clim_base, "CRU_TS_v4.09")
ndep_root <- path(clim_base, "ISIMIP3a", "ndep", "histsoc")

base_out_dir <- here::here("output", current_name)
dir_data     <- path(base_out_dir, "data")
dir_tables   <- path(base_out_dir, "tables")
dir_figures  <- path(base_out_dir, "figures")

dir_create(c(dir_data, dir_tables, dir_figures))

cat("=== 02b Prepare BAI Climate for LMM ===\n\n")
cat("BAI data:     ", bai_data_path, "\n")
cat("Climate data: ", clim_base, "\n")
cat("Output:       ", base_out_dir, "\n")
cat("Time window:  ", year_min, "-", year_max, "\n\n")

# --- 3) Load BAI Data --------------------------------------------------------

if (!file_exists(bai_data_path)) {
  stop("BAI data not found. Run 01e first.")
}

if (!file_exists(bai_site_path)) {
  stop("Site summary not found. Run 01e first.")
}

cat("Loading BAI data...\n")

bai_long <- read_csv(bai_data_path, show_col_types = FALSE)
site_meta <- read_csv(bai_site_path, show_col_types = FALSE)

cat("  BAI observations: ", nrow(bai_long), "\n")
cat("  Unique sites:     ", n_distinct(bai_long$site_code), "\n\n")

# --- 4) Aggregate BAI to Site-Year Level ------------------------------------

cat("Aggregating BAI to site-year level...\n")

bai_long <- bai_long %>%
  filter(year >= year_min, year <= year_max)

bai_site_year <- bai_long %>%
  group_by(site_code, year) %>%
  summarise(
    bai_cm2_mean = mean(bai_cm2, na.rm = TRUE),
    bai_cm2_sum  = sum(bai_cm2, na.rm = TRUE),
    bai_cm2_sd   = sd(bai_cm2, na.rm = TRUE),
    n_series     = n_distinct(series_id),
    .groups = "drop"
  )

cat("  Site-year combinations: ", nrow(bai_site_year), "\n\n")

sites_coords <- bai_long %>%
  select(site_code, lat, lon) %>%
  distinct()

cat("  Unique sites for extraction: ", nrow(sites_coords), "\n\n")

# --- 5) Helper Functions -----------------------------------------------------

build_year_index <- function(nlayers, start_year) {
  idx  <- seq_len(nlayers)
  year <- start_year + (idx - 1L) %/% 12L
  year
}

extract_annual_at_sites <- function(rast_ann, sites_df, var_name) {
  years <- names(rast_ann)
  if (is.null(years) || any(years == "")) {
    years <- as.character(seq_len(terra::nlyr(rast_ann)))
  }
  
  vals <- terra::extract(
    rast_ann,
    cbind(sites_df$lon, sites_df$lat)
  )
  
  df_long <- vals %>%
    mutate(ID = row_number()) %>%
    select(ID, everything()) %>%
    pivot_longer(
      cols      = -ID,
      names_to  = "layer",
      values_to = var_name
    )
  
  id_map <- tibble(
    ID        = seq_len(nrow(sites_df)),
    site_code = sites_df$site_code
  )
  
  df_long <- df_long %>%
    left_join(id_map, by = "ID") %>%
    select(site_code, layer, all_of(var_name)) %>%
    mutate(
      year = suppressWarnings(as.integer(str_extract(layer, "\\d+")))
    ) %>%
    filter(!is.na(year), year >= year_min, year <= year_max) %>%
    select(site_code, year, all_of(var_name))
  
  df_long
}

ensure_nc <- function(path_gz) {
  if (!str_detect(path_gz, "\\.gz$")) return(path_gz)
  nc_path <- str_remove(path_gz, "\\.gz$")
  if (!file_exists(nc_path)) {
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      install.packages("R.utils")
    }
    R.utils::gunzip(path_gz, destname = nc_path, overwrite = FALSE, remove = FALSE)
  }
  nc_path
}

# --- 6) Extract CRU TS Temperature & Precipitation ---------------------------

cat("=== EXTRACTING CRU TS v4.09 ===\n")

tmp_gz <- path(cru_dir, "tmp", "cru_ts4.09.1901.2024.tmp.dat.nc.gz")
pre_gz <- path(cru_dir, "pre", "cru_ts4.09.1901.2024.pre.dat.nc.gz")

if (!file_exists(tmp_gz) || !file_exists(pre_gz)) {
  stop("CRU TS files not found. Run 02a first.")
}

tmp_nc <- ensure_nc(tmp_gz)
pre_nc <- ensure_nc(pre_gz)

cat("Reading temperature...\n")
cru_tmp_monthly <- terra::rast(tmp_nc)

cat("Reading precipitation...\n")
cru_pre_monthly <- terra::rast(pre_nc)

bbox_boreal <- terra::ext(-180, 180, 50, 90)
cru_tmp_monthly <- terra::crop(cru_tmp_monthly, bbox_boreal)
cru_pre_monthly <- terra::crop(cru_pre_monthly, bbox_boreal)

year_index_cru <- build_year_index(terra::nlyr(cru_tmp_monthly), start_year = 1901)

cat("Computing annual temperature (mean)...\n")
cru_tmp_annual <- terra::tapp(cru_tmp_monthly, index = year_index_cru, fun = mean)
names(cru_tmp_annual) <- sort(unique(year_index_cru))

cat("Computing annual precipitation (sum)...\n")
cru_pre_annual <- terra::tapp(cru_pre_monthly, index = year_index_cru, fun = sum)
names(cru_pre_annual) <- sort(unique(year_index_cru))

cat("Extracting site values...\n")
tmp_ann_df <- extract_annual_at_sites(cru_tmp_annual, sites_coords, "tmp_C")
pre_ann_df <- extract_annual_at_sites(cru_pre_annual, sites_coords, "pre_mm")

cat("  Temperature site-years: ", nrow(tmp_ann_df), "\n")
cat("  Precipitation site-years: ", nrow(pre_ann_df), "\n\n")

# --- 7) Extract ISIMIP3a N-Deposition ----------------------------------------

cat("=== EXTRACTING ISIMIP3a N-DEPOSITION ===\n")

ndep_var_to_start_year <- list(
  "ndep-nhx" = 1850L,
  "ndep-noy" = 1850L
)

extract_ndep_var <- function(var_name) {
  var_dir <- path(ndep_root, var_name)
  
  if (!dir_exists(var_dir)) {
    warning("N-deposition directory not found: ", var_dir)
    return(NULL)
  }
  
  nc_files <- dir_ls(var_dir, glob = "*.nc")
  
  if (length(nc_files) == 0) {
    warning("No N-deposition files found for: ", var_name)
    return(NULL)
  }
  
  nc_files <- sort(nc_files)
  
  cat("  Variable: ", var_name, "\n")
  for (f in nc_files) cat("    - ", path_file(f), "\n")
  
  ras_list <- list()
  for (f in nc_files) {
    r <- try(terra::rast(f), silent = TRUE)
    if (inherits(r, "try-error")) {
      warning("    Skipping file: ", path_file(f))
      next
    }
    ras_list[[length(ras_list) + 1L]] <- r
  }
  
  if (length(ras_list) == 0) {
    warning("No valid N-deposition rasters for: ", var_name)
    return(NULL)
  }
  
  r_monthly <- do.call(c, ras_list)
  r_monthly <- terra::crop(r_monthly, bbox_boreal)
  
  nl <- terra::nlyr(r_monthly)
  y0 <- ndep_var_to_start_year[[var_name]]
  idxY <- build_year_index(nlayers = nl, start_year = y0)
  
  cat("  Computing annual sums...\n")
  r_annual <- terra::tapp(r_monthly, index = idxY, fun = sum)
  names(r_annual) <- sort(unique(idxY))
  
  df <- extract_annual_at_sites(r_annual, sites_coords, paste0(var_name, "_ann"))
  df
}

ndep_nhx_df <- extract_ndep_var("ndep-nhx")
ndep_noy_df <- extract_ndep_var("ndep-noy")

if (!is.null(ndep_nhx_df)) {
  cat("  NHx site-years: ", nrow(ndep_nhx_df), "\n")
} else {
  cat("  NHx: Not available\n")
}

if (!is.null(ndep_noy_df)) {
  cat("  NOy site-years: ", nrow(ndep_noy_df), "\n")
} else {
  cat("  NOy: Not available\n")
}
cat("\n")

# --- 8) Combine All Data -----------------------------------------------------

cat("=== COMBINING BAI + CLIMATE + N-DEPOSITION ===\n")

site_meta_clean <- site_meta %>%
  select(site_code, site_name, region, lat, lon, elevation_m,
         genus, scientific_name, species_code, year_min, year_max, n_series) %>%
  distinct()

df_lmm <- bai_site_year %>%
  left_join(tmp_ann_df, by = c("site_code", "year")) %>%
  left_join(pre_ann_df, by = c("site_code", "year"))

if (!is.null(ndep_nhx_df)) {
  df_lmm <- df_lmm %>% left_join(ndep_nhx_df, by = c("site_code", "year"))
}

if (!is.null(ndep_noy_df)) {
  df_lmm <- df_lmm %>% left_join(ndep_noy_df, by = c("site_code", "year"))
}

df_lmm <- df_lmm %>%
  left_join(site_meta_clean, by = "site_code") %>%
  arrange(site_code, year) %>%
  filter(year >= year_min, year <= year_max, !is.na(bai_cm2_mean))

cat("  Combined site-year rows: ", nrow(df_lmm), "\n")
cat("  Unique sites: ", n_distinct(df_lmm$site_code), "\n")
cat("  Year range: ", min(df_lmm$year), "-", max(df_lmm$year), "\n\n")

# --- 9) Save Output ----------------------------------------------------------

out_parquet <- path(dir_data, "02b_BAI_climate_annual.parquet")
out_csv     <- path(dir_data, "02b_BAI_climate_annual.csv")

arrow::write_parquet(df_lmm, out_parquet)
write_csv(df_lmm, out_csv)

summary_stats <- df_lmm %>%
  summarise(
    n_observations = n(),
    n_sites = n_distinct(site_code),
    year_min = min(year),
    year_max = max(year),
    mean_bai_cm2 = mean(bai_cm2_mean, na.rm = TRUE),
    mean_tmp_C = mean(tmp_C, na.rm = TRUE),
    mean_pre_mm = mean(pre_mm, na.rm = TRUE)
  )

write_csv(summary_stats, path(dir_tables, "02b_summary_statistics.csv"))

cat("=== OUTPUT FILES ===\n")
cat("  Parquet:   ", out_parquet, "\n")
cat("  CSV:       ", out_csv, "\n")
cat("  Summary:   ", path(dir_tables, "02b_summary_statistics.csv"), "\n")

cat("\n=== COLUMNS ===\n")
cat("  site_code, year\n")
cat("  bai_cm2_mean, bai_cm2_sum, bai_cm2_sd\n")
cat("  tmp_C (annual mean °C)\n")
cat("  pre_mm (annual sum mm)\n")
if (!is.null(ndep_nhx_df)) cat("  ndep-nhx_ann (kg N/ha/yr)\n")
if (!is.null(ndep_noy_df)) cat("  ndep-noy_ann (kg N/ha/yr)\n")
cat("  site_name, region, lat, lon, elevation_m\n")
cat("  genus, scientific_name, species_code\n")

cat("\n=== 02b Prepare BAI Climate for LMM COMPLETE ===\n")
cat(sprintf("Successfully prepared %d site-year observations\n", nrow(df_lmm)))