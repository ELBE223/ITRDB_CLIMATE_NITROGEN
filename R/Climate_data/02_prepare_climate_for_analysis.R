############################################################
# 02_prepare_BAI_climate_for_LMM.R
#
# Goal:
#  - Combine boreal BAI (site-year) with:
#      * CRU TS v4.09 (tmp, pre, 1901–2024)
#      * ISIMIP3a N deposition (ndep-nhx, ndep-noy, 1850–2021, histsoc)
#  - Restrict all data to the common window 1901–2010
#  - Output an annual site-level table for LMMs:
#      site_code × year × (BAI, tmp, pre, ndep_nhx, ndep_noy, metadata)
#
# Inputs:
#  - itrdb_global_example/data_prepared/itrdb_boreal_site_metadata.parquet
#  - itrdb_global_example/data_prepared/itrdb_boreal_BAI_long.parquet
#  - Input/climate/CRU_TS_v4.09/tmp/cru_ts4.09.1901.2024.tmp.dat.nc.gz
#  - Input/climate/CRU_TS_v4.09/pre/cru_ts4.09.1901.2024.pre.dat.nc.gz
#  - Input/climate/ISIMIP3a/ndep/histsoc/ndep-nhx_*.nc
#  - Input/climate/ISIMIP3a/ndep/histsoc/ndep-noy_*.nc
#
# Outputs:
#  - itrdb_global_example/data_prepared/itrdb_boreal_BAI_climate_annual.parquet
#  - itrdb_global_example/data_prepared/itrdb_boreal_BAI_climate_annual.csv
############################################################

## --- 0) Packages ----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  here,
  dplyr,
  tidyr,
  tibble,
  arrow,
  terra,
  stringr
)

## --- 0b) Common time window ----------------------------------------

year_min <- 1901L
year_max <- 2010L

## --- 1) Paths -------------------------------------------------------

# BAI / tree-ring project
bai_base <- here::here("itrdb_global_example")
prep_dir <- file.path(bai_base, "data_prepared")

boreal_sites_parquet <- file.path(prep_dir, "itrdb_boreal_site_metadata.parquet")
boreal_bai_parquet   <- file.path(prep_dir, "itrdb_boreal_BAI_long.parquet")

if (!file.exists(boreal_sites_parquet)) {
  stop("Boreal site metadata parquet not found: ", boreal_sites_parquet,
       "\nRun 05_extract_BAI_boreal_global.R first.")
}
if (!file.exists(boreal_bai_parquet)) {
  stop("Boreal BAI long parquet not found: ", boreal_bai_parquet,
       "\nRun 05_extract_BAI_boreal_global.R first.")
}

# Climate / N deposition
clim_base <- here::here("Input", "climate")

cru_dir   <- file.path(clim_base, "CRU_TS_v4.09")
ndep_root <- file.path(clim_base, "ISIMIP3a", "ndep", "histsoc")

## --- 2) Load boreal BAI data ---------------------------------------

cat("Loading boreal site metadata and BAI long table...\n")

boreal_sites <- arrow::read_parquet(boreal_sites_parquet) %>%
  as.data.frame()

boreal_bai_long <- arrow::read_parquet(boreal_bai_parquet) %>%
  as.data.frame()

cat("  Boreal sites:  ", nrow(boreal_sites), "\n")
cat("  BAI long rows: ", nrow(boreal_bai_long), "\n\n")

## --- 3) Aggregate BAI to site-year level ----------------------------

cat("Aggregating BAI to site-year level...\n")

# restrict BAI to common window first
boreal_bai_long <- boreal_bai_long %>%
  dplyr::filter(year >= year_min, year <= year_max)

bai_site_year <- boreal_bai_long %>%
  dplyr::group_by(site_code, year) %>%
  dplyr::summarise(
    BAI_mm2_mean = mean(BAI_mm2, na.rm = TRUE),
    BAI_mm2_sum  = sum(BAI_mm2,  na.rm = TRUE),
    .groups = "drop"
  )

cat("  Site-year combinations (BAI, ", year_min, "-", year_max, "): ",
    nrow(bai_site_year), "\n\n")

# unique site list for extraction
sites_coords <- boreal_sites %>%
  dplyr::select(site_code, lat, lon) %>%
  dplyr::distinct()

cat("  Unique boreal sites (for climate extraction): ",
    nrow(sites_coords), "\n\n")

## --- 4) Helper: build annual raster & extract at sites --------------

# create year index for monthly raster
build_year_index <- function(nlayers, start_year) {
  idx  <- seq_len(nlayers)
  year <- start_year + (idx - 1L) %/% 12L
  year
}

# extract annual values at site coordinates, return long table
extract_annual_at_sites <- function(rast_ann, sites_df, var_name) {
  years <- names(rast_ann)
  if (is.null(years) || any(years == "")) {
    years <- as.character(seq_len(terra::nlyr(rast_ann)))
  }
  
  vals <- terra::extract(
    rast_ann,
    cbind(sites_df$lon, sites_df$lat)
  )
  
  # vals: data.frame with ID + layers
  df_long <- vals %>%
    dplyr::mutate(ID = dplyr::row_number()) %>%
    dplyr::select(ID, dplyr::everything()) %>%
    tidyr::pivot_longer(
      cols      = -ID,
      names_to  = "layer",
      values_to = var_name
    )
  
  id_map <- tibble::tibble(
    ID        = seq_len(nrow(sites_df)),
    site_code = sites_df$site_code
  )
  
  df_long <- df_long %>%
    dplyr::left_join(id_map, by = "ID") %>%
    dplyr::select(site_code, layer, dplyr::all_of(var_name)) %>%
    dplyr::mutate(
      year = suppressWarnings(
        as.integer(stringr::str_extract(layer, "\\d+"))
      )
    ) %>%
    dplyr::filter(
      !is.na(year),
      year >= year_min,
      year <= year_max
    ) %>%
    dplyr::select(site_code, year, dplyr::all_of(var_name))
  
  df_long
}

## --- 5) CRU TS: annual temp & precip at sites -----------------------

cat("Reading CRU TS v4.09 tmp & pre and computing annual site values...\n")

# helper to ensure .nc (decompress .gz if needed)
ensure_nc <- function(path_gz) {
  if (!grepl("\\.gz$", path_gz)) return(path_gz)
  nc_path <- sub("\\.gz$", "", path_gz)
  if (!file.exists(nc_path)) {
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      install.packages("R.utils")
    }
    R.utils::gunzip(path_gz, destname = nc_path, overwrite = FALSE, remove = FALSE)
  }
  nc_path
}

tmp_gz <- file.path(cru_dir, "tmp", "cru_ts4.09.1901.2024.tmp.dat.nc.gz")
pre_gz <- file.path(cru_dir, "pre", "cru_ts4.09.1901.2024.pre.dat.nc.gz")

if (!file.exists(tmp_gz) || !file.exists(pre_gz)) {
  stop("CRU TS files not found in ", cru_dir,
       "\nRun 06_download_climate_and_ndep.R first.")
}

tmp_nc <- ensure_nc(tmp_gz)
pre_nc <- ensure_nc(pre_gz)

cat("  Reading CRU tmp from:\n    ", tmp_nc, "\n")
cru_tmp_monthly <- terra::rast(tmp_nc)

cat("  Reading CRU pre from:\n    ", pre_nc, "\n")
cru_pre_monthly <- terra::rast(pre_nc)

# optional crop to boreal band to save memory
bbox_boreal <- terra::ext(-180, 180, 50, 90)
cru_tmp_monthly <- terra::crop(cru_tmp_monthly, bbox_boreal)
cru_pre_monthly <- terra::crop(cru_pre_monthly, bbox_boreal)

# build year index (CRU TS 4.09 monthly 1901–2024)
year_index_cru <- build_year_index(terra::nlyr(cru_tmp_monthly), start_year = 1901)

# annual mean temperature (°C)
cru_tmp_annual <- terra::tapp(
  cru_tmp_monthly,
  index = year_index_cru,
  fun   = mean
)
names(cru_tmp_annual) <- sort(unique(year_index_cru))

# annual precipitation sum (mm/year)
cru_pre_annual <- terra::tapp(
  cru_pre_monthly,
  index = year_index_cru,
  fun   = sum
)
names(cru_pre_annual) <- sort(unique(year_index_cru))

# extract site-year values (and clip to 1901–2010 inside helper)
tmp_ann_df <- extract_annual_at_sites(cru_tmp_annual, sites_coords, "tmp_C")
pre_ann_df <- extract_annual_at_sites(cru_pre_annual, sites_coords, "pre_mm")

cat("  CRU tmp site-year rows: ", nrow(tmp_ann_df), "\n")
cat("  CRU pre site-year rows: ", nrow(pre_ann_df), "\n\n")

## --- 6) ISIMIP3a N-deposition: annual sums at sites -----------------

cat("Reading ISIMIP3a N-deposition (histsoc) and computing annual sums...\n")

ndep_var_to_start_year <- list(
  "ndep-nhx" = 1850L,
  "ndep-noy" = 1850L
)

extract_ndep_var <- function(var_name) {
  var_dir <- file.path(ndep_root, var_name)
  nc_files <- list.files(var_dir, pattern = "\\.nc$", full.names = TRUE)
  if (!length(nc_files)) {
    warning("No N-deposition files found for var: ", var_name,
            " in ", var_dir)
    return(NULL)
  }
  nc_files <- sort(nc_files)
  
  cat("  Variable ", var_name, " files:\n")
  for (f in nc_files) cat("    - ", basename(f), "\n")
  
  # read all files as SpatRaster, skip broken/HTML files
  ras_list <- list()
  for (f in nc_files) {
    r <- try(terra::rast(f), silent = TRUE)
    if (inherits(r, "try-error")) {
      warning("    Skipping file (cannot open as SpatRaster): ", f)
      next
    }
    ras_list[[length(ras_list) + 1L]] <- r
  }
  
  if (!length(ras_list)) {
    warning("No valid N-deposition rasters for var: ", var_name,
            " (all files failed).")
    return(NULL)
  }
  
  r_monthly <- do.call(c, ras_list)
  
  # crop to boreal band
  r_monthly <- terra::crop(r_monthly, bbox_boreal)
  
  nl   <- terra::nlyr(r_monthly)
  y0   <- ndep_var_to_start_year[[var_name]]
  idxY <- build_year_index(nlayers = nl, start_year = y0)
  
  # annual sum
  r_annual <- terra::tapp(
    r_monthly,
    index = idxY,
    fun   = sum
  )
  names(r_annual) <- sort(unique(idxY))
  
  df <- extract_annual_at_sites(r_annual, sites_coords, paste0(var_name, "_ann"))
  df
}

ndep_nhx_df <- extract_ndep_var("ndep-nhx")
ndep_noy_df <- extract_ndep_var("ndep-noy")

if (!is.null(ndep_nhx_df)) {
  cat("  Ndep NHx site-year rows: ", nrow(ndep_nhx_df), "\n")
} else {
  cat("  Ndep NHx could not be read -> will be omitted.\n")
}

if (!is.null(ndep_noy_df)) {
  cat("  Ndep NOy site-year rows: ", nrow(ndep_noy_df), "\n")
} else {
  cat("  Ndep NOy could not be read -> will be omitted.\n")
}
cat("\n")

## --- 7) Combine BAI + climate + Ndep -------------------------------

cat("Combining BAI, CRU, and N-deposition into one table...\n")

# site-level metadata (one row per site_code)
site_meta_small <- boreal_sites %>%
  dplyr::select(
    site_code, site_name, region,
    lat, lon, elevation_m,
    species_code, species_name,
    first_year, last_year, n_series
  ) %>%
  dplyr::distinct()

df_lmm <- bai_site_year %>%
  dplyr::left_join(tmp_ann_df,   by = c("site_code", "year")) %>%
  dplyr::left_join(pre_ann_df,   by = c("site_code", "year")) %>%
  { if (!is.null(ndep_nhx_df)) dplyr::left_join(., ndep_nhx_df,
                                                by = c("site_code", "year")) else . } %>%
  { if (!is.null(ndep_noy_df)) dplyr::left_join(., ndep_noy_df,
                                                by = c("site_code", "year")) else . } %>%
  dplyr::left_join(site_meta_small, by = "site_code") %>%
  dplyr::arrange(site_code, year) %>%
  dplyr::filter(
    year >= year_min, year <= year_max,
    !is.na(BAI_mm2_mean)
  )

cat("  Combined site-year rows (BAI + covariates, ",
    year_min, "-", year_max, "): ",
    nrow(df_lmm), "\n\n")

## --- 8) Save for LMM analysis --------------------------------------

out_parquet <- file.path(prep_dir, "itrdb_boreal_BAI_climate_annual.parquet")
out_csv     <- file.path(prep_dir, "itrdb_boreal_BAI_climate_annual.csv")

arrow::write_parquet(df_lmm, out_parquet)
utils::write.csv(df_lmm, out_csv, row.names = FALSE)

cat("Saved annual BAI+climate+Ndep table as:\n")
cat("  Parquet: ", normalizePath(out_parquet), "\n")
cat("  CSV:     ", normalizePath(out_csv), "\n\n")

cat("Columns (core):\n",
    "  site_code, year,\n",
    "  BAI_mm2_mean, BAI_mm2_sum,\n",
    "  tmp_C (annual mean), pre_mm (annual sum),\n",
    "  ndep-nhx_ann (annual NHx, if available),\n",
    "  ndep-noy_ann (annual NOy, if available),\n",
    "  site_name, region, lat, lon, elevation_m,\n",
    "  species_code, species_name, first_year, last_year, n_series\n")

cat("\nDone (07_prepare_BAI_climate_for_LMM.R).\n")
############################################################
