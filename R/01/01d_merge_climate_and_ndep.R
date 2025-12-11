###############################################################################
# 01d_merge_climate_and_ndep.R
#
# Description:
#   1. Load cleaned BAI data (01c).
#   2. Extract Climate Data (CRU TS) for each site (1901-2010).
#      - UNZIPS .gz files automatically if needed.
#      - Calc Seasonal Means (Winter=DJF, Spring=MAM, etc.)
#   3. Extract N-Deposition Data (ISIMIP3a) for each site (1901-2010).
#      - Calc Total N, and 3-year/5-year rolling means.
#      - Robust time handling via NetCDF "time" units.
#   4. Merge everything into one master analysis dataframe.
#
# Output:
#   - output/01d_merge_climate_and_ndep/data/01d_Analysis_Data_Merged.rds
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  tidyr,
  terra,
  zoo,
  purrr,
  stringr,
  R.utils,
  ncdf4
)

# --- 2) Directories & Paths --------------------------------------------------

step_00b <- "00b_download_ISIMIP_CRU"
step_01c <- "01c_QC_filter_merge"
current  <- "01d_merge_climate_and_ndep"

base_out_dir <- here::here("output", current)
dir_data     <- path(base_out_dir, "data")
dir_create(dir_data)

bai_rds_path <- here::here("output", step_01c, "data", "01c_BAI_long_filtered.rds")
cru_base     <- here::here("output", step_00b, "data", "climate", "CRU_TS_v4.09")
isimip_base  <- here::here("output", step_00b, "data", "climate", "ISIMIP3a", "ndep")

cat("=== 01d Merge Climate and N-Deposition ===\n\n")

# --- 3) Load BAI Data --------------------------------------------------------

safe_read_bai <- function(path) {
  if (!file_exists(path)) {
    stop("BAI data from 01c not found at: ", path)
  }
  
  fsize <- file_info(path)$size
  cat(sprintf("Checking input file: %.2f MB\n", fsize / 1024^2))
  
  if (is.na(fsize) || fsize < 1000) {
    stop("The file '01c_BAI_long_filtered.rds' is empty or corrupt (size < 1KB). Please re-run script 01c!")
  }
  
  magic <- tryCatch(tools::file_magic(path), error = function(e) NA_character_)
  if (is.na(magic)) {
    warning("tools::file_magic() could not read header of BAI file. File may be corrupt.")
  }
  
  out <- tryCatch(
    readRDS(path),
    error = function(e) {
      stop(
        "readRDS() failed for '", path,
        "'. The file is likely corrupt or not a valid RDS.\n",
        "Original error: ", conditionMessage(e),
        "\nPlease re-run 01c_QC_filter_merge to regenerate the RDS."
      )
    }
  )
  
  out
}

cat("Loading BAI data...\n")
df_bai <- safe_read_bai(bai_rds_path) %>%
  filter(year >= 1901, year <= 2010)

site_coords <- df_bai %>%
  distinct(site_code, lat, lon) %>%
  filter(!is.na(lat), !is.na(lon))

vect_sites <- terra::vect(site_coords, geom = c("lon", "lat"), crs = "EPSG:4326")

cat(sprintf("Processing %d unique sites.\n", nrow(site_coords)))

# --- 4) Helper: Unzip Function ----------------------------------------------

get_ready_nc <- function(folder, pattern_glob) {
  f_nc <- dir_ls(folder, glob = str_remove(pattern_glob, "\\.gz"))
  if (length(f_nc) > 0) return(f_nc[1])
  
  f_gz <- dir_ls(folder, glob = pattern_glob)
  if (length(f_gz) > 0) {
    cat(sprintf("  Unzipping %s ...\n", path_file(f_gz[1])))
    R.utils::gunzip(f_gz[1], remove = FALSE, skip = TRUE)
    return(str_remove(f_gz[1], "\\.gz$"))
  }
  
  return(NA_character_)
}

# --- 5) Process CRU Climate --------------------------------------------------

process_cru <- function() {
  cat("\n--- Processing CRU TS Climate ---\n")
  
  file_tmp <- get_ready_nc(path(cru_base, "tmp"), "*.nc.gz")
  file_pre <- get_ready_nc(path(cru_base, "pre"), "*.nc.gz")
  
  if (is.na(file_tmp) || is.na(file_pre)) stop("CRU NC files not found.")
  
  r_tmp <- terra::rast(file_tmp)
  r_pre <- terra::rast(file_pre)
  
  cat("Extracting CRU data...\n")
  ext_tmp <- t(terra::extract(r_tmp, vect_sites, ID = FALSE))
  ext_pre <- t(terra::extract(r_pre, vect_sites, ID = FALSE))
  
  dates  <- terra::time(r_tmp)
  years  <- as.numeric(format(dates, "%Y"))
  months <- as.numeric(format(dates, "%m"))
  
  cat("Calculating seasonal means per site...\n")
  
  climate_list <- vector("list", length = nrow(site_coords))
  
  for (i in seq_len(nrow(site_coords))) {
    site <- site_coords$site_code[i]
    
    df_site <- tibble::tibble(
      year  = years,
      month = months,
      tmp   = ext_tmp[, i],
      pre   = ext_pre[, i]
    ) %>%
      filter(year >= 1900, year <= 2010)
    
    df_seasonal <- df_site %>%
      arrange(year, month) %>%
      mutate(
        tmp_prev = dplyr::lag(tmp, 1),
        pre_prev = dplyr::lag(pre, 1),
        season   = dplyr::case_when(
          month %in% c(3, 4, 5)   ~ "spring",
          month %in% c(6, 7, 8)   ~ "summer",
          month %in% c(9, 10, 11) ~ "autumn",
          TRUE                    ~ "winter_calc"
        )
      ) %>%
      group_by(year) %>%
      summarise(
        temp_annual = mean(tmp, na.rm = TRUE),
        prec_annual = sum(pre,  na.rm = TRUE),
        
        temp_spring = mean(tmp[season == "spring"], na.rm = TRUE),
        prec_spring = sum(pre[season == "spring"],  na.rm = TRUE),
        temp_summer = mean(tmp[season == "summer"], na.rm = TRUE),
        prec_summer = sum(pre[season == "summer"],  na.rm = TRUE),
        temp_autumn = mean(tmp[season == "autumn"], na.rm = TRUE),
        prec_autumn = sum(pre[season == "autumn"],  na.rm = TRUE),
        
        t_jan      = mean(tmp[month == 1],      na.rm = TRUE),
        t_feb      = mean(tmp[month == 2],      na.rm = TRUE),
        t_dec_prev = mean(tmp_prev[month == 1], na.rm = TRUE),
        p_jan      = sum(pre[month == 1],       na.rm = TRUE),
        p_feb      = sum(pre[month == 2],       na.rm = TRUE),
        p_dec_prev = sum(pre_prev[month == 1],  na.rm = TRUE),
        .groups    = "drop"
      ) %>%
      mutate(
        temp_winter = (t_dec_prev + t_jan + t_feb) / 3,
        prec_winter = p_dec_prev + p_jan + p_feb
      ) %>%
      select(year, starts_with("temp_"), starts_with("prec_")) %>%
      mutate(site_code = site)
    
    climate_list[[i]] <- df_seasonal
  }
  
  bind_rows(climate_list)
}

df_climate <- process_cru()

# --- 6) Process ISIMIP N-Deposition (ncdf4) ----------------------------------

process_isimip <- function() {
  cat("\n--- Processing ISIMIP N-Deposition (ncdf4) ---\n")
  
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for ISIMIP processing. Please install.packages('ncdf4').")
  }
  
  vars <- c("ndep-nhx", "ndep-noy")
  all_ndep_list <- list()
  
  for (v in vars) {
    cat(sprintf("  Variable pattern: %s\n", v))
    
    pattern <- paste0("*", v, "*histsoc*.nc*")
    files_found <- dir_ls(isimip_base, recurse = TRUE, glob = pattern)
    
    if (length(files_found) == 0) {
      stop("No ISIMIP files found for pattern: ", pattern)
    }
    
    files_1901 <- files_found[grepl("1901_2021", files_found)]
    f <- if (length(files_1901) > 0) files_1901[1] else files_found[1]
    
    cat("    Using file: ", f, "\n")
    
    nc <- ncdf4::nc_open(f)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    
    dimnames <- names(nc$dim)
    lon_name <- if ("lon" %in% dimnames) "lon" else if ("x" %in% dimnames) "x" else
      stop("Could not find lon/x dimension in NetCDF: ", f)
    lat_name <- if ("lat" %in% dimnames) "lat" else if ("y" %in% dimnames) "y" else
      stop("Could not find lat/y dimension in NetCDF: ", f)
    
    lon <- ncdf4::ncvar_get(nc, lon_name)
    lat <- ncdf4::ncvar_get(nc, lat_name)
    
    time_vals  <- ncdf4::ncvar_get(nc, "time")
    time_units <- ncdf4::ncatt_get(nc, "time", "units")$value
    
    # --- ROBUST TIME HANDLING -------------------------------------------------
    # Beispiel units:
    # "days since 1901-01-01 00:00:00", "months since 1901-01-01", etc.
    origin_raw <- sub("^.*since\\s+", "", time_units)
    origin_raw <- strsplit(origin_raw, "\n")[[1]][1]
    origin_date_str <- strsplit(origin_raw, " ")[[1]][1]
    
    origin_date <- tryCatch(
      as.Date(origin_date_str),
      error = function(e) NA
    )
    
    if (!is.na(origin_date) && grepl("day", time_units, ignore.case = TRUE)) {
      dates_raw <- origin_date + time_vals
      years     <- as.numeric(format(dates_raw, "%Y"))
    } else {
      origin_year <- suppressWarnings(as.integer(substr(origin_date_str, 1, 4)))
      if (is.na(origin_year)) {
        stop(
          "Could not parse time units for file: ", f,
          " (units: ", time_units, ")"
        )
      }
      
      if (grepl("month", time_units, ignore.case = TRUE)) {
        years <- origin_year + floor(time_vals / 12)
      } else if (grepl("day", time_units, ignore.case = TRUE)) {
        years <- origin_year + floor(time_vals / 365.25)
      } else {
        stop(
          "Unsupported time units in file: ", f,
          " (units: ", time_units, ")"
        )
      }
    }
    # --------------------------------------------------------------------------
    
    nc_vars     <- names(nc$var)
    var_name_nc <- if ("nhx" %in% nc_vars) {
      "nhx"
    } else if ("noy" %in% nc_vars) {
      "noy"
    } else if (length(nc_vars) == 1) {
      nc_vars[1]
    } else {
      stop("Could not detect NHx/NOy variable in NetCDF: ", f)
    }
    
    cat("    Found data variable: ", var_name_nc, "\n")
    
    get_index <- function(val, grid) which.min(abs(grid - val))
    
    n_sites <- nrow(site_coords)
    site_lon_idx <- integer(n_sites)
    site_lat_idx <- integer(n_sites)
    
    for (i in seq_len(n_sites)) {
      site_lon_idx[i] <- get_index(site_coords$lon[i], lon)
      site_lat_idx[i] <- get_index(site_coords$lat[i], lat)
    }
    
    chunk_list <- vector("list", n_sites)
    
    for (i in seq_len(n_sites)) {
      ts_vals <- ncdf4::ncvar_get(
        nc,
        var_name_nc,
        start = c(site_lon_idx[i], site_lat_idx[i], 1),
        count = c(1, 1, -1)
      )
      
      df_raw <- tibble::tibble(
        year = years,
        val  = as.numeric(ts_vals)
      )
      
      df_ann <- df_raw %>%
        group_by(year) %>%
        summarise(
          annual_sum = sum(val, na.rm = TRUE),
          .groups    = "drop"
        ) %>%
        mutate(site_code = site_coords$site_code[i])
      
      chunk_list[[i]] <- df_ann
    }
    
    full_var_df <- bind_rows(chunk_list) %>%
      arrange(site_code, year) %>%
      distinct()
    
    all_ndep_list[[v]] <- full_var_df
  }
  
  df_nhx <- all_ndep_list[["ndep-nhx"]] %>% rename(nhx = annual_sum)
  df_noy <- all_ndep_list[["ndep-noy"]] %>% rename(noy = annual_sum)
  
  df_ndep <- df_nhx %>%
    full_join(df_noy, by = c("site_code", "year")) %>%
    mutate(ndep_total = nhx + noy) %>%
    arrange(site_code, year) %>%
    group_by(site_code) %>%
    mutate(
      ndep_total_3yr = zoo::rollmean(ndep_total, k = 3, fill = NA, align = "right"),
      ndep_total_5yr = zoo::rollmean(ndep_total, k = 5, fill = NA, align = "right")
    ) %>%
    ungroup() %>%
    filter(year >= 1901, year <= 2010)
  
  df_ndep
}

df_ndep <- process_isimip()

# --- 7) Merge & Save ---------------------------------------------------------

cat("\n--- Merging all data ---\n")

df_final <- df_bai %>%
  left_join(df_climate, by = c("site_code", "year")) %>%
  left_join(df_ndep,     by = c("site_code", "year")) %>%
  filter(!is.na(temp_annual))

cat(sprintf("Final dimensions: %d rows, %d cols.\n", nrow(df_final), ncol(df_final)))

out_file <- path(dir_data, "01d_Analysis_Data_Merged.rds")
saveRDS(df_final, out_file)

cat("\n=== 01d Processing Complete ===\n")
cat("Merged Dataset saved to: ", out_file, "\n")
