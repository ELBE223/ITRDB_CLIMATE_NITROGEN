###############################################################################
# 02a_download_climate_and_ndep.R
#
# Description:
#   Download climate and nitrogen deposition data for environmental analysis.
#   - CRU TS v4.09 (temperature + precipitation, 1901-2024)
#   - ISIMIP3a nitrogen deposition (ndep-nhx, ndep-noy, 1850-2021)
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, fs)

options(timeout = max(3600, getOption("timeout")))

# --- 2) Directories & Paths --------------------------------------------------

current_name <- "02a_download_climate_and_ndep"

base_out_dir <- here::here("output", current_name)
dir_data     <- path(base_out_dir, "data")
dir_tables   <- path(base_out_dir, "tables")
dir_figures  <- path(base_out_dir, "figures")

clim_dir   <- path(dir_data, "climate")
cru_dir    <- path(clim_dir, "CRU_TS_v4.09")
isimip_dir <- path(clim_dir, "ISIMIP3a")

dir_create(c(dir_tables, dir_figures, clim_dir, cru_dir, isimip_dir))

cat("=== 02a Download Climate and N-deposition ===\n\n")
cat("Output: ", base_out_dir, "\n")
cat("Climate data: ", clim_dir, "\n\n")

# --- 3) Safe Download Function -----------------------------------------------

safe_download <- function(url, dest, tries = 3, sleep_sec = 5) {
  if (file_exists(dest)) {
    cat("File exists, skipping: ", path_file(dest), "\n")
    return(invisible(TRUE))
  }
  
  dir_create(path_dir(dest))
  
  sys <- Sys.info()[["sysname"]]
  method <- ifelse(sys %in% c("Darwin", "Linux"), "curl", "auto")
  
  for (k in seq_len(tries)) {
    cat(sprintf("  Downloading (try %d/%d): %s\n", k, tries, path_file(dest)))
    
    ok <- try(
      utils::download.file(url, destfile = dest, mode = "wb", 
                           quiet = FALSE, method = method),
      silent = TRUE
    )
    
    if (!inherits(ok, "try-error") && file_exists(dest) && 
        file.info(dest)$size > 1e5) {
      size_mb <- round(file.info(dest)$size / 1024^2, 1)
      cat(sprintf("  -> OK: %s (%.1f MB)\n", path_file(dest), size_mb))
      return(invisible(TRUE))
    }
    
    cat("  -> Failed\n")
    if (k < tries) Sys.sleep(sleep_sec)
  }
  
  warning(sprintf("Failed to download after %d tries: %s", tries, url))
  if (file_exists(dest)) file_delete(dest)
  invisible(FALSE)
}

# --- 4) Download CRU TS v4.09 ------------------------------------------------

cat("=== DOWNLOADING CRU TS v4.09 ===\n\n")

cru_root_url <- "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.09/"
cru_ver_dir <- "cruts.2503051245.v4.09"
cru_vars <- c("tmp", "pre")

cru_downloads <- tibble::tibble(
  variable = character(),
  filename = character(),
  success = logical(),
  size_mb = numeric()
)

for (v in cru_vars) {
  cat("Variable: ", v, "\n")
  
  nc_name <- sprintf("cru_ts4.09.1901.2024.%s.dat.nc.gz", v)
  url <- paste0(cru_root_url, cru_ver_dir, "/", v, "/", nc_name)
  
  var_dir <- path(cru_dir, v)
  dir_create(var_dir)
  dest <- path(var_dir, nc_name)
  
  success <- safe_download(url, dest)
  
  size_mb <- if (file_exists(dest)) {
    round(file.info(dest)$size / 1024^2, 1)
  } else {
    NA_real_
  }
  
  cru_downloads <- bind_rows(
    cru_downloads,
    tibble::tibble(
      variable = v,
      filename = nc_name,
      success = success,
      size_mb = size_mb
    )
  )
  
  cat("\n")
}

cat("CRU TS v4.09 download complete.\n\n")

# --- 5) Download ISIMIP3a N-deposition ---------------------------------------

cat("=== DOWNLOADING ISIMIP3a N-DEPOSITION ===\n\n")

isimip_ndep_base <- "https://files.isimip.org/ISIMIP3a/InputData/socioeconomic/n-deposition"
isimip_ndep_scenarios <- c("histsoc")
isimip_ndep_vars <- c("ndep-nhx", "ndep-noy")
isimip_ndep_year_blocks <- list(c(1850, 1900), c(1901, 2021))

ndep_root <- path(isimip_dir, "ndep")
dir_create(ndep_root)

ndep_downloads <- tibble::tibble(
  scenario = character(),
  variable = character(),
  year_range = character(),
  filename = character(),
  success = logical(),
  size_mb = numeric()
)

for (scen in isimip_ndep_scenarios) {
  cat("Scenario: ", scen, "\n")
  
  for (var in isimip_ndep_vars) {
    cat("  Variable: ", var, "\n")
    
    for (rng in isimip_ndep_year_blocks) {
      y1 <- rng[1]
      y2 <- rng[2]
      
      fname <- sprintf("%s_%s_monthly_%d_%d.nc", var, scen, y1, y2)
      url <- sprintf("%s/%s/%s", isimip_ndep_base, scen, fname)
      dest <- path(ndep_root, scen, var, fname)
      
      success <- safe_download(url, dest)
      
      size_mb <- if (file_exists(dest)) {
        round(file.info(dest)$size / 1024^2, 1)
      } else {
        NA_real_
      }
      
      ndep_downloads <- bind_rows(
        ndep_downloads,
        tibble::tibble(
          scenario = scen,
          variable = var,
          year_range = sprintf("%d-%d", y1, y2),
          filename = fname,
          success = success,
          size_mb = size_mb
        )
      )
    }
  }
  cat("\n")
}

cat("ISIMIP3a N-deposition download complete.\n\n")

# --- 6) Save Download Summary ------------------------------------------------

cat("=== DOWNLOAD SUMMARY ===\n\n")

cru_summary <- cru_downloads %>%
  summarise(
    total_files = n(),
    successful = sum(success),
    failed = sum(!success),
    total_size_mb = sum(size_mb, na.rm = TRUE)
  )

cat("CRU TS v4.09:\n")
cat(sprintf("  Total files:    %d\n", cru_summary$total_files))
cat(sprintf("  Successful:     %d\n", cru_summary$successful))
cat(sprintf("  Failed:         %d\n", cru_summary$failed))
cat(sprintf("  Total size:     %.1f MB\n\n", cru_summary$total_size_mb))

ndep_summary <- ndep_downloads %>%
  summarise(
    total_files = n(),
    successful = sum(success),
    failed = sum(!success),
    total_size_mb = sum(size_mb, na.rm = TRUE)
  )

cat("ISIMIP3a N-deposition:\n")
cat(sprintf("  Total files:    %d\n", ndep_summary$total_files))
cat(sprintf("  Successful:     %d\n", ndep_summary$successful))
cat(sprintf("  Failed:         %d\n", ndep_summary$failed))
cat(sprintf("  Total size:     %.1f MB\n\n", ndep_summary$total_size_mb))

write_csv(cru_downloads, path(dir_tables, "02a_CRU_download_summary.csv"))
write_csv(ndep_downloads, path(dir_tables, "02a_ISIMIP3a_download_summary.csv"))

# --- 7) File Inventory -------------------------------------------------------

all_files <- tibble::tibble(
  type = character(),
  path = character(),
  size_mb = numeric()
)

cru_files <- dir_ls(cru_dir, recurse = TRUE, type = "file")
if (length(cru_files) > 0) {
  all_files <- bind_rows(
    all_files,
    tibble::tibble(
      type = "CRU_TS",
      path = as.character(cru_files),
      size_mb = round(file.info(cru_files)$size / 1024^2, 1)
    )
  )
}

ndep_files <- dir_ls(ndep_root, recurse = TRUE, type = "file", glob = "*.nc")
if (length(ndep_files) > 0) {
  all_files <- bind_rows(
    all_files,
    tibble::tibble(
      type = "ISIMIP3a",
      path = as.character(ndep_files),
      size_mb = round(file.info(ndep_files)$size / 1024^2, 1)
    )
  )
}

write_csv(all_files, path(dir_tables, "02a_file_inventory.csv"))

cat("=== OUTPUT FILES ===\n")
cat("  CRU summary:      ", path(dir_tables, "02a_CRU_download_summary.csv"), "\n")
cat("  ISIMIP3a summary: ", path(dir_tables, "02a_ISIMIP3a_download_summary.csv"), "\n")
cat("  File inventory:   ", path(dir_tables, "02a_file_inventory.csv"), "\n")
cat("  Climate data:     ", clim_dir, "\n")

cat("\n=== 02a Download Climate and N-deposition COMPLETE ===\n")
cat(sprintf("Total files downloaded: %d\n", nrow(all_files)))
cat(sprintf("Total data size: %.1f MB\n", sum(all_files$size_mb, na.rm = TRUE))) 
