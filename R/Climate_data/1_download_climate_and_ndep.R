############################################################
# 06_download_climate_and_ndep.R
#
# Download:
#  - CRU TS v4.09 (temperature + precipitation, 1901–2024)
#  - ISIMIP3a nitrogen deposition (ndep-nhx, ndep-noy, 1850–2021)
#
# Output structure:
#  Input/
#    climate/
#      CRU_TS_v4.09/
#        tmp/cru_ts4.09.1901.2024.tmp.dat.nc.gz
#        pre/cru_ts4.09.1901.2024.pre.dat.nc.gz
#      ISIMIP3a/
#        ndep/<scenario>/<var>/...nc
############################################################

## --- 0) Packages ----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(here)

# increase timeout for large files
options(timeout = max(3600, getOption("timeout")))

## --- 1) Paths -------------------------------------------------------

# you can change "Input" if you prefer another root
base_dir   <- here::here("Input")

clim_dir   <- file.path(base_dir, "climate")
cru_dir    <- file.path(clim_dir, "CRU_TS_v4.09")
isimip_dir <- file.path(clim_dir, "ISIMIP3a")

dir.create(clim_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(cru_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(isimip_dir, recursive = TRUE, showWarnings = FALSE)

cat("Base climate directory:\n  ", normalizePath(clim_dir), "\n\n")

## --- 2) Helper: safe download with retries --------------------------

safe_download <- function(url, dest, tries = 3, sleep_sec = 5) {
  if (file.exists(dest)) {
    message("File already exists, skipping: ", dest)
    return(invisible(TRUE))
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  
  sys <- Sys.info()[["sysname"]]
  method <- ifelse(sys %in% c("Darwin", "Linux"), "curl", "auto")
  
  for (k in seq_len(tries)) {
    message("Downloading (try ", k, "/", tries, "): ", url)
    ok <- try(
      utils::download.file(
        url,
        destfile = dest,
        mode     = "wb",
        quiet    = FALSE,
        method   = method
      ),
      silent = TRUE
    )
    if (!inherits(ok, "try-error") &&
        file.exists(dest) &&
        is.finite(file.info(dest)$size) &&
        file.info(dest)$size > 1e5) {  # > 0.1 MB sanity check
      message("  -> OK: ", dest,
              " (", round(file.info(dest)$size / 1024^2, 1), " MB)")
      return(invisible(TRUE))
    }
    
    # on failure
    cond_msg <- tryCatch(conditionMessage(attr(ok, "condition")),
                         error = function(e) "unknown error")
    message("  -> failed: ", cond_msg)
    if (k < tries) Sys.sleep(sleep_sec)
  }
  
  warning("Failed to download after ", tries, " tries: ", url)
  if (file.exists(dest)) unlink(dest)
  invisible(FALSE)
}

# ====================================================================
# PART A: CRU TS (v4.09, temperature + precipitation)
# ====================================================================

cat("=== CRU TS v4.09 download ======================================\n")

# Base directory for CRU TS v4.09
cru_root_url <- "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.09/"

# Version subdirectory for v4.09 (see CRU TS 4.09 page)
cru_ver_dir <- "cruts.2503051245.v4.09"

# Variables: mean temperature ("tmp") and precipitation ("pre")
cru_vars <- c("tmp", "pre")

for (v in cru_vars) {
  cat("---- Variable:", v, "-----------------------------------------\n")
  
  # use full-period NetCDF (1901–2024)
  nc_name <- sprintf("cru_ts4.09.1901.2024.%s.dat.nc.gz", v)
  
  url  <- paste0(cru_root_url, cru_ver_dir, "/", v, "/", nc_name)
  var_dir <- file.path(cru_dir, v)
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(var_dir, nc_name)
  
  safe_download(url, dest)
  cat("\n")
}

cat("CRU TS v4.09 download finished. Files in:\n  ",
    normalizePath(cru_dir), "\n\n")

# ====================================================================
# PART B: ISIMIP3a nitrogen deposition (ndep-nhx, ndep-noy)
# ====================================================================

cat("=== ISIMIP3a N-deposition download ==============================\n")


# Base URL for ISIMIP3a N-deposition (Yang & Tian 2023)
isimip_ndep_base <- "https://files.isimip.org/ISIMIP3a/InputData/socioeconomic/n-deposition"

# Scenarios you want (you can add "1901soc", "2015soc")
isimip_ndep_scenarios <- c("histsoc")  # c("histsoc", "1901soc", "2015soc")

# Variables: reduced nitrogen (NHx) + oxidized nitrogen (NOy)
isimip_ndep_vars <- c("ndep-nhx", "ndep-noy")

# Year blocks from the dataset description
isimip_ndep_year_blocks <- list(
  c(1850, 1900),
  c(1901, 2021)
)

ndep_root <- file.path(isimip_dir, "ndep")
dir.create(ndep_root, recursive = TRUE, showWarnings = FALSE)

n_ok   <- 0L
n_fail <- 0L

for (scen in isimip_ndep_scenarios) {
  cat("\n--- scenario:", scen, "--------------------------------------\n")
  for (var in isimip_ndep_vars) {
    cat("  variable:", var, "\n")
    for (rng in isimip_ndep_year_blocks) {
      y1 <- rng[1]
      y2 <- rng[2]
      
      # File pattern:
      #   ndep-nhx_histsoc_monthly_1850_1900.nc
      #   ndep-nhx_histsoc_monthly_1901_2021.nc
      fname <- sprintf("%s_%s_monthly_%d_%d.nc", var, scen, y1, y2)
      url   <- sprintf("%s/%s/%s", isimip_ndep_base, scen, fname)
      
      dest  <- file.path(ndep_root, scen, var, fname)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      
      ok <- safe_download(url, dest)
      if (isTRUE(ok)) n_ok <- n_ok + 1L else n_fail <- n_fail + 1L
    }
  }
}

# List all downloaded NetCDFs
ndep_files <- list.files(ndep_root, pattern = "\\.nc$", recursive = TRUE,
                         full.names = TRUE)

cat("\nISIMIP3a N-deposition summary:\n")
cat("  Successfully downloaded files:", n_ok, "\n")
cat("  Failed downloads:             ", n_fail, "\n")
cat("  Total .nc files on disk:      ", length(ndep_files), "\n")
if (length(ndep_files)) {
  cat("  Example file:\n   ",
      normalizePath(ndep_files[1]), "\n")
}

cat("\nDone (06_download_climate_and_ndep.R).\n")
############################################################
