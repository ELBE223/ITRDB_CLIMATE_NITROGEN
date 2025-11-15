############################################################
# 01_download_ITRDB_RWL_global.R
# Download ITRDB ring-width (.rwl) files by region
# and build simple region + global index tables.
#
# Regions (NOAA directory structure):
#   africa, asia, australia, europe,
#   northamerica/canada, northamerica/mexico,
#   northamerica/usa, southamerica
############################################################

# --- 0) Packages & options ---------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(here)  # project root helper

# Increase timeout for large regions (e.g. Canada, USA)
options(timeout = 600)

# --- 1) Paths & folders ------------------------------------------

base_dir  <- here::here("itrdb_global_example")
rwl_root  <- file.path(base_dir, "rwl")
meta_dir  <- file.path(base_dir, "metadata")

dir.create(base_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(rwl_root,  recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir,  recursive = TRUE, showWarnings = FALSE)

cat("Base directory:\n  ", normalizePath(base_dir), "\n\n")

# --- 2) Global ITRDB metadata text file (optional) ---------------

meta_url  <- "https://paleo.ltrr.arizona.edu/treering/ITRDBmetadata2018Dec17.txt"
meta_path <- file.path(meta_dir, "ITRDBmetadata2018Dec17.txt")

if (!file.exists(meta_path)) {
  cat("Downloading global ITRDB metadata text file...\n")
  utils::download.file(meta_url, destfile = meta_path,
                       mode = "wb", quiet = TRUE)
} else {
  cat("Metadata text file already exists, skipping download.\n")
}
cat("Metadata text file:\n  ", normalizePath(meta_path), "\n\n")

# --- 3) Regions / directories ------------------------------------

regions <- c(
  "africa",
  "asia",
  "australia",
  "europe",
  "northamerica/canada",
  "northamerica/mexico",
  "northamerica/usa",
  "southamerica"
)

base_measure_url <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/"

# Helper: list .rwl files in one region directory
get_region_rwl_files <- function(region) {
  url <- paste0(base_measure_url, region, "/")
  message("Reading directory listing: ", url)
  
  html_lines <- try(readLines(url, warn = FALSE), silent = TRUE)
  if (inherits(html_lines, "try-error")) {
    warning("Could not read directory for region '", region,
            "' (URL: ", url, ")")
    return(character(0))
  }
  
  rwl_lines <- grep("\\.rwl", html_lines, value = TRUE)
  if (!length(rwl_lines)) return(character(0))
  
  rwl_files <- sub('.*href=\"([^\"]+\\.rwl)\".*', "\\1", rwl_lines)
  rwl_files <- unique(rwl_files)
  rwl_files <- rwl_files[!is.na(rwl_files) & rwl_files != ""]
  
  rwl_files
}

# --- 4) Loop over regions: download + build indices --------------

index_list <- list()

for (reg in regions) {
  cat("------------------------------------------------------------\n")
  cat("Region:", reg, "\n")
  
  # Region subfolder under rwl/
  reg_dir <- file.path(rwl_root, reg)
  dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 4a) list .rwl files
  rwl_files <- get_region_rwl_files(reg)
  cat("  Found", length(rwl_files), "RWL files.\n")
  
  if (!length(rwl_files)) next
  
  # 4b) download all .rwl files (robust against errors/timeouts)
  for (f in rwl_files) {
    dest <- file.path(reg_dir, f)
    if (!file.exists(dest)) {
      url <- paste0(base_measure_url, reg, "/", f)
      cat("   Downloading", f, "...\n")
      
      tryCatch(
        {
          utils::download.file(url, destfile = dest,
                               mode = "wb", quiet = TRUE)
        },
        error = function(e) {
          message("   ERROR downloading ", f, ": ",
                  conditionMessage(e))
        }
      )
      
    } else {
      cat("   Skipping (already exists):", f, "\n")
    }
  }
  
  # 4c) build region-level index (even if some files failed)
  reg_index <- data.frame(
    region    = reg,
    site_code = toupper(sub("\\.rwl$", "", rwl_files)),
    rwl_file  = rwl_files,
    rwl_url   = paste0(base_measure_url, reg, "/", rwl_files),
    stringsAsFactors = FALSE
  )
  
  # Save per-region index
  reg_index_path <- file.path(
    base_dir,
    paste0("itrdb_", gsub("/", "_", reg), "_rwl_index.csv")
  )
  utils::write.csv(reg_index, reg_index_path, row.names = FALSE)
  cat("  Region index written to:\n   ",
      normalizePath(reg_index_path), "\n\n")
  
  index_list[[reg]] <- reg_index
}

# --- 5) Global index (all regions combined) -----------------------

if (length(index_list) > 0) {
  global_index <- do.call(rbind, index_list)
  
  global_index_path <- file.path(base_dir, "itrdb_all_regions_rwl_index.csv")
  utils::write.csv(global_index, global_index_path, row.names = FALSE)
  
  cat("============================================================\n")
  cat("Global RWL index written to:\n  ",
      normalizePath(global_index_path), "\n")
  cat("Columns:\n  region, site_code, rwl_file, rwl_url\n")
} else {
  cat("No RWL files found for any region. Please check network/URLs.\n")
}
