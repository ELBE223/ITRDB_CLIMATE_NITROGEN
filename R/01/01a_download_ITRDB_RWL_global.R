###############################################################################
# 01a_download_ITRDB_RWL_global.R
#
# Description:
#   Download ITRDB ring-width (.rwl) files by region using parallel processing.
#   Downloads BOTH:
#     - classic Tucson .rwl files (e.g., alge001.rwl)
#     - NOAA template .rwl files (e.g., alge001-noaa.rwl)
#   Builds a structured output folder system and saves metadata indices.
#
# Structure:
#   output/01a_download_ITRDB_RWL_global/
#       ├── data/
#       │    └── rwl/        (Raw .rwl files, both Tucson & NOAA)
#       ├── tables/          (CSV indices of downloaded files)
#       └── figures/         (Placeholders for future plots)
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,      # Path management
  fs,        # File system
  dplyr,     # Data manipulation
  stringr,   # String tools
  furrr,     # Parallel processing
  readr      # Fast CSV I/O
)

options(timeout = 600)  # Allow for slow connections

# Parallel setup
num_cores <- parallel::detectCores() - 1
plan(multisession, workers = num_cores)
cat(sprintf("--- Setup: Using %d cores for parallel downloading ---\n", num_cores))

# --- 2) Directory Structure --------------------------------------------------

script_name  <- "01a_download_ITRDB_RWL_global"
base_out_dir <- here::here("output", script_name)

dir_data    <- path(base_out_dir, "data")
dir_rwl     <- path(dir_data, "rwl")
dir_tables  <- path(base_out_dir, "tables")
dir_figures <- path(base_out_dir, "figures")

dir_create(c(dir_data, dir_rwl, dir_tables, dir_figures))

cat("Directory Structure Created:\n")
cat("  Root:   ", base_out_dir, "\n")
cat("  Data:   ", dir_data, "\n")
cat("  Tables: ", dir_tables, "\n\n")

# --- 3) Global ITRDB Metadata (Reference) ------------------------------------

meta_url  <- "https://paleo.ltrr.arizona.edu/treering/ITRDBmetadata2018Dec17.txt"
meta_path <- path(dir_data, "ITRDBmetadata2018Dec17.txt")

if (!file_exists(meta_path)) {
  cat("Downloading global metadata text file...\n")
  tryCatch(
    download.file(meta_url, destfile = meta_path, mode = "wb", quiet = TRUE),
    error = function(e) warning("Could not download metadata: ", e$message)
  )
} else {
  cat("Metadata text file already exists.\n")
}

# --- 4) Scrape File Lists (Map URLs) -----------------------------------------

base_measure_url <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/"

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

# Function: list .rwl files for a region (Tucson + NOAA template)
get_region_file_list <- function(region) {
  url <- paste0(base_measure_url, region, "/")
  cat("  Scanning: ", region, " ...\n", sep = "")
  
  html_lines <- try(readLines(url, warn = FALSE), silent = TRUE)
  if (inherits(html_lines, "try-error")) {
    warning("Could not read directory: ", region)
    return(NULL)
  }
  
  # Extract all .rwl hrefs (both tucson.rwl and *-noaa.rwl)
  rwl_lines <- grep('href=\"[^\"]+\\.rwl\"', html_lines, value = TRUE)
  if (length(rwl_lines) == 0) return(NULL)
  
  filenames <- sub('.*href=\"([^\"]+\\.rwl)\".*', "\\1", rwl_lines)
  filenames <- unique(filenames)
  filenames <- filenames[!is.na(filenames) & filenames != ""]
  
  if (length(filenames) == 0) return(NULL)
  
  # Classify file type by name
  file_type <- dplyr::case_when(
    stringr::str_detect(filenames, "(?i)-noaa\\.rwl$") ~ "noaa_rwl",   # NOAA template
    TRUE                                               ~ "tucson_rwl" # classic .rwl
  )
  
  data.frame(
    region    = region,
    filename  = filenames,
    file_type = file_type,
    url       = paste0(url, filenames),
    dest_path = path(dir_rwl, region, filenames),
    stringsAsFactors = FALSE
  )
}

cat("--- Step 4: Scanning NOAA directories for .rwl files (Tucson + NOAA) ---\n")

all_files_df <- purrr::map_dfr(regions, get_region_file_list)

cat(sprintf("  Found %d total .rwl files across %d regions.\n", 
            nrow(all_files_df), length(regions)))
cat("  By file type:\n")
print(table(all_files_df$file_type))
cat("\n")

# Save master index
write_csv(all_files_df, path(dir_tables, "01a_ITRDB_master_index.csv"))

# --- 5) Parallel Download ----------------------------------------------------

cat("--- Step 5: Starting parallel download of .rwl files ---\n")

# Ensure region subfolders exist
unique_dirs <- unique(dirname(all_files_df$dest_path))
dir_create(unique_dirs)

# Safe download wrapper
download_safe <- function(url, dest) {
  if (file_exists(dest)) {
    return("skipped")
  }
  tryCatch({
    utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    "downloaded"
  }, error = function(e) {
    "failed"
  })
}

download_results <- future_map2_chr(
  .x = all_files_df$url,
  .y = all_files_df$dest_path,
  .f = download_safe,
  .options = furrr_options(seed = TRUE)
)

all_files_df$status <- download_results

# --- 6) Summary & Save Indices -----------------------------------------------

cat("\n--- Download Summary ---\n")
print(table(all_files_df$status))
cat("\nBy file type and status:\n")
print(table(all_files_df$file_type, all_files_df$status))

final_index_path <- path(dir_tables, "01a_ITRDB_download_status.csv")
write_csv(all_files_df, final_index_path)

cat("\nProcessing complete.\n")
cat("Output saved to: ", base_out_dir, "\n", sep = "")
