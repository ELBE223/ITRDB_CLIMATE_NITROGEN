###############################################################################
# 00a_download_ITRDB_RWL_global.R  
#
# Description:
#   Download ITRDB ring-width (.rwl) files by region using parallel processing.
#   Downloads:
#     - classic Tucson .rwl files (e.g., alge001.rwl)
#     - NOAA template .rwl files (e.g., alge001-noaa.rwl)
#     - NOAA metadata TXT files (e.g., alge001-rwl-noaa.txt)  <-- NEW
#
# Output structure:
#   output/00a_download_ITRDB_RWL_global/
#       ├── data/
#       │    ├── rwl/       (Raw .rwl files by region)
#       │    ├── meta_txt/  (NOAA *-rwl-noaa.txt by region)   <-- NEW
#       │    └── ITRDBmetadata2018Dec17.txt (global reference)
#       ├── tables/
#       └── figures/
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, fs, dplyr, stringr, furrr, future, readr, purrr
)

CFG <- list(
  TIMEOUT_SEC = 600,
  PARALLEL = TRUE,
  MAX_WORKERS = 4,        # keep this conservative for web servers
  RETRIES = 3,
  BACKOFF_SEC = 2
)

options(timeout = CFG$TIMEOUT_SEC)

num_cores <- max(1, parallel::detectCores() - 1)
workers <- if (CFG$PARALLEL) min(num_cores, CFG$MAX_WORKERS) else 1
future::plan(future::multisession, workers = workers)
cat(sprintf("--- Setup: Using %d workers for downloading ---\n", workers))

# --- 2) Directory Structure --------------------------------------------------

script_name  <- "00a_download_ITRDB_RWL_global"
base_out_dir <- here::here("output", script_name)

dir_data     <- fs::path(base_out_dir, "data")
dir_rwl      <- fs::path(dir_data, "rwl")
dir_meta_txt <- fs::path(dir_data, "meta_txt")   # NEW
dir_tables   <- fs::path(base_out_dir, "tables")
dir_figures  <- fs::path(base_out_dir, "figures")

fs::dir_create(c(dir_data, dir_rwl, dir_meta_txt, dir_tables, dir_figures))

cat("Directory Structure Created:\n")
cat("  Root:     ", base_out_dir, "\n")
cat("  RWL:      ", dir_rwl, "\n")
cat("  META TXT: ", dir_meta_txt, "\n")
cat("  Tables:   ", dir_tables, "\n\n")

# --- 3) Global ITRDB Metadata (Reference) ------------------------------------

meta_url  <- "https://paleo.ltrr.arizona.edu/treering/ITRDBmetadata2018Dec17.txt"
meta_path <- fs::path(dir_data, "ITRDBmetadata2018Dec17.txt")

if (!fs::file_exists(meta_path)) {
  cat("Downloading global metadata text file...\n")
  tryCatch(
    utils::download.file(meta_url, destfile = meta_path, mode = "wb", quiet = TRUE),
    error = function(e) warning("Could not download metadata: ", e$message)
  )
} else {
  cat("Global metadata text file already exists.\n")
}

# --- 4) Scrape File Lists (Map URLs) -----------------------------------------

base_measure_url <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/"

regions <- c(
  "asia",
  "europe",
  "northamerica/canada",
  "northamerica/usa"
)

# helper: extract href targets ending with pattern
extract_hrefs <- function(html_lines, pattern_regex) {
  x <- grep(pattern_regex, html_lines, value = TRUE)
  if (length(x) == 0) return(character())
  f <- sub('.*href="([^"]+)".*', "\\1", x)
  f <- unique(f)
  f <- f[!is.na(f) & f != ""]
  f
}

get_region_file_list <- function(region) {
  url <- paste0(base_measure_url, region, "/")
  cat("  Scanning: ", region, " ...\n", sep = "")
  
  html_lines <- try(readLines(url, warn = FALSE), silent = TRUE)
  if (inherits(html_lines, "try-error")) {
    warning("Could not read directory: ", region)
    return(NULL)
  }
  
  # 1) .rwl files
  rwl_files <- extract_hrefs(html_lines, 'href="[^"]+\\.rwl"')
  
  # 2) NOAA per-site metadata TXT files (only the ones we need)
  #    example: alge001-rwl-noaa.txt
  txt_files <- extract_hrefs(html_lines, 'href="[^"]+-rwl-noaa\\.txt"')
  
  if (length(rwl_files) == 0 && length(txt_files) == 0) return(NULL)
  
  # build df for rwl
  df_rwl <- tibble()
  if (length(rwl_files) > 0) {
    rwl_type <- dplyr::case_when(
      stringr::str_detect(rwl_files, "(?i)-noaa\\.rwl$") ~ "noaa_rwl",
      TRUE                                              ~ "tucson_rwl"
    )
    df_rwl <- tibble(
      region = region,
      filename = rwl_files,
      file_type = rwl_type,
      url = paste0(url, rwl_files),
      dest_path = fs::path(dir_rwl, region, rwl_files)
    )
  }
  
  # build df for txt
  df_txt <- tibble()
  if (length(txt_files) > 0) {
    df_txt <- tibble(
      region = region,
      filename = txt_files,
      file_type = "noaa_meta_txt",
      url = paste0(url, txt_files),
      dest_path = fs::path(dir_meta_txt, region, txt_files)
    )
  }
  
  bind_rows(df_rwl, df_txt)
}

cat("--- Step 4: Scanning NOAA directories for .rwl + *-rwl-noaa.txt ---\n")
all_files_df <- purrr::map_dfr(regions, get_region_file_list)

cat(sprintf("  Found %d total files across %d regions.\n",
            nrow(all_files_df), length(regions)))
cat("  By file type:\n")
print(table(all_files_df$file_type))
cat("\n")

# Derive site_base (useful for downstream joins)
all_files_df <- all_files_df %>%
  mutate(
    site_base = case_when(
      file_type == "noaa_rwl"      ~ str_to_lower(str_remove(filename, "(?i)-noaa\\.rwl$")),
      file_type == "tucson_rwl"    ~ str_to_lower(str_remove(filename, "(?i)\\.rwl$")),
      file_type == "noaa_meta_txt" ~ str_to_lower(str_remove(filename, "(?i)-rwl-noaa\\.txt$")),
      TRUE ~ NA_character_
    )
  )

# Save master index (important for step 01a / 01b)
readr::write_csv(all_files_df, fs::path(dir_tables, "00a_ITRDB_master_index.csv"))

# --- 5) Parallel Download ----------------------------------------------------

cat("--- Step 5: Starting parallel download ---\n")

# Ensure region subfolders exist locally (both for rwl and meta_txt)
unique_dirs <- unique(dirname(all_files_df$dest_path))
fs::dir_create(unique_dirs)

download_safe <- function(url, dest, retries = CFG$RETRIES, backoff = CFG$BACKOFF_SEC) {
  if (fs::file_exists(dest)) return("skipped")
  
  for (i in seq_len(retries)) {
    ok <- tryCatch({
      utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    
    if (ok && fs::file_exists(dest) && fs::file_size(dest) > 0) return("downloaded")
    
    # backoff (simple)
    Sys.sleep(backoff * i)
  }
  
  "failed"
}

download_results <- furrr::future_map2_chr(
  .x = all_files_df$url,
  .y = all_files_df$dest_path,
  .f = ~ download_safe(.x, .y),
  .options = furrr::furrr_options(seed = TRUE)
)

all_files_df$status <- download_results

# --- 6) Summary & Save Indices -----------------------------------------------

cat("\n--- Download Summary ---\n")
print(table(all_files_df$status))

cat("\nBy file type and status:\n")
print(table(all_files_df$file_type, all_files_df$status))

# Save status log
final_index_path <- fs::path(dir_tables, "00a_ITRDB_download_status.csv")
readr::write_csv(all_files_df, final_index_path)

cat("\n--- Sanity check (local counts) ---\n")
cat("RWL files local:      ", sum(all_files_df$file_type %in% c("tucson_rwl", "noaa_rwl") &
                                    all_files_df$status %in% c("downloaded", "skipped")), "\n", sep = "")
cat("META TXT files local: ", sum(all_files_df$file_type == "noaa_meta_txt" &
                                    all_files_df$status %in% c("downloaded", "skipped")), "\n", sep = "")
cat("\nProcessing complete.\n")
cat("Output saved to: ", base_out_dir, "\n", sep = "")
cat("Next: run 01a (now parsing from data/meta_txt/ ...).\n")
