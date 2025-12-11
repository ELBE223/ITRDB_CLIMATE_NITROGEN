###############################################################################
# 01a_analyse_ITRDB.R
#
# Description:
#   Parse site metadata (Lat, Lon, Species) from NOAA text headers.
#   Link metadata to BOTH local RWL file types from step 00a:
#     - tucson_rwl   (e.g., alge001.rwl)
#     - noaa_rwl     (e.g., alge001-noaa.rwl)
#
# Input:
#   - output/00a_download_ITRDB_RWL_global/tables/00a_ITRDB_download_status.csv
#   - output/00a_download_ITRDB_RWL_global/data/rwl/
#
# Output:
#   - output/01a_analyse_ITRDB/tables/01a_ITRDB_site_metadata.csv
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,     
  fs,        
  dplyr,     
  stringr,   
  readr,     
  furrr,     
  tidyr      
)

# Parallel setup
num_cores <- parallel::detectCores() - 1
plan(multisession, workers = num_cores)

# --- 2) Directories & Paths --------------------------------------------------

# Input from 00a
input_script_name <- "00a_download_ITRDB_RWL_global"
input_dir         <- here::here("output", input_script_name)
input_index_path  <- path(input_dir, "tables", "00a_ITRDB_download_status.csv")
input_rwl_root    <- path(input_dir, "data", "rwl")

# Output for 01a
current_script_name <- "01a_analyse_ITRDB"
base_out_dir        <- here::here("output", current_script_name)

dir_tables  <- path(base_out_dir, "tables")
dir_data    <- path(base_out_dir, "data")
dir_figures <- path(base_out_dir, "figures")

dir_create(c(dir_tables, dir_data, dir_figures))

cat("--- Setup ---\n")
cat("Reading input from: ", input_dir, "\n")
cat("Saving output to:   ", base_out_dir, "\n")
cat(sprintf("Using %d cores for parallel processing.\n\n", num_cores))

# --- 3) Read & Prepare Input Index -------------------------------------------

if (!file_exists(input_index_path)) {
  stop("Input index not found. Please run 00a_download_ITRDB_RWL_global.R first.")
}

files_df <- read_csv(input_index_path, show_col_types = FALSE)

# Derive site_base for both file types
files_df <- files_df %>%
  mutate(
    site_base = case_when(
      file_type == "noaa_rwl"   ~ str_to_lower(str_remove(filename, "(?i)-noaa\\.rwl$")),
      file_type == "tucson_rwl" ~ str_to_lower(str_remove(filename, "(?i)\\.rwl$")),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(site_base),
    str_detect(site_base, "^[a-z]{2,4}[0-9]{3}$")  # standard ITRDB site codes
  )

# Build per-site file mapping (Tucson + NOAA)
site_files_df <- files_df %>%
  group_by(region, site_base) %>%
  summarise(
    rwl_noaa_file   = filename[file_type == "noaa_rwl"][1],
    rwl_tucson_file = filename[file_type == "tucson_rwl"][1],
    .groups = "drop"
  )

# Target for metadata fetch: only sites with a NOAA template
target_df <- site_files_df %>%
  filter(!is.na(rwl_noaa_file))

cat(sprintf("Found %d sites with NOAA template files for metadata.\n", nrow(target_df)))

# --- 4) Helper for NOAA Metadata ---------------------------------------------

base_measure_url <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/"

# Fetch and parse metadata for ONE site
process_site_metadata <- function(region, site_base, rwl_noaa_file) {
  
  # Build metadata text URL (NOAA convention)
  txt_url <- paste0(base_measure_url, region, "/", site_base, "-rwl-noaa.txt")
  
  # Read header lines
  lines <- tryCatch(
    readLines(txt_url, warn = FALSE),
    error = function(e) return(NULL)
  )
  if (is.null(lines)) return(NULL)
  
  h <- lines[grep("^#", lines)]  # header lines only
  
  get_field <- function(label) {
    m <- grep(paste0(label, "\\s*:"), h, value = TRUE)
    if (!length(m)) return(NA_character_)
    val <- sub(".*:\\s*", "", m[1])
    trimws(val)
  }
  
  # Coordinates
  north <- suppressWarnings(as.numeric(get_field("Northernmost_Latitude")))
  south <- suppressWarnings(as.numeric(get_field("Southernmost_Latitude")))
  east  <- suppressWarnings(as.numeric(get_field("Easternmost_Longitude")))
  west  <- suppressWarnings(as.numeric(get_field("Westernmost_Longitude")))
  
  lat <- if (all(is.na(c(north, south)))) NA_real_ else mean(c(north, south), na.rm = TRUE)
  lon <- if (all(is.na(c(east, west))))  NA_real_ else mean(c(east, west), na.rm = TRUE)
  
  # Other fields
  elev         <- suppressWarnings(as.numeric(get_field("Elevation_m")))
  fyear        <- suppressWarnings(as.integer(get_field("First_Year")))
  lyear        <- suppressWarnings(as.integer(get_field("Last_Year")))
  site_name    <- get_field("Site_Name")
  species_name <- get_field("Species_Name")
  species_code <- get_field("Tree_Species_Code")
  
  # Site code from Collection_Name, fallback to SITE_BASE
  site_code_extracted <- toupper(get_field("Collection_Name"))
  if (is.na(site_code_extracted) || site_code_extracted == "") {
    site_code_extracted <- toupper(site_base)
  }
  
  tibble(
    region        = region,
    site_code     = site_code_extracted,
    site_base     = site_base,          # lower-case internal key
    site_name     = site_name,
    lat           = lat,
    lon           = lon,
    elevation_m   = elev,
    first_year    = fyear,
    last_year     = lyear,
    species_name  = species_name,
    species_code  = species_code,
    rwl_noaa_file = rwl_noaa_file,
    meta_url      = txt_url
  )
}

# --- 5) Fetch Metadata in Parallel -------------------------------------------

cat("\n--- Fetching metadata from NOAA (parallel) ---\n")
cat("This may take a few minutes...\n")

meta_results <- target_df %>%
  select(region, site_base, rwl_noaa_file) %>%
  future_pmap_dfr(
    process_site_metadata,
    .progress = TRUE,
    .options  = furrr_options(seed = TRUE)
  )

# --- 6) Link Metadata with Local Files & Save --------------------------------

cat("\nLinking metadata with local file info...\n")

final_meta_df <- meta_results %>%
  # add tucson_rwl information
  left_join(
    site_files_df %>% select(region, site_base, rwl_tucson_file),
    by = c("region", "site_base")
  ) %>%
  # build full local paths as character
  mutate(
    full_rwl_noaa_path = ifelse(
      is.na(rwl_noaa_file),
      NA_character_,
      as.character(path(input_rwl_root, region, rwl_noaa_file))
    ),
    full_rwl_tucson_path = ifelse(
      is.na(rwl_tucson_file),
      NA_character_,
      as.character(path(input_rwl_root, region, rwl_tucson_file))
    )
  ) %>%
  mutate(
    has_rwl_noaa_local   = !is.na(full_rwl_noaa_path)   & file_exists(full_rwl_noaa_path),
    has_rwl_tucson_local = !is.na(full_rwl_tucson_path) & file_exists(full_rwl_tucson_path)
  ) %>%
  # optional: store site_base also in upper-case for readability
  mutate(site_base_upper = toupper(site_base)) %>%
  relocate(site_base_upper, .after = site_base)

# Save to CSV
output_file <- path(dir_tables, "01a_ITRDB_site_metadata.csv")
write_csv(final_meta_df, output_file)

cat("\n--- Complete ---\n")
cat("Successfully parsed metadata for ", nrow(final_meta_df), " sites.\n", sep = "")
cat("Metadata saved to:\n  ", output_file, "\n", sep = "")
