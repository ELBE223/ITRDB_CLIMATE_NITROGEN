############################################################
# 02_build_BAI_base_global.R
# Build global site-level metadata + BAI-ready info
# Uses:
#   - itrdb_all_regions_rwl_index.csv (from 01)
#   - NOAA template headers (*-rwl-noaa.txt) online
#   - classic *.rwl files locally for dplR / BAI
############################################################

# --- 0) Packages ---------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(here)  # project root helper

# --- 1) Paths ------------------------------------------------------

base_dir   <- here::here("itrdb_global_example")  # same as in 01
rwl_root   <- file.path(base_dir, "rwl")
meta_dir   <- file.path(base_dir, "metadata")
index_path <- file.path(base_dir, "itrdb_all_regions_rwl_index.csv")

if (!file.exists(index_path)) {
  stop("Global index file not found: ", index_path,
       "\nRun 01_download_ITRDB_RWL_global.R first.")
}

dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2) Read global index of all RWL files ------------------------

index_df <- read.csv(index_path, stringsAsFactors = FALSE)

# Required columns from script 01
required_cols <- c("region", "rwl_file")
if (!all(required_cols %in% names(index_df))) {
  stop("Index file must contain columns: ",
       paste(required_cols, collapse = ", "))
}

# Keep only NOAA template RWL files (*-noaa.rwl)
noaa_idx <- grepl("-noaa\\.rwl$", index_df$rwl_file, ignore.case = TRUE)
noaa_df  <- index_df[noaa_idx, , drop = FALSE]

if (nrow(noaa_df) == 0L) {
  stop("No '*-noaa.rwl' files found in global index. Check script 01.")
}

# Base codes without '-noaa.rwl' (e.g. 'ausl039-noaa.rwl' -> 'ausl039')
base_codes <- tolower(sub("-noaa\\.rwl$", "", noaa_df$rwl_file))
regions    <- noaa_df$region  # e.g. "australia", "northamerica/usa"

# Only ring-width base files: letters(2–4) + digits(3)
# Includes US codes like 'ak001', excludes 'newz076e' etc.
ring_mask  <- grepl("^[a-z]{2,4}[0-9]{3}$", base_codes)
noaa_df    <- noaa_df[ring_mask, , drop = FALSE]
base_codes <- base_codes[ring_mask]
regions    <- regions[ring_mask]

if (nrow(noaa_df) == 0L) {
  stop("No ring-width '*-noaa.rwl' files found (pattern [a-z]{2,4}[0-9]{3}-noaa.rwl).")
}

# --- 3) Helper: fetch NOAA template header for one site -----------

base_measure_url <- "https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/"

get_noaa_header_lines <- function(region_path, site_base) {
  # region_path e.g. "australia", "northamerica/usa"
  # site_base   e.g. "ausl039"
  url_txt <- paste0(base_measure_url, region_path, "/", site_base, "-rwl-noaa.txt")
  message("Fetching NOAA template header: ", url_txt)
  
  out <- try(readLines(url_txt, warn = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) {
    warning("Could not read NOAA template for ", site_base,
            " in region ", region_path,
            " (URL: ", url_txt, ")")
    return(NULL)
  }
  out
}

# --- 4) Helper: parse metadata from NOAA header -------------------

parse_noaa_metadata <- function(lines, region_path, site_base, rwl_noaa_file) {
  # Keep only header/comment lines
  h <- lines[grepl("^#", lines)]
  
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
  
  # Simple averages (most sites are points with N=S, E=W)
  lat <- if (all(is.na(c(north, south)))) NA_real_ else
    mean(c(north, south), na.rm = TRUE)
  lon <- if (all(is.na(c(east, west)))) NA_real_ else
    mean(c(east, west), na.rm = TRUE)
  
  elev  <- suppressWarnings(as.numeric(get_field("Elevation_m")))
  fyear <- suppressWarnings(as.integer(get_field("First_Year")))
  lyear <- suppressWarnings(as.integer(get_field("Last_Year")))
  
  site_code <- toupper(get_field("Collection_Name"))
  if (is.na(site_code) || site_code == "") {
    site_code <- toupper(site_base)
  }
  
  data.frame(
    region           = region_path,
    site_code        = site_code,
    site_base        = toupper(site_base),
    site_name        = get_field("Site_Name"),
    location         = get_field("Location"),
    lat_north        = north,
    lat_south        = south,
    lon_east         = east,
    lon_west         = west,
    lat              = lat,
    lon              = lon,
    elevation_m      = elev,
    first_year       = fyear,
    last_year        = lyear,
    species_name     = get_field("Species_Name"),
    common_name      = get_field("Common_Name"),
    species_code     = get_field("Tree_Species_Code"),
    parameter        = get_field("Parameter_Keywords"),
    rwl_noaa_file    = rwl_noaa_file,
    stringsAsFactors = FALSE
  )
}

# --- 5) Loop over all ring-width NOAA files ----------------------

meta_list <- vector("list", length(base_codes))
keep_flag <- rep(FALSE, length(base_codes))

for (i in seq_along(base_codes)) {
  base_i        <- base_codes[i]
  region_i      <- regions[i]
  rwl_noaa_file <- noaa_df$rwl_file[i]
  
  header_lines <- get_noaa_header_lines(region_i, base_i)
  if (is.null(header_lines)) next
  
  meta_list[[i]] <- parse_noaa_metadata(
    lines         = header_lines,
    region_path   = region_i,
    site_base     = base_i,
    rwl_noaa_file = rwl_noaa_file
  )
  keep_flag[i] <- TRUE
}

if (!any(keep_flag)) {
  stop("No metadata could be parsed from NOAA template headers.")
}

# Combine successfully parsed sites
meta_df <- do.call(rbind, meta_list[keep_flag])

if (is.null(meta_df) || nrow(meta_df) == 0L) {
  stop("No metadata could be parsed from NOAA template headers.")
}

# --- 6) Add local classic RWL file info --------------------------

# 'classic' Tucson-format RWL name from NOAA file name:
# e.g. 'ausl039-noaa.rwl' -> 'ausl039.rwl'
meta_df$classic_rwl_file <- sub("-noaa\\.rwl$", ".rwl", meta_df$rwl_noaa_file)

# Full local path: rwl_root / region / file
meta_df$classic_rwl_path <- file.path(
  rwl_root,
  meta_df$region,
  meta_df$classic_rwl_file
)

meta_df$has_classic_rwl <- file.exists(meta_df$classic_rwl_path)

# Full URL for classic RWL
meta_df$classic_rwl_url <- paste0(
  base_measure_url,
  meta_df$region,
  "/",
  meta_df$classic_rwl_file
)

# --- 7) Save global site-level metadata table --------------------

site_meta_path <- file.path(meta_dir, "itrdb_global_site_metadata.csv")

write.csv(
  meta_df,
  site_meta_path,
  row.names = FALSE
)

cat("Global site-level metadata written to:\n",
    normalizePath(site_meta_path), "\n")
cat("Columns include (among others):\n",
    "  region, site_code, site_name, lat, lon, elevation_m,\n",
    "  species_code, species_name, first_year, last_year,\n",
    "  rwl_noaa_file, classic_rwl_file, classic_rwl_path, has_classic_rwl\n\n")
