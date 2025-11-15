############################################################
# 05_extract_BAI_boreal_global.R
#
# Goals:
#  - Restrict ITRDB global BAI data to the boreal zone
#    (lat >= 50°N)
#  - Inspect species composition and number of data points
#  - Save boreal BAI data as Parquet for later LMM analysis
#  - Summarise how species and leaf-type composition change
#    with increasing latitude cutoff (50,55,60,65,70 °N)
#
# Uses:
#  - metadata/itrdb_global_site_metadata.csv
#  - metadata/itrdb_global_BAI_site_summary.csv
#  - bai/<region>/<SITE_CODE>_BAI.csv   (from script 03)
############################################################

## --- 0) Packages ---------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  here,
  dplyr,
  tidyr,
  arrow
)

## --- 1) Paths ------------------------------------------------------

base_dir <- here::here("itrdb_global_example")

meta_dir <- file.path(base_dir, "metadata")
bai_root <- file.path(base_dir, "bai")

site_meta_path <- file.path(meta_dir, "itrdb_global_site_metadata.csv")
bai_sum_path   <- file.path(meta_dir, "itrdb_global_BAI_site_summary.csv")

if (!file.exists(site_meta_path)) {
  stop("Site metadata not found: ", site_meta_path,
       "\nRun 02_build_BAI_base_global.R first.")
}
if (!file.exists(bai_sum_path)) {
  stop("BAI site summary not found: ", bai_sum_path,
       "\nRun 03_compute_BAI_all_sites_global.R first.")
}

# Folder for prepared data (for modelling & summaries)
prep_dir <- file.path(base_dir, "data_prepared")
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

## --- 2) Load metadata and BAI summary ------------------------------

site_meta_all <- read.csv(site_meta_path, stringsAsFactors = FALSE)
bai_sum       <- read.csv(bai_sum_path,   stringsAsFactors = FALSE)

# required columns in metadata
req_cols_meta <- c(
  "region", "site_code", "site_name",
  "lat", "lon", "elevation_m",
  "species_code", "species_name",
  "classic_rwl_path", "has_classic_rwl",
  "first_year", "last_year"
)
missing_meta <- setdiff(req_cols_meta, names(site_meta_all))
if (length(missing_meta) > 0) {
  stop("Missing columns in site metadata: ",
       paste(missing_meta, collapse = ", "))
}

# required columns in BAI summary
req_cols_sum <- c(
  "region", "site_code",
  "bai_mean", "bai_min", "bai_max",
  "n_series"
)
missing_sum <- setdiff(req_cols_sum, names(bai_sum))
if (length(missing_sum) > 0) {
  stop("Missing columns in BAI summary: ",
       paste(missing_sum, collapse = ", "))
}

## --- 3) Restrict to sites with coords + RWL -----------------------

site_meta <- site_meta_all %>%
  dplyr::filter(!is.na(lat), !is.na(lon), has_classic_rwl)

## --- 4) Boreal subset (lat >= 50°N) -------------------------------

boreal_sites <- site_meta %>%
  dplyr::filter(lat >= 50)

cat("Total sites with RWL and coordinates:", nrow(site_meta),  "\n")
cat("Boreal sites (lat >= 50°N):          ", nrow(boreal_sites), "\n\n")

if (nrow(boreal_sites) == 0L) {
  stop("No boreal sites found (lat >= 50).")
}

## --- 5) Species composition in boreal zone ------------------------

species_counts <- boreal_sites %>%
  dplyr::count(species_code, species_name, sort = TRUE)

cat("Number of species in boreal subset:", nrow(species_counts), "\n\n")
cat("Top 20 species by number of sites:\n")
print(head(species_counts, 20))
cat("\n")

# save boreal species summary (for inspection / later use)
species_summary_path <- file.path(
  prep_dir,
  "itrdb_boreal_species_summary.csv"
)
write.csv(species_counts, species_summary_path, row.names = FALSE)

cat("Boreal species summary saved to:\n  ",
    normalizePath(species_summary_path), "\n\n")

## --- 6) Join BAI summary to boreal sites --------------------------

bai_sum_small <- bai_sum %>%
  dplyr::select(region, site_code,
                bai_mean, bai_min, bai_max,
                n_series)

boreal_sites_full <- boreal_sites %>%
  dplyr::left_join(bai_sum_small,
                   by = c("region", "site_code"))

## --- 7) Build boreal BAI long table (site-year-series) ------------

boreal_long_list <- vector("list", nrow(boreal_sites_full))

for (i in seq_len(nrow(boreal_sites_full))) {
  row <- boreal_sites_full[i, ]
  
  bai_file <- file.path(
    bai_root,
    row$region,
    paste0(row$site_code, "_BAI.csv")
  )
  
  if (!file.exists(bai_file)) {
    message("Skipping (no BAI file): ", bai_file)
    next
  }
  
  bai_df <- try(
    read.csv(bai_file, stringsAsFactors = FALSE),
    silent = TRUE
  )
  
  if (inherits(bai_df, "try-error") ||
      !"year" %in% names(bai_df)) {
    message("Skipping (could not read BAI): ", bai_file)
    next
  }
  
  series_cols <- setdiff(names(bai_df), "year")
  if (!length(series_cols)) {
    message("Skipping (no series columns): ", bai_file)
    next
  }
  
  # long format: year x series
  bai_long <- bai_df %>%
    tidyr::pivot_longer(
      cols      = all_of(series_cols),
      names_to  = "series_id",
      values_to = "BAI_mm2"
    )
  
  # enforce types and drop NA
  bai_long$year    <- as.integer(bai_long$year)
  bai_long$BAI_mm2 <- as.numeric(bai_long$BAI_mm2)
  bai_long <- bai_long[!is.na(bai_long$BAI_mm2), , drop = FALSE]
  
  # attach site-level info
  bai_long$region        <- row$region
  bai_long$site_code     <- row$site_code
  bai_long$site_name     <- row$site_name
  bai_long$lat           <- row$lat
  bai_long$lon           <- row$lon
  bai_long$elevation_m   <- row$elevation_m
  bai_long$species_code  <- row$species_code
  bai_long$species_name  <- row$species_name
  bai_long$first_year    <- row$first_year
  bai_long$last_year     <- row$last_year
  bai_long$n_series_site <- row$n_series
  
  boreal_long_list[[i]] <- bai_long
}

boreal_long_list <- boreal_long_list[!vapply(boreal_long_list, is.null, logical(1))]

if (!length(boreal_long_list)) {
  stop("No boreal BAI time series could be loaded.")
}

boreal_bai_long <- dplyr::bind_rows(boreal_long_list)

cat("Total BAI data points in boreal subset (rows in long table): ",
    nrow(boreal_bai_long), "\n\n")

## --- 8) Save boreal data as Parquet -------------------------------

# Site-level metadata (one row per site)
boreal_sites_parquet_path <- file.path(
  prep_dir,
  "itrdb_boreal_site_metadata.parquet"
)

arrow::write_parquet(
  boreal_sites_full,
  boreal_sites_parquet_path
)

# Long BAI dataset (site-year-series)
boreal_bai_parquet_path <- file.path(
  prep_dir,
  "itrdb_boreal_BAI_long.parquet"
)

arrow::write_parquet(
  boreal_bai_long,
  boreal_bai_parquet_path
)

cat("Boreal site metadata saved to:\n  ",
    normalizePath(boreal_sites_parquet_path), "\n")
cat("Boreal BAI long data saved to:\n  ",
    normalizePath(boreal_bai_parquet_path), "\n\n")

## --- 9) Leaf-type mapping (conifer vs broadleaf) ------------------

# build leaf-type mapping from boreal species list
species_leaf <- species_counts %>%
  dplyr::mutate(
    genus = sub(" .*", "", species_name),
    leaf_type = dplyr::case_when(
      genus %in% c("Pinus", "Picea", "Tsuga",
                   "Abies", "Larix", "Pseudotsuga",
                   "Chamaecyparis", "Thuja") ~ "conifer",
      genus %in% c("Quercus", "Fagus", "Populus",
                   "Betula", "Salix", "Fraxinus",
                   "Alnus", "Tilia") ~ "broadleaf",
      TRUE ~ "unknown"
    )
  ) %>%
  dplyr::select(species_code, species_name, genus, leaf_type)

species_leaf_path <- file.path(
  prep_dir,
  "itrdb_boreal_species_leaf_type.csv"
)
write.csv(species_leaf, species_leaf_path, row.names = FALSE)

cat("Species leaf-type mapping saved to:\n  ",
    normalizePath(species_leaf_path), "\n\n")

## --- 10) Latitude gradient summary (50,55,60,65,70 °N) ------------

# base set: all sites with coords + RWL, lat >= 50
site_lat <- site_meta %>%
  dplyr::filter(lat >= 50) %>%
  dplyr::left_join(
    species_leaf %>% dplyr::select(species_code, leaf_type),
    by = "species_code"
  )

cat("Sites with lat >= 50°N and RWL (for gradient):",
    nrow(site_lat), "\n\n")

lat_thresholds <- c(50, 55, 60, 65, 70)

lat_stats_list <- lapply(lat_thresholds, function(th) {
  dat <- dplyr::filter(site_lat, lat >= th)
  if (nrow(dat) == 0L) return(NULL)
  
  summ <- dat %>%
    dplyr::summarise(
      lat_min              = th,
      n_sites              = dplyr::n(),
      n_sites_with_type    = sum(!is.na(leaf_type)),
      n_species_total      = dplyr::n_distinct(species_code),
      n_species_conifer    = dplyr::n_distinct(species_code[leaf_type == "conifer"]),
      n_species_broadleaf  = dplyr::n_distinct(species_code[leaf_type == "broadleaf"]),
      n_species_unknown    = dplyr::n_distinct(species_code[is.na(leaf_type) | leaf_type == "unknown"]),
      n_conifer_sites      = sum(leaf_type == "conifer",   na.rm = TRUE),
      n_broadleaf_sites    = sum(leaf_type == "broadleaf", na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      n_unknown_sites        = n_sites - n_sites_with_type,
      prop_conifer_species   = ifelse(n_species_total > 0,
                                      n_species_conifer / n_species_total,
                                      NA_real_),
      prop_broadleaf_species = ifelse(n_species_total > 0,
                                      n_species_broadleaf / n_species_total,
                                      NA_real_),
      prop_conifer_sites     = ifelse(n_sites_with_type > 0,
                                      n_conifer_sites / n_sites_with_type,
                                      NA_real_),
      prop_broadleaf_sites   = ifelse(n_sites_with_type > 0,
                                      n_broadleaf_sites / n_sites_with_type,
                                      NA_real_)
    )
  
  summ
})

lat_stats <- dplyr::bind_rows(lat_stats_list)

lat_stats_path <- file.path(
  prep_dir,
  "itrdb_boreal_latitude_species_gradient.csv"
)
write.csv(lat_stats, lat_stats_path, row.names = FALSE)

cat("Latitude gradient table saved to:\n  ",
    normalizePath(lat_stats_path), "\n\n")

print(lat_stats)

############################################################
# Result:
#  - Boreal site metadata:
#      data_prepared/itrdb_boreal_site_metadata.parquet
#  - Boreal BAI long table (for LMM etc.):
#      data_prepared/itrdb_boreal_BAI_long.parquet
#  - Boreal species summary:
#      data_prepared/itrdb_boreal_species_summary.csv
#  - Species leaf-type mapping:
#      data_prepared/itrdb_boreal_species_leaf_type.csv
#  - Latitude gradient (50,55,60,65,70 °N):
#      data_prepared/itrdb_boreal_latitude_species_gradient.csv
############################################################
