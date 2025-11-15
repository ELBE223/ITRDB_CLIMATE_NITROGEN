############################################################
# 04_QC_BAI_global.R
# Quality check for global BAI:
#  - recompute BAI from RWL
#  - compare with stored BAI CSVs
#  - flag extreme mean BAI sites
#  - create a world map, colored by continent group
#
# Uses:
#   - metadata/itrdb_global_site_metadata.csv
#   - metadata/itrdb_global_BAI_site_summary.csv
#   - bai/<region>/<SITE_CODE>_BAI.csv
############################################################

## --- 0) Packages ---------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  here,
  dplR,
  sf,
  ggplot2,
  rnaturalearth,
  rnaturalearthdata
)

## --- 1) Paths ------------------------------------------------------

base_dir       <- here::here("itrdb_global_example")
rwl_root       <- file.path(base_dir, "rwl")
meta_dir       <- file.path(base_dir, "metadata")
bai_root       <- file.path(base_dir, "bai")

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

dir.create(bai_root, recursive = TRUE, showWarnings = FALSE)

## --- 2) Read metadata and summary ---------------------------------

meta_df <- read.csv(site_meta_path, stringsAsFactors = FALSE)
sum_df  <- read.csv(bai_sum_path,   stringsAsFactors = FALSE)

# keep only sites with RWL
meta_df <- meta_df[meta_df$has_classic_rwl, , drop = FALSE]

# simple check: now we use 'region' instead of 'region_path'
if (!all(c("region", "site_code", "classic_rwl_path") %in% names(meta_df))) {
  stop("Metadata is missing required columns (region, site_code, classic_rwl_path).")
}

## --- 3) QC: recompute BAI and compare ------------------------------

qc_list <- list()

for (i in seq_len(nrow(meta_df))) {
  row <- meta_df[i, ]
  
  msg <- paste0("[", i, "/", nrow(meta_df), "] ",
                row$site_code, " - ", row$site_name,
                " (", row$region, ")")
  cat(msg, "\n")
  
  # Paths
  rwl_path <- row$classic_rwl_path
  bai_file <- file.path(
    bai_root,
    row$region,
    paste0(row$site_code, "_BAI.csv")
  )
  
  # Check presence
  if (!file.exists(rwl_path)) {
    cat("  Missing RWL file, skip.\n")
    next
  }
  if (!file.exists(bai_file)) {
    cat("  Missing stored BAI file, skip.\n")
    next
  }
  
  # --- Read RWL and recompute BAI ----------------------------------
  rwl <- try(read.rwl(rwl_path), silent = TRUE)
  if (inherits(rwl, "try-error") || !is.data.frame(rwl) ||
      nrow(rwl) == 0L || ncol(rwl) == 0L) {
    cat("  Could not read RWL, skip.\n")
    next
  }
  
  bai_new <- try(bai.in(rwl), silent = TRUE)
  if (inherits(bai_new, "try-error") || !is.data.frame(bai_new) ||
      nrow(bai_new) == 0L || ncol(bai_new) == 0L) {
    cat("  BAI computation failed, skip.\n")
    next
  }
  
  # --- Read stored BAI ---------------------------------------------
  bai_stored_df <- try(read.csv(bai_file, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(bai_stored_df, "try-error") ||
      !"year" %in% names(bai_stored_df)) {
    cat("  Could not read stored BAI CSV, skip.\n")
    next
  }
  
  # Years
  years_new    <- suppressWarnings(as.integer(rownames(bai_new)))
  years_stored <- suppressWarnings(as.integer(bai_stored_df$year))
  
  if (all(is.na(years_new)) || all(is.na(years_stored))) {
    cat("  Could not parse years, skip.\n")
    next
  }
  
  common_years <- intersect(years_new, years_stored)
  if (!length(common_years)) {
    cat("  No overlapping years, skip.\n")
    next
  }
  
  # Columns (series)
  stored_cols <- setdiff(colnames(bai_stored_df), "year")
  common_cols <- intersect(colnames(bai_new), stored_cols)
  if (!length(common_cols)) {
    cat("  No common series columns, skip.\n")
    next
  }
  
  # Align by years and columns
  idx_new    <- match(common_years, years_new)
  idx_stored <- match(common_years, years_stored)
  
  bai_new_sub    <- bai_new[idx_new, common_cols, drop = FALSE]
  bai_stored_sub <- bai_stored_df[idx_stored, c("year", common_cols), drop = FALSE]
  
  # matrices
  m_new    <- as.matrix(bai_new_sub)
  m_stored <- as.matrix(bai_stored_sub[, common_cols, drop = FALSE])
  
  if (!all(dim(m_new) == dim(m_stored))) {
    cat("  Dimension mismatch after alignment, skip.\n")
    next
  }
  
  # Differences
  diffs <- as.numeric(m_new - m_stored)
  diffs <- diffs[!is.na(diffs)]
  
  if (!length(diffs)) {
    cat("  No numeric differences (all NA), skip.\n")
    next
  }
  
  max_abs_diff <- max(abs(diffs), na.rm = TRUE)
  mean_bai_rec <- mean(as.numeric(m_new), na.rm = TRUE)
  
  qc_list[[length(qc_list) + 1L]] <- data.frame(
    region              = row$region,
    site_code           = row$site_code,
    site_name           = row$site_name,
    species_code        = row$species_code,
    species_name        = row$species_name,
    lat                 = row$lat,
    lon                 = row$lon,
    elevation_m         = row$elevation_m,
    n_years_common      = length(common_years),
    n_series_common     = length(common_cols),
    max_abs_diff        = max_abs_diff,
    mean_bai_recomputed = mean_bai_rec,
    stringsAsFactors    = FALSE
  )
}

# Combine
if (length(qc_list) == 0L) {
  stop("No QC entries could be computed. Check inputs and paths.")
}

qc_df <- do.call(rbind, qc_list)

qc_path <- file.path(meta_dir, "itrdb_global_BAI_QC_summary.csv")
write.csv(qc_df, qc_path, row.names = FALSE)

cat("\nQC summary written to:\n  ", normalizePath(qc_path), "\n\n")

## --- 4) Quick QC stats --------------------------------------------

tol <- 1e-6

n_ok  <- sum(qc_df$max_abs_diff < tol, na.rm = TRUE)
n_bad <- sum(qc_df$max_abs_diff >= tol, na.rm = TRUE)

cat("Sites with BAI == stored BAI (up to rounding):", n_ok, "\n")
cat("Sites with deviations >", tol, ":", n_bad, "\n\n")

if (n_bad > 0) {
  cat("Examples with larger deviations:\n")
  print(head(qc_df[order(-qc_df$max_abs_diff), ], 10))
}

## --- 5) Extreme mean BAI sites ------------------------------------

# threshold can be adjusted, e.g. 20000 mm^2
high_thresh <- 20000

extreme_df <- qc_df[qc_df$mean_bai_recomputed > high_thresh, , drop = FALSE]

extreme_path <- file.path(meta_dir, "itrdb_global_BAI_extreme_sites.csv")
write.csv(extreme_df, extreme_path, row.names = FALSE)

cat("Sites with mean BAI >", high_thresh, "mm^2:\n")
if (nrow(extreme_df) == 0L) {
  cat("  None.\n\n")
} else {
  print(extreme_df[, c("region", "site_code", "site_name",
                       "species_code", "mean_bai_recomputed")])
  cat("\nList written to:\n  ", normalizePath(extreme_path), "\n\n")
}

## --- 6) World map of sites by continent group ---------------------

# Use BAI summary for plotting (has lat/lon + mean BAI etc.)
plot_df <- sum_df

# Ensure lat/lon + region exist
if (!all(c("lat", "lon", "region") %in% names(plot_df))) {
  stop("BAI summary is missing 'lat', 'lon', or 'region' columns.")
}

plot_df <- plot_df[!is.na(plot_df$lat) & !is.na(plot_df$lon), , drop = FALSE]

if (nrow(plot_df) == 0L) {
  stop("No sites with coordinates available for mapping.")
}

# Continent grouping based on region
continent_group <- function(region) {
  if (grepl("^africa", region, ignore.case = TRUE)) {
    "Africa"
  } else if (grepl("^asia", region, ignore.case = TRUE)) {
    "Asia"
  } else if (grepl("^australia", region, ignore.case = TRUE)) {
    "Australia/Oceania"
  } else if (grepl("^europe", region, ignore.case = TRUE)) {
    "Europe"
  } else if (grepl("^southamerica", region, ignore.case = TRUE)) {
    "South America"
  } else if (grepl("^mexico", region, ignore.case = TRUE)) {
    "North America"
  } else if (grepl("^northamerica", region, ignore.case = TRUE)) {
    "North America"
  } else {
    "Other"
  }
}

plot_df$continent_group <- vapply(
  plot_df$region,
  continent_group,
  FUN.VALUE = character(1)
)

# Optional: log10 mean BAI for symbol size
if (!"bai_mean" %in% names(plot_df)) {
  plot_df$bai_mean <- NA_real_
}
plot_df$bai_mean_log10 <- log10(plot_df$bai_mean + 1)

# Convert to sf points
sites_sf <- sf::st_as_sf(
  plot_df,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# World basemap
world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# Output folder for maps
maps_dir <- file.path(base_dir, "maps_global")
dir.create(maps_dir, recursive = TRUE, showWarnings = FALSE)

# Map: sites colored by continent group
p_world <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey80") +
  geom_sf(
    data = sites_sf,
    aes(color = continent_group),
    size = 1.5,
    alpha = 0.8
  ) +
  coord_sf(expand = FALSE) +
  guides(color = guide_legend(title = "Continent group")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "ITRDB sites with BAI (global)",
    subtitle = "Colored by continent group",
    x = NULL,
    y = NULL
  )

map_file <- file.path(maps_dir, "BAI_sites_world_by_continent.png")
ggsave(map_file, p_world, width = 10, height = 6, dpi = 300)

cat("World map saved to:\n  ", normalizePath(map_file), "\n\n")

############################################################
# Outputs:
#  - QC summary:
#      metadata/itrdb_global_BAI_QC_summary.csv
#  - Extreme mean BAI sites:
#      metadata/itrdb_global_BAI_extreme_sites.csv
#  - World map:
#      maps_global/BAI_sites_world_by_continent.png
############################################################
