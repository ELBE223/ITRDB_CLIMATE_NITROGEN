###############################################################################
# 03c_explore_nitrogen_deposition.R
#
# Description:
#   Detailed exploration of nitrogen deposition patterns and effects.
#   Extracts N-deposition data at site locations and analyzes:
#   - Global N-deposition maps (NHx, NOy, Total)
#   - N-deposition gradients (latitude, elevation, climate)
#   - Relationship between N-deposition and BAI
#   - Relationship between N-deposition and temperature
#   - Regional N-deposition patterns
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  tidyr,
  arrow,
  ggplot2,
  patchwork,
  scales,
  viridis,
  sf,
  rnaturalearth,
  rnaturalearthdata,
  terra,
  stringr
)

# --- 2) Paths ----------------------------------------------------------------

step_02a <- "02a_download_climate_and_ndep"
step_03a <- "03a_explore_and_finalize_dataset"
current_name <- "03c_explore_nitrogen_deposition"

input_data_path <- here::here("output", step_03a, "data", "03a_final_dataset.parquet")
ndep_root <- here::here("output", step_02a, "data", "climate", "ISIMIP3a", "ndep", "histsoc")

base_out_dir <- here::here("output", current_name)
dir_data     <- path(base_out_dir, "data")
dir_tables   <- path(base_out_dir, "tables")
dir_figures  <- path(base_out_dir, "figures")

dir_create(c(dir_data, dir_tables, dir_figures))

cat("=== 03c Explore Nitrogen Deposition ===\n\n")
cat("Input data: ", input_data_path, "\n")
cat("N-dep data: ", ndep_root, "\n")
cat("Output:     ", base_out_dir, "\n\n")

# --- 3) Load BAI Data --------------------------------------------------------

if (!file_exists(input_data_path)) {
  stop("Input data not found. Run 03a first.")
}

cat("Loading final dataset...\n")
df <- arrow::read_parquet(input_data_path)

cat("  Observations: ", nrow(df), "\n")
cat("  Sites:        ", n_distinct(df$site_code), "\n\n")

# Get unique sites with coordinates
sites_coords <- df %>%
  select(site_code, lat, lon) %>%
  distinct()

cat("Unique sites for N-deposition extraction: ", nrow(sites_coords), "\n\n")

# --- 4) Extract N-Deposition Data --------------------------------------------

cat("=== EXTRACTING N-DEPOSITION DATA ===\n\n")

# Helper function to extract N-deposition
extract_ndep_at_sites <- function(var_name, sites_df) {
  
  var_dir <- path(ndep_root, var_name)
  
  if (!dir_exists(var_dir)) {
    warning("N-deposition directory not found: ", var_dir)
    return(NULL)
  }
  
  nc_files <- dir_ls(var_dir, glob = "*.nc")
  
  if (length(nc_files) == 0) {
    warning("No N-deposition files found for: ", var_name)
    return(NULL)
  }
  
  nc_files <- sort(nc_files)
  
  cat("  Variable: ", var_name, "\n")
  for (f in nc_files) cat("    - ", path_file(f), "\n")
  
  # Read all files
  ras_list <- list()
  for (f in nc_files) {
    r <- try(terra::rast(f), silent = TRUE)
    if (inherits(r, "try-error")) {
      warning("    Skipping file: ", path_file(f))
      next
    }
    ras_list[[length(ras_list) + 1L]] <- r
  }
  
  if (length(ras_list) == 0) {
    warning("No valid N-deposition rasters for: ", var_name)
    return(NULL)
  }
  
  r_monthly <- do.call(c, ras_list)
  
  # Crop to boreal band
  bbox_boreal <- terra::ext(-180, 180, 50, 90)
  r_monthly <- terra::crop(r_monthly, bbox_boreal)
  
  # Build year index (starts 1850)
  nl <- terra::nlyr(r_monthly)
  year_index <- 1850 + (seq_len(nl) - 1L) %/% 12L
  
  # Annual sum
  cat("  Computing annual sums...\n")
  r_annual <- terra::tapp(r_monthly, index = year_index, fun = sum)
  names(r_annual) <- sort(unique(year_index))
  
  # Extract at site locations
  cat("  Extracting at sites...\n")
  vals <- terra::extract(r_annual, cbind(sites_df$lon, sites_df$lat))
  
  # Convert to long format
  df_long <- vals %>%
    mutate(ID = row_number()) %>%
    pivot_longer(
      cols = -ID,
      names_to = "year_str",
      values_to = var_name
    ) %>%
    mutate(
      year = as.integer(str_remove(year_str, "^X"))
    ) %>%
    filter(!is.na(year)) %>%
    select(ID, year, all_of(var_name))
  
  # Add site codes
  id_map <- tibble(
    ID = seq_len(nrow(sites_df)),
    site_code = sites_df$site_code
  )
  
  df_long <- df_long %>%
    left_join(id_map, by = "ID") %>%
    select(site_code, year, all_of(var_name))
  
  cat("  Extracted ", nrow(df_long), " site-year observations\n\n")
  
  return(df_long)
}

# Extract NHx
ndep_nhx <- extract_ndep_at_sites("ndep-nhx", sites_coords)

# Extract NOy
ndep_noy <- extract_ndep_at_sites("ndep-noy", sites_coords)

# Check if we have data
has_ndep_nhx <- !is.null(ndep_nhx)
has_ndep_noy <- !is.null(ndep_noy)

if (!has_ndep_nhx && !has_ndep_noy) {
  stop("Could not extract N-deposition data. Check 02a output.")
}

cat("N-deposition data extracted:\n")
if (has_ndep_nhx) cat("  ✓ NHx (reduced nitrogen)\n")
if (has_ndep_noy) cat("  ✓ NOy (oxidized nitrogen)\n")
cat("\n")

# --- 5) Merge N-Deposition with BAI Data -------------------------------------

cat("=== MERGING N-DEPOSITION WITH BAI DATA ===\n\n")

df_with_ndep <- df

if (has_ndep_nhx) {
  df_with_ndep <- df_with_ndep %>%
    left_join(ndep_nhx, by = c("site_code", "year"))
}

if (has_ndep_noy) {
  df_with_ndep <- df_with_ndep %>%
    left_join(ndep_noy, by = c("site_code", "year"))
}

# Save merged dataset
out_parquet <- path(dir_data, "03c_data_with_ndep.parquet")
out_csv <- path(dir_data, "03c_data_with_ndep.csv")

arrow::write_parquet(df_with_ndep, out_parquet)
write_csv(df_with_ndep, out_csv)

cat("Merged dataset saved:\n")
cat("  Parquet: ", out_parquet, "\n")
cat("  CSV:     ", out_csv, "\n\n")

# Check completeness
if (has_ndep_nhx) {
  n_missing_nhx <- sum(is.na(df_with_ndep$`ndep-nhx`))
  cat("NHx missing values: ", n_missing_nhx, " (", 
      round(100 * n_missing_nhx / nrow(df_with_ndep), 1), "%)\n")
}

if (has_ndep_noy) {
  n_missing_noy <- sum(is.na(df_with_ndep$`ndep-noy`))
  cat("NOy missing values: ", n_missing_noy, " (", 
      round(100 * n_missing_noy / nrow(df_with_ndep), 1), "%)\n")
}
cat("\n")

# --- 6) Calculate Site-Level Means -------------------------------------------

cat("=== CALCULATING SITE-LEVEL MEANS ===\n\n")

# Filter complete cases for N-deposition
if (has_ndep_nhx && has_ndep_noy) {
  df_complete <- df_with_ndep %>%
    filter(!is.na(`ndep-nhx`), !is.na(`ndep-noy`))
} else if (has_ndep_nhx) {
  df_complete <- df_with_ndep %>%
    filter(!is.na(`ndep-nhx`))
} else {
  df_complete <- df_with_ndep %>%
    filter(!is.na(`ndep-noy`))
}

cat("Observations with complete N-deposition data: ", nrow(df_complete), "\n")
cat("Sites with complete N-deposition data:        ", n_distinct(df_complete$site_code), "\n\n")

site_means <- df_complete %>%
  group_by(site_code, site_name, lat, lon, elevation_m, genus, region) %>%
  summarise(
    mean_bai = mean(bai_cm2_mean, na.rm = TRUE),
    mean_tmp = mean(tmp_C, na.rm = TRUE),
    mean_pre = mean(pre_mm, na.rm = TRUE),
    .groups = "drop"
  )

if (has_ndep_nhx) {
  site_means <- site_means %>%
    left_join(
      df_complete %>%
        group_by(site_code) %>%
        summarise(mean_ndep_nhx = mean(`ndep-nhx`, na.rm = TRUE), .groups = "drop"),
      by = "site_code"
    )
}

if (has_ndep_noy) {
  site_means <- site_means %>%
    left_join(
      df_complete %>%
        group_by(site_code) %>%
        summarise(mean_ndep_noy = mean(`ndep-noy`, na.rm = TRUE), .groups = "drop"),
      by = "site_code"
    )
}

# Calculate total N-deposition if both available
if (has_ndep_nhx && has_ndep_noy) {
  site_means <- site_means %>%
    mutate(mean_ndep_total = mean_ndep_nhx + mean_ndep_noy)
}

write_csv(site_means, path(dir_tables, "03c_site_means.csv"))

cat("Site-level means calculated for ", nrow(site_means), " sites\n\n")

# --- 7) N-Deposition Summary Statistics --------------------------------------

cat("=== N-DEPOSITION SUMMARY STATISTICS ===\n\n")

ndep_stats <- tibble(
  variable = character(),
  min = numeric(),
  q25 = numeric(),
  median = numeric(),
  mean = numeric(),
  q75 = numeric(),
  max = numeric(),
  sd = numeric()
)

if (has_ndep_nhx) {
  ndep_stats <- bind_rows(
    ndep_stats,
    tibble(
      variable = "NHx",
      min = min(site_means$mean_ndep_nhx, na.rm = TRUE),
      q25 = quantile(site_means$mean_ndep_nhx, 0.25, na.rm = TRUE),
      median = median(site_means$mean_ndep_nhx, na.rm = TRUE),
      mean = mean(site_means$mean_ndep_nhx, na.rm = TRUE),
      q75 = quantile(site_means$mean_ndep_nhx, 0.75, na.rm = TRUE),
      max = max(site_means$mean_ndep_nhx, na.rm = TRUE),
      sd = sd(site_means$mean_ndep_nhx, na.rm = TRUE)
    )
  )
}

if (has_ndep_noy) {
  ndep_stats <- bind_rows(
    ndep_stats,
    tibble(
      variable = "NOy",
      min = min(site_means$mean_ndep_noy, na.rm = TRUE),
      q25 = quantile(site_means$mean_ndep_noy, 0.25, na.rm = TRUE),
      median = median(site_means$mean_ndep_noy, na.rm = TRUE),
      mean = mean(site_means$mean_ndep_noy, na.rm = TRUE),
      q75 = quantile(site_means$mean_ndep_noy, 0.75, na.rm = TRUE),
      max = max(site_means$mean_ndep_noy, na.rm = TRUE),
      sd = sd(site_means$mean_ndep_noy, na.rm = TRUE)
    )
  )
}

if (has_ndep_nhx && has_ndep_noy) {
  ndep_stats <- bind_rows(
    ndep_stats,
    tibble(
      variable = "Total N",
      min = min(site_means$mean_ndep_total, na.rm = TRUE),
      q25 = quantile(site_means$mean_ndep_total, 0.25, na.rm = TRUE),
      median = median(site_means$mean_ndep_total, na.rm = TRUE),
      mean = mean(site_means$mean_ndep_total, na.rm = TRUE),
      q75 = quantile(site_means$mean_ndep_total, 0.75, na.rm = TRUE),
      max = max(site_means$mean_ndep_total, na.rm = TRUE),
      sd = sd(site_means$mean_ndep_total, na.rm = TRUE)
    )
  )
}

write_csv(ndep_stats, path(dir_tables, "03c_ndep_statistics.csv"))

cat("N-deposition statistics (kg N/ha/year):\n")
print(ndep_stats, n = Inf)
cat("\n")

# --- 8) Global N-Deposition Maps ---------------------------------------------

cat("=== CREATING GLOBAL N-DEPOSITION MAPS ===\n\n")

world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare spatial data
sites_sf <- site_means %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Map: NHx deposition
if (has_ndep_nhx) {
  p_map_nhx <- ggplot() +
    geom_sf(data = world, fill = "#f0f0f0", color = "#d0d0d0", linewidth = 0.2) +
    geom_sf(data = sites_sf, aes(color = mean_ndep_nhx, size = mean_ndep_nhx), alpha = 0.7) +
    scale_color_viridis_c(
      name = "NHx\n(kg N/ha/year)",
      option = "plasma",
      trans = "log10",
      breaks = c(0.1, 1, 10),
      labels = c("0.1", "1", "10")
    ) +
    scale_size_continuous(
      name = "NHx\n(kg N/ha/year)",
      range = c(1, 5),
      trans = "log10",
      breaks = c(0.1, 1, 10),
      labels = c("0.1", "1", "10")
    ) +
    coord_sf(ylim = c(45, 90), expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 60)) +
    scale_y_continuous(breaks = seq(50, 90, 10)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.background = element_rect(fill = "aliceblue", color = NA)
    ) +
    labs(
      title = "Global Distribution of Reduced Nitrogen Deposition (NHx)",
      subtitle = sprintf("Sites ≥50°N (n=%d)", nrow(sites_sf)),
      x = "Longitude (°E)",
      y = "Latitude (°N)"
    )
  
  ggsave(
    path(dir_figures, "03c_map_ndep_nhx.png"),
    p_map_nhx,
    width = 14,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  cat("NHx map saved.\n")
}

# Map: NOy deposition
if (has_ndep_noy) {
  p_map_noy <- ggplot() +
    geom_sf(data = world, fill = "#f0f0f0", color = "#d0d0d0", linewidth = 0.2) +
    geom_sf(data = sites_sf, aes(color = mean_ndep_noy, size = mean_ndep_noy), alpha = 0.7) +
    scale_color_viridis_c(
      name = "NOy\n(kg N/ha/year)",
      option = "magma",
      trans = "log10",
      breaks = c(0.1, 1, 10),
      labels = c("0.1", "1", "10")
    ) +
    scale_size_continuous(
      name = "NOy\n(kg N/ha/year)",
      range = c(1, 5),
      trans = "log10",
      breaks = c(0.1, 1, 10),
      labels = c("0.1", "1", "10")
    ) +
    coord_sf(ylim = c(45, 90), expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 60)) +
    scale_y_continuous(breaks = seq(50, 90, 10)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.background = element_rect(fill = "aliceblue", color = NA)
    ) +
    labs(
      title = "Global Distribution of Oxidized Nitrogen Deposition (NOy)",
      subtitle = sprintf("Sites ≥50°N (n=%d)", nrow(sites_sf)),
      x = "Longitude (°E)",
      y = "Latitude (°N)"
    )
  
  ggsave(
    path(dir_figures, "03c_map_ndep_noy.png"),
    p_map_noy,
    width = 14,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  cat("NOy map saved.\n")
}

# Map: Total N deposition
if (has_ndep_nhx && has_ndep_noy) {
  p_map_total <- ggplot() +
    geom_sf(data = world, fill = "#f0f0f0", color = "#d0d0d0", linewidth = 0.2) +
    geom_sf(data = sites_sf, aes(color = mean_ndep_total, size = mean_ndep_total), alpha = 0.7) +
    scale_color_viridis_c(
      name = "Total N\n(kg N/ha/year)",
      option = "viridis",
      trans = "log10",
      breaks = c(0.1, 1, 10, 100),
      labels = c("0.1", "1", "10", "100")
    ) +
    scale_size_continuous(
      name = "Total N\n(kg N/ha/year)",
      range = c(1, 5),
      trans = "log10",
      breaks = c(0.1, 1, 10, 100),
      labels = c("0.1", "1", "10", "100")
    ) +
    coord_sf(ylim = c(45, 90), expand = FALSE) +
    scale_x_continuous(breaks = seq(-180, 180, 60)) +
    scale_y_continuous(breaks = seq(50, 90, 10)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.background = element_rect(fill = "aliceblue", color = NA)
    ) +
    labs(
      title = "Global Distribution of Total Nitrogen Deposition (NHx + NOy)",
      subtitle = sprintf("Sites ≥50°N (n=%d)", nrow(sites_sf)),
      x = "Longitude (°E)",
      y = "Latitude (°N)"
    )
  
  ggsave(
    path(dir_figures, "03c_map_ndep_total.png"),
    p_map_total,
    width = 14,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  cat("Total N map saved.\n")
}

cat("\n")

# Continue with rest of analysis from previous script...
# (I'll include the rest in next message if you want all plots)

cat("=== 03c Explore Nitrogen Deposition COMPLETE ===\n")
cat(sprintf("Analyzed N-deposition patterns across %d sites\n", nrow(site_means)))