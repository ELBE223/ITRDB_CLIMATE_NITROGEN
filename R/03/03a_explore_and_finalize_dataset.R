###############################################################################
# 03a_explore_and_finalize_dataset.R
#
# Description:
#   Exploratory analysis and finalization of BAI + climate dataset.
#   - Filter for complete cases (no missing data)
#   - Detailed summaries (temporal, regional, genus, species)
#   - Visualizations (maps, time series, distributions)
#   - Species richness analysis
#   - Final QC checks
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
  forcats,
  sf,
  rnaturalearth,
  rnaturalearthdata,
  patchwork,
  scales
)

# --- 2) Time Window & Paths --------------------------------------------------

year_min <- 1901L
year_max <- 2010L

step_02b <- "02b_prepare_BAI_climate_for_LMM"
current_name <- "03a_explore_and_finalize_dataset"

input_data_path <- here::here("output", step_02b, "data", "02b_BAI_climate_annual.parquet")

base_out_dir <- here::here("output", current_name)
dir_data     <- path(base_out_dir, "data")
dir_tables   <- path(base_out_dir, "tables")
dir_figures  <- path(base_out_dir, "figures")

dir_create(c(dir_data, dir_tables, dir_figures))

cat("=== 03a Explore and Finalize Dataset ===\n\n")
cat("Input data: ", input_data_path, "\n")
cat("Output:     ", base_out_dir, "\n")
cat("Time span:  ", year_min, "-", year_max, "\n\n")

# --- 3) Load Data ------------------------------------------------------------

if (!file_exists(input_data_path)) {
  stop("Input data not found. Run 02b first.")
}

cat("Loading BAI + climate data...\n")
df_raw <- arrow::read_parquet(input_data_path)

cat("  Total observations: ", nrow(df_raw), "\n")
cat("  Unique sites:       ", n_distinct(df_raw$site_code), "\n")
cat("  Year range:         ", min(df_raw$year), "-", max(df_raw$year), "\n\n")

# --- 4) Data Completeness Analysis -------------------------------------------

cat("=== DATA COMPLETENESS ANALYSIS ===\n\n")

required_cols <- c("bai_cm2_mean", "tmp_C", "pre_mm")
optional_cols <- c("ndep-nhx_ann", "ndep-noy_ann")

# Check all columns
all_check_cols <- c(required_cols, optional_cols[optional_cols %in% names(df_raw)])

completeness_summary <- tibble(
  variable = character(),
  n_total = integer(),
  n_missing = integer(),
  pct_missing = numeric(),
  n_complete = integer(),
  pct_complete = numeric()
)

for (col in all_check_cols) {
  n_total <- nrow(df_raw)
  n_missing <- sum(is.na(df_raw[[col]]))
  n_complete <- n_total - n_missing
  
  completeness_summary <- bind_rows(
    completeness_summary,
    tibble(
      variable = col,
      n_total = n_total,
      n_missing = n_missing,
      pct_missing = round(100 * n_missing / n_total, 2),
      n_complete = n_complete,
      pct_complete = round(100 * n_complete / n_total, 2)
    )
  )
}

write_csv(completeness_summary, path(dir_tables, "03a_data_completeness.csv"))

cat("Data completeness by variable:\n")
print(completeness_summary, n = Inf)
cat("\n")

# --- 5) Filter for Complete Cases --------------------------------------------

cat("=== FILTERING FOR COMPLETE CASES ===\n\n")

has_ndep_nhx <- "ndep-nhx_ann" %in% names(df_raw)
has_ndep_noy <- "ndep-noy_ann" %in% names(df_raw)

filter_cols <- required_cols
if (has_ndep_nhx) filter_cols <- c(filter_cols, "ndep-nhx_ann")
if (has_ndep_noy) filter_cols <- c(filter_cols, "ndep-noy_ann")

df_complete <- df_raw %>%
  filter(
    year >= year_min, 
    year <= year_max,
    if_all(all_of(filter_cols), ~ !is.na(.))
  )

cat("Dataset filtering results:\n")
cat("  Before:       ", nrow(df_raw), " observations from ", 
    n_distinct(df_raw$site_code), " sites\n")
cat("  After:        ", nrow(df_complete), " observations from ", 
    n_distinct(df_complete$site_code), " sites\n")
cat("  Lost:         ", nrow(df_raw) - nrow(df_complete), 
    sprintf(" observations (%.1f%%)\n", 100 * (nrow(df_raw) - nrow(df_complete)) / nrow(df_raw)))
cat("  Sites lost:   ", n_distinct(df_raw$site_code) - n_distinct(df_complete$site_code), "\n\n")

# --- 6) Temporal Coverage Analysis -------------------------------------------

cat("=== TEMPORAL COVERAGE ANALYSIS ===\n\n")

temporal_summary <- df_complete %>%
  group_by(year) %>%
  summarise(
    n_sites = n_distinct(site_code),
    n_observations = n(),
    mean_bai = mean(bai_cm2_mean, na.rm = TRUE),
    sd_bai = sd(bai_cm2_mean, na.rm = TRUE),
    mean_tmp = mean(tmp_C, na.rm = TRUE),
    mean_pre = mean(pre_mm, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(temporal_summary, path(dir_tables, "03a_temporal_coverage.csv"))

cat("Temporal coverage:\n")
cat("  Year range:           ", min(temporal_summary$year), "-", max(temporal_summary$year), "\n")
cat("  Sites per year (min): ", min(temporal_summary$n_sites), "\n")
cat("  Sites per year (max): ", max(temporal_summary$n_sites), "\n")
cat("  Sites per year (mean):", round(mean(temporal_summary$n_sites), 1), "\n\n")

# Temporal coverage plot
p_temporal <- ggplot(temporal_summary, aes(x = year)) +
  geom_line(aes(y = n_sites), color = "#2c3e50", linewidth = 1) +
  geom_point(aes(y = n_sites), color = "#3498db", size = 2, alpha = 0.6) +
  scale_x_continuous(breaks = seq(year_min, year_max, 10)) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Temporal Coverage: Number of Sites per Year",
    subtitle = sprintf("%d-%d (complete cases only)", year_min, year_max),
    x = "Year",
    y = "Number of Sites"
  )

# BAI temporal trend
p_bai_temporal <- ggplot(temporal_summary, aes(x = year, y = mean_bai)) +
  geom_line(color = "#27ae60", linewidth = 1) +
  geom_ribbon(aes(ymin = mean_bai - sd_bai, ymax = mean_bai + sd_bai),
              alpha = 0.2, fill = "#27ae60") +
  scale_x_continuous(breaks = seq(year_min, year_max, 10)) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Mean BAI Over Time",
    subtitle = "Mean ± SD across all sites",
    x = "Year",
    y = "BAI (cm²/year)"
  )

ggsave(
  path(dir_figures, "03a_bai_temporal_trend.png"),
  p_bai_temporal,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

cat("Temporal plots saved.\n\n")

# --- 7) Site-Level Summary ---------------------------------------------------

cat("=== SITE-LEVEL SUMMARY ===\n\n")

site_summary <- df_complete %>%
  group_by(site_code, site_name, region, lat, lon, elevation_m, 
           genus, scientific_name, species_code) %>%
  summarise(
    n_years = n(),
    year_min_data = min(year),
    year_max_data = max(year),
    mean_bai_cm2 = mean(bai_cm2_mean, na.rm = TRUE),
    sd_bai_cm2 = sd(bai_cm2_mean, na.rm = TRUE),
    cv_bai = sd_bai_cm2 / mean_bai_cm2,
    mean_tmp_C = mean(tmp_C, na.rm = TRUE),
    mean_pre_mm = mean(pre_mm, na.rm = TRUE),
    .groups = "drop"
  )

if (has_ndep_nhx) {
  site_summary <- site_summary %>%
    left_join(
      df_complete %>%
        group_by(site_code) %>%
        summarise(mean_ndep_nhx = mean(`ndep-nhx_ann`, na.rm = TRUE), .groups = "drop"),
      by = "site_code"
    )
}

if (has_ndep_noy) {
  site_summary <- site_summary %>%
    left_join(
      df_complete %>%
        group_by(site_code) %>%
        summarise(mean_ndep_noy = mean(`ndep-noy_ann`, na.rm = TRUE), .groups = "drop"),
      by = "site_code"
    )
}

write_csv(site_summary, path(dir_tables, "03a_site_summary.csv"))

cat("Site-level statistics:\n")
cat("  Total sites:       ", nrow(site_summary), "\n")
cat("  Years per site:\n")
cat("    Min:             ", min(site_summary$n_years), "\n")
cat("    Max:             ", max(site_summary$n_years), "\n")
cat("    Mean:            ", round(mean(site_summary$n_years), 1), "\n")
cat("    Median:          ", median(site_summary$n_years), "\n\n")

# --- 8) Regional Summary -----------------------------------------------------

cat("=== REGIONAL SUMMARY ===\n\n")

regional_summary <- df_complete %>%
  group_by(region) %>%
  summarise(
    n_sites = n_distinct(site_code),
    n_observations = n(),
    n_genera = n_distinct(genus),
    mean_bai_cm2 = mean(bai_cm2_mean, na.rm = TRUE),
    sd_bai_cm2 = sd(bai_cm2_mean, na.rm = TRUE),
    mean_tmp_C = mean(tmp_C, na.rm = TRUE),
    mean_pre_mm = mean(pre_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sites))

write_csv(regional_summary, path(dir_tables, "03a_regional_summary.csv"))

cat("Regional distribution:\n")
print(regional_summary, n = Inf)
cat("\n")

# Regional distribution plot
p_regional <- ggplot(regional_summary, aes(x = fct_reorder(region, n_sites), y = n_sites)) +
  geom_col(fill = "#3498db", alpha = 0.8) +
  geom_text(aes(label = comma(n_sites)), hjust = -0.2, size = 4) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11)
  ) +
  labs(
    title = "Regional Distribution of Sites",
    x = NULL,
    y = "Number of Sites"
  )

ggsave(
  path(dir_figures, "03a_regional_distribution.png"),
  p_regional,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

# --- 9) Genus Summary --------------------------------------------------------

cat("=== GENUS SUMMARY ===\n\n")

genus_summary <- df_complete %>%
  group_by(genus) %>%
  summarise(
    n_sites = n_distinct(site_code),
    n_observations = n(),
    mean_bai_cm2 = mean(bai_cm2_mean, na.rm = TRUE),
    sd_bai_cm2 = sd(bai_cm2_mean, na.rm = TRUE),
    mean_tmp_C = mean(tmp_C, na.rm = TRUE),
    mean_pre_mm = mean(pre_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sites))

write_csv(genus_summary, path(dir_tables, "03a_genus_summary.csv"))

cat("Top 10 genera:\n")
print(genus_summary %>% head(10), n = 10)
cat("\n")

# Genus distribution plot
p_genus <- genus_summary %>%
  slice_max(n_sites, n = 15) %>%
  ggplot(aes(x = fct_reorder(genus, n_sites), y = n_sites)) +
  geom_col(fill = "#27ae60", alpha = 0.8) +
  geom_text(aes(label = comma(n_sites)), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11)
  ) +
  labs(
    title = "Top 15 Genera by Number of Sites",
    x = NULL,
    y = "Number of Sites"
  )

ggsave(
  path(dir_figures, "03a_genus_distribution.png"),
  p_genus,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

# --- 10) Species Richness Analysis -------------------------------------------

cat("=== SPECIES RICHNESS ANALYSIS ===\n\n")

# Species richness by region
species_richness_region <- df_complete %>%
  group_by(region) %>%
  summarise(
    n_genera = n_distinct(genus),
    n_species = n_distinct(scientific_name),
    n_sites = n_distinct(site_code),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species))

write_csv(species_richness_region, path(dir_tables, "03a_species_richness_by_region.csv"))

cat("Species richness by region:\n")
print(species_richness_region, n = Inf)
cat("\n")

# Species list with site counts
species_list <- df_complete %>%
  group_by(genus, scientific_name, species_code) %>%
  summarise(
    n_sites = n_distinct(site_code),
    n_observations = n(),
    mean_bai = mean(bai_cm2_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sites))

write_csv(species_list, path(dir_tables, "03a_species_list.csv"))

cat("Total diversity:\n")
cat("  Genera:     ", n_distinct(df_complete$genus), "\n")
cat("  Species:    ", n_distinct(df_complete$scientific_name), "\n\n")

# --- 11) Spatial Distribution Map --------------------------------------------

cat("=== CREATING SPATIAL MAP ===\n\n")

world <- ne_countries(scale = "medium", returnclass = "sf")

sites_sf <- site_summary %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Map with sites colored by genus
top_genera <- genus_summary %>%
  slice_max(n_sites, n = 8) %>%
  pull(genus)

sites_sf <- sites_sf %>%
  mutate(
    genus_display = if_else(genus %in% top_genera, genus, "Other"),
    genus_display = factor(genus_display, levels = c(sort(top_genera), "Other"))
  )

genus_colors <- c(
  "Pinus"      = "#2ca02c",
  "Picea"      = "#1f77b4",
  "Larix"      = "#ff7f0e",
  "Abies"      = "#9467bd",
  "Quercus"    = "#8c564b",
  "Tsuga"      = "#e377c2",
  "Betula"     = "#bcbd22",
  "Juniperus"  = "#17becf",
  "Other"      = "#7f7f7f"
)

p_map <- ggplot() +
  geom_sf(data = world, fill = "#f0f0f0", color = "#d0d0d0", linewidth = 0.2) +
  geom_sf(data = sites_sf, aes(color = genus_display, size = mean_bai_cm2), alpha = 0.7) +
  scale_color_manual(values = genus_colors, name = "Genus") +
  scale_size_continuous(name = "Mean BAI (cm²)", range = c(1, 4), labels = comma) +
  coord_sf(ylim = c(45, 90), expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 60)) +
  scale_y_continuous(breaks = seq(50, 90, 10)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "aliceblue", color = NA)
  ) +
  labs(
    title = "Spatial Distribution of Sites (≥50°N)",
    subtitle = sprintf("Complete dataset: %s sites, %d-%d", 
                       comma(nrow(sites_sf)), year_min, year_max),
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  )

ggsave(
  path(dir_figures, "03a_spatial_map.png"),
  p_map,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("Spatial map saved.\n\n")

# --- 12) Variable Distributions ----------------------------------------------

cat("=== VARIABLE DISTRIBUTIONS ===\n\n")

# BAI distribution
p_bai_dist <- ggplot(df_complete, aes(x = bai_cm2_mean)) +
  geom_histogram(bins = 50, fill = "#27ae60", alpha = 0.7) +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("0.1", "1", "10", "100")
  ) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Distribution of BAI Values",
    x = "BAI (cm²/year, log scale)",
    y = "Count"
  )

# Temperature distribution
p_tmp_dist <- ggplot(df_complete, aes(x = tmp_C)) +
  geom_histogram(bins = 50, fill = "#e74c3c", alpha = 0.7) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Distribution of Temperature",
    x = "Temperature (°C)",
    y = "Count"
  )

# Precipitation distribution
p_pre_dist <- ggplot(df_complete, aes(x = pre_mm)) +
  geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Distribution of Precipitation",
    x = "Precipitation (mm/year)",
    y = "Count"
  )

# Combine plots
p_distributions <- (p_bai_dist | p_tmp_dist | p_pre_dist) +
  plot_annotation(
    title = "Distributions of Key Variables",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(
  path(dir_figures, "03a_variable_distributions.png"),
  p_distributions,
  width = 14,
  height = 5,
  dpi = 300,
  bg = "white"
)

cat("Distribution plots saved.\n\n")

# --- 13) Save Final Dataset --------------------------------------------------

cat("=== SAVING FINAL DATASET ===\n\n")

out_parquet <- path(dir_data, "03a_final_dataset.parquet")
out_csv     <- path(dir_data, "03a_final_dataset.csv")

arrow::write_parquet(df_complete, out_parquet)
write_csv(df_complete, out_csv)

cat("Final dataset saved:\n")
cat("  Parquet: ", out_parquet, "\n")
cat("  CSV:     ", out_csv, "\n\n")

# --- 14) Final Summary Statistics --------------------------------------------

cat("=== FINAL SUMMARY ===\n\n")

final_stats <- tibble(
  metric = c(
    "Total observations",
    "Unique sites",
    "Time span",
    "Years covered",
    "Mean obs per site",
    "Mean BAI (cm²/year)",
    "SD BAI (cm²/year)",
    "Mean temperature (°C)",
    "Mean precipitation (mm/year)",
    "Number of genera",
    "Number of species",
    "Number of regions"
  ),
  value = c(
    format(nrow(df_complete), big.mark = ","),
    format(n_distinct(df_complete$site_code), big.mark = ","),
    paste(year_min, "-", year_max),
    as.character(year_max - year_min + 1),
    format(round(nrow(df_complete) / n_distinct(df_complete$site_code), 1), nsmall = 1),
    format(round(mean(df_complete$bai_cm2_mean, na.rm = TRUE), 3), nsmall = 3),
    format(round(sd(df_complete$bai_cm2_mean, na.rm = TRUE), 3), nsmall = 3),
    format(round(mean(df_complete$tmp_C, na.rm = TRUE), 2), nsmall = 2),
    format(round(mean(df_complete$pre_mm, na.rm = TRUE), 1), nsmall = 1),
    as.character(n_distinct(df_complete$genus)),
    as.character(n_distinct(df_complete$scientific_name)),
    as.character(n_distinct(df_complete$region))
  )
)

print(final_stats, n = Inf)
write_csv(final_stats, path(dir_tables, "03a_final_summary.csv"))

cat("\n=== OUTPUT FILES ===\n")
cat("Data:\n")
cat("  Final dataset (parquet): ", out_parquet, "\n")
cat("  Final dataset (CSV):     ", out_csv, "\n\n")
cat("Tables:\n")
cat("  Data completeness:       ", path(dir_tables, "03a_data_completeness.csv"), "\n")
cat("  Site summary:            ", path(dir_tables, "03a_site_summary.csv"), "\n")
cat("  Temporal coverage:       ", path(dir_tables, "03a_temporal_coverage.csv"), "\n")
cat("  Regional summary:        ", path(dir_tables, "03a_regional_summary.csv"), "\n")
cat("  Genus summary:           ", path(dir_tables, "03a_genus_summary.csv"), "\n")
cat("  Species richness:        ", path(dir_tables, "03a_species_richness_by_region.csv"), "\n")
cat("  Species list:            ", path(dir_tables, "03a_species_list.csv"), "\n")
cat("  Final summary:           ", path(dir_tables, "03a_final_summary.csv"), "\n\n")
cat("Figures:\n")
cat("  ", dir_figures, "\n")

cat("\n=== 03a Explore and Finalize Dataset COMPLETE ===\n")
cat(sprintf("Successfully created final dataset with %s observations from %s sites\n",
            format(nrow(df_complete), big.mark = ","), 
            format(n_distinct(df_complete$site_code), big.mark = ",")))
