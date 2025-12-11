###############################################################################
# 02a_descriptive_plots.R
#
# Description:
#   Exploratory Data Analysis (EDA) and descriptive plots for the merged
#   BAI-Climate-NDeposition dataset from 01d.
#
# Tasks:
#   1. Load merged analysis data.
#   2. Print comprehensive EDA summaries to console.
#   3. Create descriptive plots (saved as PNG).
#
# Input:
#   - output/01d_merge_climate_and_ndep/data/01d_Analysis_Data_Merged.rds
#
# Output:
#   - output/02a_descriptive_plots/figures/*.png
#   - output/02a_descriptive_plots/data/02a_summary_statistics.csv
###############################################################################

# --- 1) Packages -------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  tidyr,
  ggplot2,
  sf,
  rnaturalearth,
  rnaturalearthdata,
  
  viridis,
  patchwork,
  scales,
  corrplot
)

# --- 2) Paths ----------------------------------------------------------------

step_01d <- "01d_merge_climate_and_ndep"
current  <- "02a_descriptive_plots"

in_file <- here::here("output", step_01d, "data", "01d_Analysis_Data_Merged.rds")

base_out_dir <- here::here("output", current)
out_dir_fig  <- fs::path(base_out_dir, "figures")
out_dir_data <- fs::path(base_out_dir, "data")
fs::dir_create(out_dir_fig)
fs::dir_create(out_dir_data)

cat("###############################################################################\n")
cat("# 02a DESCRIPTIVE PLOTS & EXPLORATORY DATA ANALYSIS\n")
cat("###############################################################################\n\n")
cat("Input file:  ", in_file, "\n")
cat("Output dir:  ", base_out_dir, "\n\n")

# --- 3) Load Data ------------------------------------------------------------

if (!file.exists(in_file)) {
  stop("Merged file from 01d not found. Expected at: ", in_file)
}

df <- readRDS(in_file)

cat("=== DATA LOADED ===\n")
cat(sprintf("Dimensions: %d rows x %d columns\n\n", nrow(df), ncol(df)))

# --- 4) Basic Data Overview --------------------------------------------------

cat("###############################################################################\n")
cat("# EXPLORATORY DATA ANALYSIS - CONSOLE OUTPUT\n")
cat("###############################################################################\n\n")

cat("--- COLUMN NAMES ---\n")
print(names(df))

cat("\n--- DATA STRUCTURE ---\n")
str(df, max.level = 1)

cat("\n--- FIRST 10 ROWS ---\n")
print(head(df, 10))

cat("\n--- YEAR RANGE ---\n")
yr_range <- range(df$year, na.rm = TRUE)
cat(sprintf("Years covered: %d to %d (%d years)\n", yr_range[1], yr_range[2], diff(yr_range) + 1))

cat("\n--- UNIQUE COUNTS ---\n")
cat(sprintf("Unique sites:        %d\n", n_distinct(df$site_code)))
cat(sprintf("Unique species:      %d\n", n_distinct(df$species_code)))
cat(sprintf("Unique tree IDs:     %d\n", n_distinct(paste(df$site_code, df$tree_id))))
cat(sprintf("Total observations:  %d\n", nrow(df)))

# --- 5) Summary Statistics ---------------------------------------------------

cat("\n--- SUMMARY STATISTICS FOR KEY VARIABLES ---\n")

key_vars <- c("bai", "temp_annual", "prec_annual", "ndep_total", 
              "ndep_total_3yr", "ndep_total_5yr", "lat", "lon")

summary_stats <- df %>%
  summarise(across(
    all_of(key_vars[key_vars %in% names(df)]),
    list(
      n      = ~sum(!is.na(.)),
      n_na   = ~sum(is.na(.)),
      mean   = ~mean(., na.rm = TRUE),
      sd     = ~sd(., na.rm = TRUE),
      min    = ~min(., na.rm = TRUE),
      q25    = ~quantile(., 0.25, na.rm = TRUE),
      median = ~median(., na.rm = TRUE),
      q75    = ~quantile(., 0.75, na.rm = TRUE),
      max    = ~max(., na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("variable", "statistic"),
    names_pattern = "(.*)_(n|n_na|mean|sd|min|q25|median|q75|max)$"
  ) %>%
  pivot_wider(names_from = statistic, values_from = value)

print(summary_stats, n = 20)

write_csv(summary_stats, path(out_dir_data, "02a_summary_statistics.csv"))
cat("\nSummary statistics saved to: ", path(out_dir_data, "02a_summary_statistics.csv"), "\n")

# --- 6) Missingness Analysis -------------------------------------------------

cat("\n--- MISSINGNESS BY VARIABLE ---\n")

missing_df <- df %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    pct_missing = round(n_missing / nrow(df) * 100, 2)
  ) %>%
  arrange(desc(pct_missing))

print(missing_df, n = 30)

# --- 7) Site-Level Summary ---------------------------------------------------

cat("\n--- SITE-LEVEL SUMMARY ---\n")

site_summary <- df %>%
  group_by(site_code, lat, lon, species_code) %>%
  summarise(
    n_obs        = n(),
    n_trees      = n_distinct(tree_id),
    n_years      = n_distinct(year),
    year_min     = min(year, na.rm = TRUE),
    year_max     = max(year, na.rm = TRUE),
    mean_bai     = mean(bai, na.rm = TRUE),
    sd_bai       = sd(bai, na.rm = TRUE),
    mean_ndep5   = mean(ndep_total_5yr, na.rm = TRUE),
    mean_temp    = mean(temp_annual, na.rm = TRUE),
    mean_prec    = mean(prec_annual, na.rm = TRUE),
    .groups      = "drop"
  )

cat(sprintf("Number of sites: %d\n", nrow(site_summary)))
cat("\nSite summary (first 15 rows):\n")
print(head(site_summary, 15))

write_csv(site_summary, path(out_dir_data, "02a_site_summary.csv"))

# --- 8) Species Distribution -------------------------------------------------

cat("\n--- SPECIES DISTRIBUTION ---\n")

species_counts <- df %>%
  group_by(species_code) %>%
  summarise(
    n_sites = n_distinct(site_code),
    n_trees = n_distinct(paste(site_code, tree_id)),
    n_obs   = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sites))

cat(sprintf("Number of unique species codes: %d\n\n", nrow(species_counts)))
cat("Top 20 species by number of sites:\n")
print(head(species_counts, 20))

# --- 9) Temporal Coverage ----------------------------------------------------

cat("\n--- TEMPORAL COVERAGE ---\n")

yearly_summary <- df %>%
  group_by(year) %>%
  summarise(
    n_obs      = n(),
    n_sites    = n_distinct(site_code),
    n_trees    = n_distinct(paste(site_code, tree_id)),
    mean_bai   = mean(bai, na.rm = TRUE),
    mean_ndep5 = mean(ndep_total_5yr, na.rm = TRUE),
    mean_temp  = mean(temp_annual, na.rm = TRUE),
    .groups    = "drop"
  )

cat("Yearly coverage summary (first and last 10 years):\n")
print(head(yearly_summary, 10))
cat("...\n")
print(tail(yearly_summary, 10))

# --- 10) Correlation Matrix --------------------------------------------------

cat("\n--- CORRELATION MATRIX (KEY VARIABLES) ---\n")

cor_vars <- c("bai", "temp_annual", "prec_annual", "ndep_total_5yr", 
              "temp_summer", "prec_summer", "lat")
cor_vars <- cor_vars[cor_vars %in% names(df)]

cor_data <- df %>%
  select(all_of(cor_vars)) %>%
  filter(complete.cases(.))

cor_matrix <- cor(cor_data)
print(round(cor_matrix, 3))

cat("\n###############################################################################\n")
cat("# GENERATING PLOTS\n")
cat("###############################################################################\n\n")

# --- 11) PLOT 1: World Map of Sites ------------------------------------------

cat("Creating: World map of sites...\n")

site_sf <- st_as_sf(
  site_summary %>% filter(!is.na(lat), !is.na(lon)),
  coords = c("lon", "lat"),
  crs    = 4326
)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p_map <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = site_sf, aes(color = mean_ndep5), size = 2, alpha = 0.8) +
  scale_color_viridis(
    option = "plasma",
    na.value = "black",
    name = expression(paste("Mean N dep. (5-yr, g N ", m^-2, " ", yr^-1, ")"))
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(45, 90), expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(linewidth = 0.1, colour = "grey80"),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  ) +
  labs(
    title = "Spatial Distribution of Study Sites",
    subtitle = "Northern Hemisphere (≥50°N) colored by mean 5-year N deposition",
    x = NULL, y = NULL
  )

ggsave(path(out_dir_fig, "02a_map_sites_ndep.png"), p_map, width = 12, height = 6, dpi = 300)

# --- 12) PLOT 2: BAI Distribution --------------------------------------------

cat("Creating: BAI distribution plots...\n")

p_bai_hist <- ggplot(df, aes(x = bai)) +
  
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "steelblue", alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 1) +
  scale_x_continuous(limits = c(0, quantile(df$bai, 0.99, na.rm = TRUE))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of Basal Area Increment (BAI)",
    subtitle = "Truncated at 99th percentile for visibility",
    x = expression(paste("BAI (", cm^2, " ", yr^-1, ")")),
    y = "Density"
  )

p_bai_log <- ggplot(df, aes(x = log10(bai))) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "steelblue", alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of log10(BAI)",
    x = expression(paste("log10(BAI, ", cm^2, " ", yr^-1, ")")),
    y = "Density"
  )

p_bai_combined <- p_bai_hist + p_bai_log + plot_layout(ncol = 2)
ggsave(path(out_dir_fig, "02a_bai_distribution.png"), p_bai_combined, width = 12, height = 5, dpi = 300)

# --- 13) PLOT 3: N Deposition Distribution -----------------------------------

cat("Creating: N deposition distribution...\n")

p_ndep_hist <- ggplot(df, aes(x = ndep_total_5yr)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "darkgreen", alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of 5-Year Mean N Deposition",
    x = expression(paste("N deposition (g N ", m^-2, " ", yr^-1, ")")),
    y = "Density"
  )

ggsave(path(out_dir_fig, "02a_ndep_distribution.png"), p_ndep_hist, width = 8, height = 5, dpi = 300)

# --- 14) PLOT 4: Climate Distributions ---------------------------------------

cat("Creating: Climate distribution plots...\n")

p_temp <- ggplot(df, aes(x = temp_annual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "coral", alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(title = "Annual Temperature", x = "Temperature (°C)", y = "Density")

p_prec <- ggplot(df, aes(x = prec_annual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "dodgerblue", alpha = 0.7) +
  geom_density(color = "darkblue", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(title = "Annual Precipitation", x = "Precipitation (mm)", y = "Density")

p_climate <- p_temp + p_prec + plot_layout(ncol = 2)
ggsave(path(out_dir_fig, "02a_climate_distributions.png"), p_climate, width = 12, height = 5, dpi = 300)

# --- 15) PLOT 5: Temporal Trends ---------------------------------------------

cat("Creating: Temporal trend plots...\n")

p_time_bai <- ggplot(yearly_summary, aes(x = year, y = mean_bai)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred", fill = "pink", alpha = 0.3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean BAI Over Time",
    x = "Year",
    y = expression(paste("Mean BAI (", cm^2, " ", yr^-1, ")"))
  )

p_time_ndep <- ggplot(yearly_summary, aes(x = year, y = mean_ndep5)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred", fill = "pink", alpha = 0.3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean N Deposition (5-yr) Over Time",
    x = "Year",
    y = expression(paste("N deposition (g N ", m^-2, " ", yr^-1, ")"))
  )

p_time_sites <- ggplot(yearly_summary, aes(x = year, y = n_sites)) +
  geom_line(color = "purple", linewidth = 0.8) +
  geom_area(fill = "purple", alpha = 0.2) +
  theme_minimal(base_size = 12) +
  labs(title = "Number of Sites Per Year", x = "Year", y = "Number of Sites")

p_temporal <- p_time_bai / p_time_ndep / p_time_sites
ggsave(path(out_dir_fig, "02a_temporal_trends.png"), p_temporal, width = 10, height = 10, dpi = 300)

# --- 16) PLOT 6: BAI vs N Deposition Scatter ---------------------------------

cat("Creating: BAI vs N deposition scatter...\n")

set.seed(42)
df_sample <- df %>%
  sample_n(size = min(50000, nrow(df))) %>%
  filter(!is.na(bai), !is.na(ndep_total_5yr))

p_scatter <- ggplot(df_sample, aes(x = ndep_total_5yr, y = bai)) +
  geom_point(alpha = 0.1, size = 0.5, color = "steelblue") +
  geom_smooth(method = "gam", se = TRUE, color = "darkred", fill = "pink") +
  scale_y_continuous(limits = c(0, quantile(df_sample$bai, 0.95, na.rm = TRUE))) +
  theme_minimal(base_size = 12) +
  labs(
    title = "BAI vs. N Deposition (5-year mean)",
    subtitle = "Random sample of 50,000 observations, truncated at 95th percentile BAI",
    x = expression(paste("N deposition (5-yr, g N ", m^-2, " ", yr^-1, ")")),
    y = expression(paste("BAI (", cm^2, " ", yr^-1, ")"))
  )

ggsave(path(out_dir_fig, "02a_scatter_bai_ndep.png"), p_scatter, width = 8, height = 6, dpi = 300)

# --- 17) PLOT 7: BAI by Species (Top 10) -------------------------------------

cat("Creating: BAI by species boxplot...\n")

top_species <- species_counts %>%
  slice_head(n = 10) %>%
  pull(species_code)

df_top_species <- df %>%
  filter(species_code %in% top_species)

p_species <- ggplot(df_top_species, aes(x = reorder(species_code, bai, FUN = median, na.rm = TRUE), y = bai)) +
  geom_boxplot(fill = "lightgreen", outlier.alpha = 0.1, outlier.size = 0.5) +
  scale_y_continuous(limits = c(0, quantile(df_top_species$bai, 0.95, na.rm = TRUE))) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "BAI Distribution by Species (Top 10 by Site Count)",
    subtitle = "Truncated at 95th percentile for visibility",
    x = "Species Code",
    y = expression(paste("BAI (", cm^2, " ", yr^-1, ")"))
  )

ggsave(path(out_dir_fig, "02a_bai_by_species.png"), p_species, width = 10, height = 6, dpi = 300)

# --- 18) PLOT 8: Correlation Heatmap -----------------------------------------

cat("Creating: Correlation heatmap...\n")

png(path(out_dir_fig, "02a_correlation_heatmap.png"), width = 800, height = 800, res = 150)
corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  title = "Correlation Matrix of Key Variables",
  mar = c(0, 0, 2, 0)
)
dev.off()

# --- 19) PLOT 9: Site Mean BAI vs Climate ------------------------------------

cat("Creating: Site-level BAI vs climate plots...\n")

p_site_temp <- ggplot(site_summary, aes(x = mean_temp, y = mean_bai)) +
  geom_point(aes(color = mean_ndep5), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  scale_color_viridis(option = "plasma", name = "N dep.") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Site Mean BAI vs. Temperature",
    x = "Mean Annual Temperature (°C)",
    y = expression(paste("Mean BAI (", cm^2, " ", yr^-1, ")"))
  )

p_site_prec <- ggplot(site_summary, aes(x = mean_prec, y = mean_bai)) +
  geom_point(aes(color = mean_ndep5), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  scale_color_viridis(option = "plasma", name = "N dep.") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Site Mean BAI vs. Precipitation",
    x = "Mean Annual Precipitation (mm)",
    y = expression(paste("Mean BAI (", cm^2, " ", yr^-1, ")"))
  )

p_site_climate <- p_site_temp + p_site_prec + plot_layout(ncol = 2, guides = "collect")
ggsave(path(out_dir_fig, "02a_site_bai_vs_climate.png"), p_site_climate, width = 14, height = 6, dpi = 300)

# --- 20) PLOT 10: Decadal Patterns -------------------------------------------

cat("Creating: Decadal pattern plots...\n")

df_decade <- df %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade) %>%
  summarise(
    n_obs      = n(),
    n_sites    = n_distinct(site_code),
    mean_bai   = mean(bai, na.rm = TRUE),
    sd_bai     = sd(bai, na.rm = TRUE),
    mean_ndep5 = mean(ndep_total_5yr, na.rm = TRUE),
    mean_temp  = mean(temp_annual, na.rm = TRUE),
    .groups    = "drop"
  )

p_decade_bai <- ggplot(df_decade, aes(x = factor(decade), y = mean_bai)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_bai - sd_bai/sqrt(n_obs), 
                    ymax = mean_bai + sd_bai/sqrt(n_obs)), width = 0.3) +
  theme_minimal(base_size = 12) +
  labs(title = "Mean BAI by Decade", x = "Decade", 
       y = expression(paste("Mean BAI (", cm^2, " ", yr^-1, ")")))

p_decade_ndep <- ggplot(df_decade, aes(x = factor(decade), y = mean_ndep5)) +
  geom_col(fill = "darkgreen", alpha = 0.8) +
  theme_minimal(base_size = 12) +
  labs(title = "Mean N Deposition by Decade", x = "Decade",
       y = expression(paste("N dep. (g N ", m^-2, " ", yr^-1, ")")))

p_decade <- p_decade_bai + p_decade_ndep + plot_layout(ncol = 2)
ggsave(path(out_dir_fig, "02a_decadal_patterns.png"), p_decade, width = 14, height = 5, dpi = 300)

# --- 21) Summary Output ------------------------------------------------------

cat("\n###############################################################################\n")
cat("# 02a DESCRIPTIVE PLOTS COMPLETE\n")
cat("###############################################################################\n\n")

cat("PLOTS SAVED:\n")
plot_files <- dir_ls(out_dir_fig, glob = "*.png")
for (f in plot_files) {
  cat("  - ", path_file(f), "\n")
}

cat("\nDATA FILES SAVED:\n")
data_files <- dir_ls(out_dir_data)
for (f in data_files) {
  cat("  - ", path_file(f), "\n")
}

cat("\n=== 02a COMPLETE ===\n")
