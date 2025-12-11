###############################################################################
# 01e_QC_01d_explore.R
#
# Quick QC and exploration of:
#   output/01d_merge_climate_and_ndep/data/01d_Analysis_Data_Merged.rds
#
# Tasks:
#   1. Basic checks (dimensions, year range, missingness).
#   2. Site-level summaries (mean BAI, mean N deposition, climate).
#   3. World map of BAI sites colored by mean N deposition.
#   4. Simple distribution plots (N deposition, BAI).
###############################################################################

# --- 1) Packages -------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  ggplot2,
  sf,
  rnaturalearth,
  rnaturalearthdata,
  viridis,     # for colours
  patchwork    # for combining plots (optional)
)

# --- 2) Paths ----------------------------------------------------------------

step_01d <- "01d_merge_climate_and_ndep"
in_file  <- here::here("output", step_01d, "data", "01d_Analysis_Data_Merged.rds")

out_dir_fig <- here::here("output", step_01d, "figures")
fs::dir_create(out_dir_fig)

cat("=== 01e QC & Exploration for 01d ===\n\n")
cat("Input file: ", in_file, "\n")

# --- 3) Load data ------------------------------------------------------------

if (!file.exists(in_file)) {
  stop("Merged file from 01d not found. Expected at: ", in_file)
}

df <- readRDS(in_file)

cat(sprintf("Data loaded. %d rows, %d cols.\n", nrow(df), ncol(df)))

# --- 4) Basic QC summaries ---------------------------------------------------

cat("\n--- BASIC STRUCTURE ---\n")
print(str(df, max.level = 1))

cat("\n--- YEAR RANGE ---\n")
yr_range <- range(df$year, na.rm = TRUE)
cat(sprintf("Years: %d – %d\n", yr_range[1], yr_range[2]))

cat("\n--- KEY VARIABLE MISSINGNESS ---\n")

qc_vars <- c("bai", "temp_annual", "prec_annual",
             "ndep_total", "ndep_total_3yr", "ndep_total_5yr")

missing_summary <- df %>%
  summarise(across(all_of(qc_vars),
                   ~ sum(is.na(.)),
                   .names = "na_{.col}")) %>%
  tidyr::pivot_longer(everything(),
                      names_to = "var",
                      values_to = "n_na") %>%
  mutate(
    var    = sub("^na_", "", var),
    pct_na = n_na / nrow(df) * 100
  )

print(missing_summary)

# --- 5) Site-level summary ---------------------------------------------------

cat("\n--- SITE-LEVEL SUMMARY ---\n")

site_summary <- df %>%
  group_by(site_code, lat, lon) %>%
  summarise(
    n_obs      = n(),
    n_years    = n_distinct(year),
    year_min   = min(year, na.rm = TRUE),
    year_max   = max(year, na.rm = TRUE),
    mean_bai   = mean(bai, na.rm = TRUE),
    med_bai    = median(bai, na.rm = TRUE),
    mean_ndep  = mean(ndep_total, na.rm = TRUE),
    mean_ndep5 = mean(ndep_total_5yr, na.rm = TRUE),
    mean_temp  = mean(temp_annual, na.rm = TRUE),
    mean_prec  = mean(prec_annual, na.rm = TRUE),
    .groups    = "drop"
  )

cat(sprintf("Unique sites: %d\n", nrow(site_summary)))
cat("Example of site summary:\n")
print(head(site_summary, 10))

# Save site summary as CSV (optional)
site_summary_file <- file.path(out_dir_fig, "01e_site_summary.csv")
readr::write_csv(site_summary, site_summary_file)
cat("Site summary saved to: ", site_summary_file, "\n")

# --- 6) World map of sites (mean N deposition) -------------------------------

cat("\n--- PLOTTING WORLD MAP OF SITES ---\n")

# Turn sites into sf object
site_sf <- st_as_sf(
  site_summary %>% filter(!is.na(lat), !is.na(lon)),
  coords = c("lon", "lat"),
  crs    = 4326
)

# Base world
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p_map <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = site_sf,
          aes(color = mean_ndep5),
          size = 1.8,
          alpha = 0.8) +
  scale_color_viridis(
    option = "C",
    na.value = "black",
    name = expression(paste("Mean N deposition (5-yr, ", g ~ N ~ m^-2 ~ yr^-1, ")"))
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(linewidth = 0.1, colour = "grey80"),
    legend.position  = "bottom"
  ) +
  labs(
    title = "Global distribution of BAI sites",
    subtitle = "Points colored by mean 5-year total N deposition (1901–2010)",
    x = NULL, y = NULL
  )

map_file <- file.path(out_dir_fig, "01e_worldmap_sites_mean_ndep5.png")
ggsave(map_file, p_map, width = 10, height = 5, dpi = 300)
cat("World map saved to: ", map_file, "\n")

# --- 7) Simple distributions: BAI and N deposition ---------------------------

cat("\n--- PLOTTING DISTRIBUTIONS ---\n")

# For distributions we work on sample of rows (to speed up)
set.seed(123)
df_sample <- df %>%
  sample_n(size = min(200000, nrow(df))) %>%
  filter(!is.na(bai), !is.na(ndep_total_5yr))

p_hist_ndep <- ggplot(df_sample, aes(x = ndep_total_5yr)) +
  geom_histogram(bins = 50) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Distribution of 5-year mean N deposition",
    x = expression(paste("N deposition (5-yr, ", g ~ N ~ m^-2 ~ yr^-1, ")")),
    y = "Count"
  )

p_hist_bai <- ggplot(df_sample, aes(x = bai)) +
  geom_histogram(bins = 50) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Distribution of Basal Area Increment (BAI)",
    x = expression(paste("BAI (", mm^2 ~ yr^-1, "?)")),
    y = "Count"
  )

# Combine with patchwork (if installed)
p_dist <- p_hist_ndep + p_hist_bai + plot_layout(ncol = 2)

dist_file <- file.path(out_dir_fig, "01e_distributions_ndep5_bai.png")
ggsave(dist_file, p_dist, width = 10, height = 4, dpi = 300)
cat("Distribution plots saved to: ", dist_file, "\n")

# --- 8) Simple time-coverage check -------------------------------------------

cat("\n--- TIME COVERAGE (PER DECADE) ---\n")

df_decade <- df %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade) %>%
  summarise(
    n_obs       = n(),
    n_sites     = n_distinct(site_code),
    mean_ndep5  = mean(ndep_total_5yr, na.rm = TRUE),
    mean_bai    = mean(bai, na.rm = TRUE),
    .groups     = "drop"
  )

print(df_decade)

p_decade <- ggplot(df_decade, aes(x = factor(decade), y = mean_ndep5, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 11) +
  labs(
    title = "Mean 5-year N deposition over decades",
    x = "Decade",
    y = expression(paste("Mean N deposition (5-yr, ", g ~ N ~ m^-2 ~ yr^-1, ")"))
  )

decade_file <- file.path(out_dir_fig, "01e_decadal_ndep5.png")
ggsave(decade_file, p_decade, width = 7, height = 4, dpi = 300)
cat("Decadal N deposition plot saved to: ", decade_file, "\n")

cat("\n=== 01e QC & Exploration complete ===\n")
