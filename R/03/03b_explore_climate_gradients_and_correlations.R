###############################################################################
# 03b_explore_climate_gradients_and_correlations.R
#
# Description:
#   Exploratory analysis of climate and nitrogen deposition gradients.
#   - Latitudinal and elevational gradients
#   - Climate space analysis
#   - Correlation analysis between variables
#   - Outlier detection and diagnostics
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
  corrplot,
  scales,
  viridis
)

# --- 2) Paths ----------------------------------------------------------------

step_03a <- "03a_explore_and_finalize_dataset"
current_name <- "03b_explore_climate_gradients_and_correlations"

input_data_path <- here::here("output", step_03a, "data", "03a_final_dataset.parquet")

base_out_dir <- here::here("output", current_name)
dir_data     <- path(base_out_dir, "data")
dir_tables   <- path(base_out_dir, "tables")
dir_figures  <- path(base_out_dir, "figures")

dir_create(c(dir_data, dir_tables, dir_figures))

cat("=== 03b Explore Climate Gradients and Correlations ===\n\n")
cat("Input data: ", input_data_path, "\n")
cat("Output:     ", base_out_dir, "\n\n")

# --- 3) Load Data ------------------------------------------------------------

if (!file_exists(input_data_path)) {
  stop("Input data not found. Run 03a first.")
}

cat("Loading final dataset...\n")
df <- arrow::read_parquet(input_data_path)

cat("  Observations: ", nrow(df), "\n")
cat("  Sites:        ", n_distinct(df$site_code), "\n\n")

# Check for N-deposition columns
has_ndep_nhx <- "ndep-nhx_ann" %in% names(df)
has_ndep_noy <- "ndep-noy_ann" %in% names(df)

# --- 4) Site-Level Means for Gradient Analysis ------------------------------

cat("=== CALCULATING SITE-LEVEL MEANS ===\n\n")

site_means <- df %>%
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
      df %>%
        group_by(site_code) %>%
        summarise(mean_ndep_nhx = mean(`ndep-nhx_ann`, na.rm = TRUE), .groups = "drop"),
      by = "site_code"
    )
}

if (has_ndep_noy) {
  site_means <- site_means %>%
    left_join(
      df %>%
        group_by(site_code) %>%
        summarise(mean_ndep_noy = mean(`ndep-noy_ann`, na.rm = TRUE), .groups = "drop"),
      by = "site_code"
    )
}

write_csv(site_means, path(dir_tables, "03b_site_means.csv"))

cat("Site-level means calculated for ", nrow(site_means), " sites\n\n")

# --- 5) Latitudinal Gradients ------------------------------------------------

cat("=== LATITUDINAL GRADIENTS ===\n\n")

# Temperature vs latitude
p_lat_tmp <- ggplot(site_means, aes(x = lat, y = mean_tmp)) +
  geom_point(alpha = 0.5, size = 2, color = "#e74c3c") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
  scale_x_continuous(breaks = seq(50, 90, 5)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Temperature vs Latitude",
    x = "Latitude (°N)",
    y = "Mean Temperature (°C)"
  )

# Precipitation vs latitude
p_lat_pre <- ggplot(site_means, aes(x = lat, y = mean_pre)) +
  geom_point(alpha = 0.5, size = 2, color = "#3498db") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
  scale_x_continuous(breaks = seq(50, 90, 5)) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Precipitation vs Latitude",
    x = "Latitude (°N)",
    y = "Mean Precipitation (mm/year)"
  )

# BAI vs latitude
p_lat_bai <- ggplot(site_means, aes(x = lat, y = mean_bai)) +
  geom_point(alpha = 0.5, size = 2, color = "#27ae60") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
  scale_x_continuous(breaks = seq(50, 90, 5)) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("0.1", "1", "10", "100")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "BAI vs Latitude",
    x = "Latitude (°N)",
    y = "Mean BAI (cm²/year, log scale)"
  )

# N-deposition vs latitude
if (has_ndep_nhx) {
  p_lat_ndep_nhx <- ggplot(site_means, aes(x = lat, y = mean_ndep_nhx)) +
    geom_point(alpha = 0.5, size = 2, color = "#9b59b6") +
    geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
    scale_x_continuous(breaks = seq(50, 90, 5)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "N-deposition (NHx) vs Latitude",
      x = "Latitude (°N)",
      y = "Mean NHx (kg N/ha/year)"
    )
}

if (has_ndep_noy) {
  p_lat_ndep_noy <- ggplot(site_means, aes(x = lat, y = mean_ndep_noy)) +
    geom_point(alpha = 0.5, size = 2, color = "#e67e22") +
    geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
    scale_x_continuous(breaks = seq(50, 90, 5)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "N-deposition (NOy) vs Latitude",
      x = "Latitude (°N)",
      y = "Mean NOy (kg N/ha/year)"
    )
}

# Combine latitudinal plots
if (has_ndep_nhx && has_ndep_noy) {
  p_lat_combined <- (p_lat_tmp | p_lat_pre) / (p_lat_bai | p_lat_ndep_nhx) / p_lat_ndep_noy +
    plot_annotation(
      title = "Latitudinal Gradients",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
} else if (has_ndep_nhx) {
  p_lat_combined <- (p_lat_tmp | p_lat_pre) / (p_lat_bai | p_lat_ndep_nhx) +
    plot_annotation(
      title = "Latitudinal Gradients",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
} else {
  p_lat_combined <- (p_lat_tmp | p_lat_pre) / p_lat_bai +
    plot_annotation(
      title = "Latitudinal Gradients",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}

ggsave(
  path(dir_figures, "03b_latitudinal_gradients.png"),
  p_lat_combined,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("Latitudinal gradient plots saved.\n\n")

# --- 6) Elevational Gradients ------------------------------------------------

cat("=== ELEVATIONAL GRADIENTS ===\n\n")

# Remove sites with missing elevation
site_means_elev <- site_means %>%
  filter(!is.na(elevation_m))

cat("Sites with elevation data: ", nrow(site_means_elev), "\n\n")

if (nrow(site_means_elev) > 10) {
  
  # Temperature vs elevation
  p_elev_tmp <- ggplot(site_means_elev, aes(x = elevation_m, y = mean_tmp)) +
    geom_point(alpha = 0.5, size = 2, color = "#e74c3c") +
    geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
    scale_x_continuous(labels = comma) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "Temperature vs Elevation",
      x = "Elevation (m)",
      y = "Mean Temperature (°C)"
    )
  
  # Precipitation vs elevation
  p_elev_pre <- ggplot(site_means_elev, aes(x = elevation_m, y = mean_pre)) +
    geom_point(alpha = 0.5, size = 2, color = "#3498db") +
    geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "Precipitation vs Elevation",
      x = "Elevation (m)",
      y = "Mean Precipitation (mm/year)"
    )
  
  # BAI vs elevation
  p_elev_bai <- ggplot(site_means_elev, aes(x = elevation_m, y = mean_bai)) +
    geom_point(alpha = 0.5, size = 2, color = "#27ae60") +
    geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", linewidth = 1) +
    scale_x_continuous(labels = comma) +
    scale_y_log10(
      breaks = c(0.1, 1, 10, 100),
      labels = c("0.1", "1", "10", "100")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "BAI vs Elevation",
      x = "Elevation (m)",
      y = "Mean BAI (cm²/year, log scale)"
    )
  
  p_elev_combined <- (p_elev_tmp | p_elev_pre) / p_elev_bai +
    plot_annotation(
      title = "Elevational Gradients",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  ggsave(
    path(dir_figures, "03b_elevational_gradients.png"),
    p_elev_combined,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  cat("Elevational gradient plots saved.\n\n")
}

# --- 7) Climate Space Analysis -----------------------------------------------

cat("=== CLIMATE SPACE ANALYSIS ===\n\n")

# Temperature-Precipitation space colored by BAI
p_climate_space_bai <- ggplot(site_means, aes(x = mean_tmp, y = mean_pre, color = mean_bai)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c(
    trans = "log10",
    name = "Mean BAI\n(cm²/year)",
    breaks = c(0.1, 1, 10, 100),
    labels = c("0.1", "1", "10", "100")
  ) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Climate Space: BAI Distribution",
    x = "Mean Temperature (°C)",
    y = "Mean Precipitation (mm/year)"
  )

ggsave(
  path(dir_figures, "03b_climate_space_bai.png"),
  p_climate_space_bai,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Temperature-Precipitation space colored by region
p_climate_space_region <- ggplot(site_means, aes(x = mean_tmp, y = mean_pre, color = region)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set2", name = "Region") +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Climate Space: Regional Distribution",
    x = "Mean Temperature (°C)",
    y = "Mean Precipitation (mm/year)"
  )

ggsave(
  path(dir_figures, "03b_climate_space_region.png"),
  p_climate_space_region,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# N-deposition vs climate if available
if (has_ndep_nhx && has_ndep_noy) {
  
  site_means <- site_means %>%
    mutate(total_ndep = mean_ndep_nhx + mean_ndep_noy)
  
  p_climate_ndep <- ggplot(site_means, aes(x = mean_tmp, y = mean_pre, color = total_ndep)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_viridis_c(name = "Total N-dep\n(kg N/ha/year)") +
    scale_y_continuous(labels = comma) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "Climate Space: N-deposition Distribution",
      x = "Mean Temperature (°C)",
      y = "Mean Precipitation (mm/year)"
    )
  
  ggsave(
    path(dir_figures, "03b_climate_space_ndep.png"),
    p_climate_ndep,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

cat("Climate space plots saved.\n\n")

# --- 8) Correlation Analysis -------------------------------------------------

cat("=== CORRELATION ANALYSIS ===\n\n")

# Select numeric variables for correlation
cor_vars <- c("mean_bai", "mean_tmp", "mean_pre", "lat", "elevation_m")

if (has_ndep_nhx) cor_vars <- c(cor_vars, "mean_ndep_nhx")
if (has_ndep_noy) cor_vars <- c(cor_vars, "mean_ndep_noy")

# Filter complete cases
site_means_cor <- site_means %>%
  select(all_of(cor_vars)) %>%
  drop_na()

cat("Sites with complete data for correlation: ", nrow(site_means_cor), "\n\n")

# Calculate correlation matrix
cor_matrix <- cor(site_means_cor, use = "complete.obs")

# Save correlation matrix
cor_df <- as.data.frame(cor_matrix)
cor_df <- tibble::rownames_to_column(cor_df, "variable")
write_csv(cor_df, path(dir_tables, "03b_correlation_matrix.csv"))

# Correlation plot
png(
  path(dir_figures, "03b_correlation_matrix.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300,
  bg = "white"
)

corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  number.cex = 0.8,
  cl.cex = 0.8,
  title = "Correlation Matrix: Site-Level Means",
  mar = c(0, 0, 2, 0)
)

dev.off()

cat("Correlation analysis saved.\n\n")

# Pairwise correlations table
cor_pairs <- tibble::tibble(
  var1 = character(),
  var2 = character(),
  correlation = numeric(),
  abs_correlation = numeric()
)

for (i in 1:(length(cor_vars) - 1)) {
  for (j in (i + 1):length(cor_vars)) {
    cor_pairs <- bind_rows(
      cor_pairs,
      tibble(
        var1 = cor_vars[i],
        var2 = cor_vars[j],
        correlation = cor_matrix[i, j],
        abs_correlation = abs(cor_matrix[i, j])
      )
    )
  }
}

cor_pairs <- cor_pairs %>%
  arrange(desc(abs_correlation))

write_csv(cor_pairs, path(dir_tables, "03b_correlation_pairs.csv"))

cat("Top 10 strongest correlations:\n")
print(cor_pairs %>% head(10), n = 10)
cat("\n")

# --- 9) Outlier Detection ----------------------------------------------------

cat("=== OUTLIER DETECTION ===\n\n")

# Function to detect outliers using IQR method
detect_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 3 * iqr
  upper <- q3 + 3 * iqr
  x < lower | x > upper
}

# Detect outliers in key variables
outlier_flags <- site_means %>%
  select(site_code, mean_bai, mean_tmp, mean_pre) %>%
  mutate(
    outlier_bai = detect_outliers(mean_bai),
    outlier_tmp = detect_outliers(mean_tmp),
    outlier_pre = detect_outliers(mean_pre),
    outlier_any = outlier_bai | outlier_tmp | outlier_pre
  )

outliers <- outlier_flags %>%
  filter(outlier_any) %>%
  left_join(site_means, by = "site_code")

write_csv(outliers, path(dir_tables, "03b_outliers.csv"))

cat("Outliers detected:\n")
cat("  BAI outliers:          ", sum(outlier_flags$outlier_bai), "\n")
cat("  Temperature outliers:  ", sum(outlier_flags$outlier_tmp), "\n")
cat("  Precipitation outliers:", sum(outlier_flags$outlier_pre), "\n")
cat("  Total unique sites:    ", sum(outlier_flags$outlier_any), "\n\n")

# Visualize outliers
p_outliers_bai <- ggplot(outlier_flags, aes(x = mean_bai, fill = outlier_bai)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("0.1", "1", "10", "100")
  ) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
    name = "Outlier"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "BAI Distribution with Outliers",
    x = "Mean BAI (cm²/year, log scale)",
    y = "Count"
  )

p_outliers_tmp <- ggplot(outlier_flags, aes(x = mean_tmp, fill = outlier_tmp)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
    name = "Outlier"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Temperature Distribution with Outliers",
    x = "Mean Temperature (°C)",
    y = "Count"
  )

p_outliers_pre <- ggplot(outlier_flags, aes(x = mean_pre, fill = outlier_pre)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
    name = "Outlier"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "Precipitation Distribution with Outliers",
    x = "Mean Precipitation (mm/year)",
    y = "Count"
  )

p_outliers_combined <- p_outliers_bai / p_outliers_tmp / p_outliers_pre +
  plot_annotation(
    title = "Outlier Detection (3×IQR method)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(
  path(dir_figures, "03b_outlier_detection.png"),
  p_outliers_combined,
  width = 10,
  height = 12,
  dpi = 300,
  bg = "white"
)

cat("Outlier detection plots saved.\n\n")

# --- 10) Summary Statistics --------------------------------------------------

cat("=== GRADIENT SUMMARY STATISTICS ===\n\n")

gradient_stats <- tibble(
  variable = character(),
  min = numeric(),
  q25 = numeric(),
  median = numeric(),
  mean = numeric(),
  q75 = numeric(),
  max = numeric(),
  sd = numeric()
)

for (var in cor_vars) {
  if (var %in% names(site_means_cor)) {
    x <- site_means_cor[[var]]
    gradient_stats <- bind_rows(
      gradient_stats,
      tibble(
        variable = var,
        min = min(x, na.rm = TRUE),
        q25 = quantile(x, 0.25, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE),
        q75 = quantile(x, 0.75, na.rm = TRUE),
        max = max(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE)
      )
    )
  }
}

write_csv(gradient_stats, path(dir_tables, "03b_gradient_statistics.csv"))

print(gradient_stats)
cat("\n")

# --- 11) Output Summary ------------------------------------------------------

cat("=== OUTPUT FILES ===\n")
cat("Tables:\n")
cat("  Site means:          ", path(dir_tables, "03b_site_means.csv"), "\n")
cat("  Correlation matrix:  ", path(dir_tables, "03b_correlation_matrix.csv"), "\n")
cat("  Correlation pairs:   ", path(dir_tables, "03b_correlation_pairs.csv"), "\n")
cat("  Outliers:            ", path(dir_tables, "03b_outliers.csv"), "\n")
cat("  Gradient statistics: ", path(dir_tables, "03b_gradient_statistics.csv"), "\n\n")
cat("Figures:\n")
cat("  ", dir_figures, "\n")

cat("\n=== 03b Explore Climate Gradients and Correlations COMPLETE ===\n")