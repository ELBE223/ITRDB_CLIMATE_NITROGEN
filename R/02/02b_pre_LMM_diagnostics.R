###############################################################################
# 02b_pre_LMM_diagnostics.R
#
# Description:
#   Statistical diagnostics before fitting Linear Mixed Models (LMM).
#   Tests assumptions and explores data structure.
#
# Tasks:
#   1. Distributional checks (normality, transformation assessment)
#   2. Variance structure exploration (Levene's test, ICC)
#   3. Bivariate relationships (correlations, simple regressions)
#   4. Collinearity diagnostics (VIF, correlation significance)
#   5. Group-level sample size checks (balance across sites/species)
#
# Input:
#   - output/01d_merge_climate_and_ndep/data/01d_Analysis_Data_Merged.rds
#
# Output:
#   - output/02b_pre_LMM_diagnostics/figures/*.png
#   - output/02b_pre_LMM_diagnostics/data/*.csv
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
  car,
  lme4,
  performance,
  patchwork
)

# --- 2) Paths ----------------------------------------------------------------

step_01d <- "01d_merge_climate_and_ndep"
current  <- "02b_pre_LMM_diagnostics"

in_file <- here::here("output", step_01d, "data", "01d_Analysis_Data_Merged.rds")

base_out_dir <- here::here("output", current)
out_dir_fig  <- fs::path(base_out_dir, "figures")
out_dir_data <- fs::path(base_out_dir, "data")
fs::dir_create(out_dir_fig)
fs::dir_create(out_dir_data)

cat("###############################################################################\n")
cat("# 02b PRE-LMM DIAGNOSTICS\n")
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

cat("--- FIRST 10 ROWS ---\n")
print(head(df, 10))

# Remove rows with missing key variables
df_complete <- df %>%
  filter(!is.na(bai), !is.na(ndep_total_5yr), !is.na(temp_annual), !is.na(prec_annual))

cat(sprintf("\nComplete cases for analysis: %d (%.1f%% of total)\n", 
            nrow(df_complete), 100 * nrow(df_complete) / nrow(df)))

cat("\n###############################################################################\n")
cat("# 1. DISTRIBUTIONAL CHECKS\n")
cat("###############################################################################\n\n")

# --- 1.1 Shapiro-Wilk Test (on sample due to size limits) --------------------

cat("--- SHAPIRO-WILK NORMALITY TEST ---\n")
cat("(Using random sample of 5000 observations due to test limitations)\n\n")

set.seed(123)
sample_bai <- sample(df_complete$bai, min(5000, nrow(df_complete)))
sample_log_bai <- log(sample_bai)

sw_raw <- shapiro.test(sample_bai)
sw_log <- shapiro.test(sample_log_bai)

cat("BAI (raw):\n")
cat(sprintf("  W = %.6f, p-value = %.2e\n", sw_raw$statistic, sw_raw$p.value))
cat(sprintf("  Interpretation: %s\n\n", 
            ifelse(sw_raw$p.value < 0.05, "NOT normally distributed (p < 0.05)", "Approximately normal")))

cat("BAI (log-transformed):\n")
cat(sprintf("  W = %.6f, p-value = %.2e\n", sw_log$statistic, sw_log$p.value))
cat(sprintf("  Interpretation: %s\n\n", 
            ifelse(sw_log$p.value < 0.05, "NOT normally distributed (p < 0.05)", "Approximately normal")))

# --- 1.2 Skewness and Kurtosis -----------------------------------------------

cat("--- SKEWNESS AND KURTOSIS ---\n")

calc_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  (sum((x - m)^3) / n) / (s^3)
}

calc_kurtosis <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  (sum((x - m)^4) / n) / (s^4) - 3
}

skew_raw <- calc_skewness(df_complete$bai)
kurt_raw <- calc_kurtosis(df_complete$bai)
skew_log <- calc_skewness(log(df_complete$bai))
kurt_log <- calc_kurtosis(log(df_complete$bai))

cat(sprintf("BAI (raw):      Skewness = %.3f, Kurtosis = %.3f\n", skew_raw, kurt_raw))
cat(sprintf("BAI (log):      Skewness = %.3f, Kurtosis = %.3f\n", skew_log, kurt_log))
cat("(Ideal: Skewness ≈ 0, Kurtosis ≈ 0 for normal distribution)\n")
cat(sprintf("\nRecommendation: %s\n\n", 
            ifelse(abs(skew_log) < abs(skew_raw), 
                   "Log-transformation IMPROVES normality", 
                   "Log-transformation does NOT improve normality")))

# --- 1.3 Q-Q Plots -----------------------------------------------------------

cat("Creating: Q-Q plots for normality assessment...\n")

p_qq_raw <- ggplot(data.frame(bai = sample_bai), aes(sample = bai)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal(base_size = 12) +
  labs(title = "Q-Q Plot: BAI (raw)", x = "Theoretical Quantiles", y = "Sample Quantiles")

p_qq_log <- ggplot(data.frame(bai = sample_log_bai), aes(sample = bai)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  theme_minimal(base_size = 12) +
  labs(title = "Q-Q Plot: BAI (log-transformed)", x = "Theoretical Quantiles", y = "Sample Quantiles")

p_qq <- p_qq_raw + p_qq_log
ggsave(path(out_dir_fig, "02b_qq_plots.png"), p_qq, width = 12, height = 5, dpi = 300)

# --- 1.4 Distribution Comparison Plot ----------------------------------------

cat("Creating: Distribution comparison plots...\n")

p_dist_raw <- ggplot(df_complete, aes(x = bai)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "steelblue", alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 1) +
  scale_x_continuous(limits = c(0, quantile(df_complete$bai, 0.99))) +
  theme_minimal(base_size = 12) +
  labs(title = "BAI (raw)", x = "BAI", y = "Density")

p_dist_log <- ggplot(df_complete, aes(x = log(bai))) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "steelblue", alpha = 0.7) +
  geom_density(color = "darkred", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(title = "BAI (log-transformed)", x = "log(BAI)", y = "Density")

p_dist <- p_dist_raw + p_dist_log
ggsave(path(out_dir_fig, "02b_distribution_comparison.png"), p_dist, width = 12, height = 5, dpi = 300)

cat("\n###############################################################################\n")
cat("# 2. VARIANCE STRUCTURE EXPLORATION\n")
cat("###############################################################################\n\n")

# --- 2.1 Levene's Test for Homogeneity of Variance ---------------------------

cat("--- LEVENE'S TEST FOR HOMOGENEITY OF VARIANCE ---\n\n")

cat("Testing variance homogeneity across SITES (using sample of 50 sites):\n")
top_sites <- df_complete %>%
  count(site_code) %>%
  slice_max(n, n = 50) %>%
  pull(site_code)

df_levene_site <- df_complete %>% filter(site_code %in% top_sites)
levene_site <- car::leveneTest(bai ~ factor(site_code), data = df_levene_site)
cat(sprintf("  F = %.3f, df1 = %d, df2 = %d, p-value = %.2e\n", 
            levene_site$`F value`[1], levene_site$Df[1], levene_site$Df[2], levene_site$`Pr(>F)`[1]))
cat(sprintf("  Interpretation: %s\n\n",
            ifelse(levene_site$`Pr(>F)`[1] < 0.05, 
                   "Variances are HETEROGENEOUS across sites (p < 0.05)", 
                   "Variances are homogeneous across sites")))

cat("Testing variance homogeneity across SPECIES (top 10):\n")
top_species <- df_complete %>%
  count(species_code) %>%
  slice_max(n, n = 10) %>%
  pull(species_code)

df_levene_sp <- df_complete %>% filter(species_code %in% top_species)
levene_sp <- car::leveneTest(bai ~ factor(species_code), data = df_levene_sp)
cat(sprintf("  F = %.3f, df1 = %d, df2 = %d, p-value = %.2e\n", 
            levene_sp$`F value`[1], levene_sp$Df[1], levene_sp$Df[2], levene_sp$`Pr(>F)`[1]))
cat(sprintf("  Interpretation: %s\n\n",
            ifelse(levene_sp$`Pr(>F)`[1] < 0.05, 
                   "Variances are HETEROGENEOUS across species (p < 0.05)", 
                   "Variances are homogeneous across species")))

# --- 2.2 Intraclass Correlation Coefficient (ICC) ----------------------------

cat("--- INTRACLASS CORRELATION COEFFICIENT (ICC) ---\n")
cat("Fitting null model: bai ~ 1 + (1|site_code)\n")

model_null <- lmer(bai ~ 1 + (1 | site_code), data = df_complete, REML = TRUE)

icc_result <- performance::icc(model_null)
cat(sprintf("\nICC (site_code): %.4f\n", icc_result$ICC_adjusted))
cat(sprintf("  - %.1f%% of variance is BETWEEN sites\n", icc_result$ICC_adjusted * 100))
cat(sprintf("  - %.1f%% of variance is WITHIN sites\n", (1 - icc_result$ICC_adjusted) * 100))
cat(sprintf("\nInterpretation: %s\n",
            ifelse(icc_result$ICC_adjusted > 0.05,
                   "ICC > 0.05: Random effects for site ARE justified",
                   "ICC < 0.05: Random effects may not be necessary")))

cat("\nFitting model with species: bai ~ 1 + (1|site_code) + (1|species_code)\n")

model_null2 <- lmer(bai ~ 1 + (1 | site_code) + (1 | species_code), data = df_complete, REML = TRUE)

var_comp <- as.data.frame(VarCorr(model_null2))
total_var <- sum(var_comp$vcov) + sigma(model_null2)^2

cat("\nVariance Components:\n")
for (i in 1:nrow(var_comp)) {
  pct <- var_comp$vcov[i] / total_var * 100
  cat(sprintf("  %s: %.4f (%.1f%%)\n", var_comp$grp[i], var_comp$vcov[i], pct))
}
cat(sprintf("  Residual: %.4f (%.1f%%)\n", sigma(model_null2)^2, sigma(model_null2)^2 / total_var * 100))

cat("\n###############################################################################\n")
cat("# 3. BIVARIATE RELATIONSHIPS\n")
cat("###############################################################################\n\n")

# --- 3.1 Spearman Correlations -----------------------------------------------

cat("--- SPEARMAN CORRELATIONS WITH BAI ---\n")
cat("(Using Spearman due to non-normal distribution of BAI)\n\n")

predictors <- c("ndep_total_5yr", "ndep_total_3yr", "ndep_total",
                "temp_annual", "temp_summer", "temp_winter",
                "prec_annual", "prec_summer", "prec_winter", "lat")

cor_results <- data.frame(
  variable = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (pred in predictors) {
  if (pred %in% names(df_complete)) {
    test <- cor.test(df_complete$bai, df_complete[[pred]], method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      variable = pred,
      rho = test$estimate,
      p_value = test$p.value
    ))
  }
}

cor_results <- cor_results %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    strength = case_when(
      abs(rho) >= 0.5 ~ "Strong",
      abs(rho) >= 0.3 ~ "Moderate",
      abs(rho) >= 0.1 ~ "Weak",
      TRUE            ~ "Negligible"
    )
  ) %>%
  arrange(desc(abs(rho)))

print(cor_results)

write_csv(cor_results, path(out_dir_data, "02b_spearman_correlations.csv"))

# --- 3.2 Simple Linear Regressions -------------------------------------------

cat("\n--- SIMPLE LINEAR REGRESSIONS (BAI ~ predictor) ---\n\n")

lm_results <- data.frame(
  predictor = character(),
  estimate = numeric(),
  std_error = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  r_squared = numeric(),
  stringsAsFactors = FALSE
)

main_preds <- c("ndep_total_5yr", "temp_annual", "prec_annual", "lat")

for (pred in main_preds) {
  if (pred %in% names(df_complete)) {
    formula <- as.formula(paste("bai ~", pred))
    lm_fit <- lm(formula, data = df_complete)
    summ <- summary(lm_fit)
    
    lm_results <- rbind(lm_results, data.frame(
      predictor = pred,
      estimate = coef(summ)[2, 1],
      std_error = coef(summ)[2, 2],
      t_value = coef(summ)[2, 3],
      p_value = coef(summ)[2, 4],
      r_squared = summ$r.squared
    ))
  }
}

lm_results <- lm_results %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

print(lm_results)

write_csv(lm_results, path(out_dir_data, "02b_simple_regressions.csv"))

# --- 3.3 Non-linearity Check -------------------------------------------------

cat("\n--- NON-LINEARITY CHECK (Linear vs Quadratic) ---\n\n")

for (pred in main_preds) {
  if (pred %in% names(df_complete)) {
    lm_lin <- lm(as.formula(paste("bai ~", pred)), data = df_complete)
    lm_quad <- lm(as.formula(paste("bai ~", pred, "+ I(", pred, "^2)")), data = df_complete)
    
    anova_test <- anova(lm_lin, lm_quad)
    p_val <- anova_test$`Pr(>F)`[2]
    
    cat(sprintf("%s:\n", pred))
    cat(sprintf("  Linear R² = %.4f, Quadratic R² = %.4f\n", 
                summary(lm_lin)$r.squared, summary(lm_quad)$r.squared))
    cat(sprintf("  ANOVA p-value = %.2e\n", p_val))
    cat(sprintf("  Interpretation: %s\n\n",
                ifelse(p_val < 0.05, "Quadratic term IS significant - consider non-linear effect",
                       "Linear relationship is adequate")))
  }
}

# --- 3.4 Scatter Plots with Regression Lines ---------------------------------

cat("Creating: Bivariate scatter plots...\n")

set.seed(42)
df_sample <- df_complete %>% sample_n(min(30000, nrow(df_complete)))

p_ndep <- ggplot(df_sample, aes(x = ndep_total_5yr, y = bai)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_smooth(method = "loess", color = "blue", se = FALSE, linetype = "dashed") +
  scale_y_continuous(limits = c(0, quantile(df_sample$bai, 0.95))) +
  theme_minimal(base_size = 11) +
  labs(title = "BAI ~ N Deposition", x = "N dep. (5-yr)", y = "BAI")

p_temp <- ggplot(df_sample, aes(x = temp_annual, y = bai)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_smooth(method = "loess", color = "blue", se = FALSE, linetype = "dashed") +
  scale_y_continuous(limits = c(0, quantile(df_sample$bai, 0.95))) +
  theme_minimal(base_size = 11) +
  labs(title = "BAI ~ Temperature", x = "Temp. annual (°C)", y = "BAI")

p_prec <- ggplot(df_sample, aes(x = prec_annual, y = bai)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_smooth(method = "loess", color = "blue", se = FALSE, linetype = "dashed") +
  scale_y_continuous(limits = c(0, quantile(df_sample$bai, 0.95))) +
  theme_minimal(base_size = 11) +
  labs(title = "BAI ~ Precipitation", x = "Precip. annual (mm)", y = "BAI")

p_lat <- ggplot(df_sample, aes(x = lat, y = bai)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_smooth(method = "loess", color = "blue", se = FALSE, linetype = "dashed") +
  scale_y_continuous(limits = c(0, quantile(df_sample$bai, 0.95))) +
  theme_minimal(base_size = 11) +
  labs(title = "BAI ~ Latitude", x = "Latitude (°N)", y = "BAI")

p_bivar <- (p_ndep + p_temp) / (p_prec + p_lat)
ggsave(path(out_dir_fig, "02b_bivariate_relationships.png"), p_bivar, width = 12, height = 10, dpi = 300)

cat("\n###############################################################################\n")
cat("# 4. COLLINEARITY DIAGNOSTICS\n")
cat("###############################################################################\n\n")

# --- 4.1 Correlation Matrix of Predictors ------------------------------------

cat("--- CORRELATION MATRIX OF PREDICTORS ---\n\n")

pred_vars <- c("ndep_total_5yr", "temp_annual", "prec_annual", 
               "temp_summer", "temp_winter", "prec_summer", "prec_winter", "lat")
pred_vars <- pred_vars[pred_vars %in% names(df_complete)]

cor_matrix <- cor(df_complete[, pred_vars], use = "complete.obs", method = "spearman")
print(round(cor_matrix, 3))

cat("\nHighly correlated pairs (|r| > 0.7):\n")
for (i in 1:(ncol(cor_matrix) - 1)) {
  for (j in (i + 1):ncol(cor_matrix)) {
    if (abs(cor_matrix[i, j]) > 0.7) {
      cat(sprintf("  %s - %s: r = %.3f\n", 
                  rownames(cor_matrix)[i], colnames(cor_matrix)[j], cor_matrix[i, j]))
    }
  }
}

# --- 4.2 VIF Analysis --------------------------------------------------------

cat("\n--- VARIANCE INFLATION FACTOR (VIF) ---\n")
cat("Fitting model: bai ~ ndep_total_5yr + temp_annual + prec_annual + lat\n\n")

vif_model <- lm(bai ~ ndep_total_5yr + temp_annual + prec_annual + lat, data = df_complete)
vif_values <- car::vif(vif_model)

vif_df <- data.frame(
  variable = names(vif_values),
  VIF = as.numeric(vif_values)
) %>%
  mutate(
    interpretation = case_when(
      VIF < 5  ~ "OK",
      VIF < 10 ~ "Moderate collinearity",
      TRUE     ~ "HIGH collinearity - consider removing"
    )
  )

print(vif_df)

write_csv(vif_df, path(out_dir_data, "02b_vif_results.csv"))

cat("\nVIF interpretation:\n")
cat("  VIF < 5:  No collinearity concern\n")
cat("  VIF 5-10: Moderate collinearity\n")
cat("  VIF > 10: Severe collinearity - predictor should be removed or combined\n")

# --- 4.3 VIF with Seasonal Variables -----------------------------------------

cat("\n--- VIF WITH SEASONAL CLIMATE VARIABLES ---\n")
cat("Fitting model: bai ~ ndep_total_5yr + temp_summer + temp_winter + prec_summer + prec_winter\n\n")

if (all(c("temp_summer", "temp_winter", "prec_summer", "prec_winter") %in% names(df_complete))) {
  vif_model2 <- lm(bai ~ ndep_total_5yr + temp_summer + temp_winter + prec_summer + prec_winter, 
                   data = df_complete)
  vif_values2 <- car::vif(vif_model2)
  
  vif_df2 <- data.frame(
    variable = names(vif_values2),
    VIF = as.numeric(vif_values2)
  ) %>%
    mutate(
      interpretation = case_when(
        VIF < 5  ~ "OK",
        VIF < 10 ~ "Moderate collinearity",
        TRUE     ~ "HIGH collinearity"
      )
    )
  
  print(vif_df2)
}

# --- 4.4 Correlation Heatmap -------------------------------------------------

cat("\nCreating: Predictor correlation heatmap...\n")

cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(Var1 = Var1, Var2 = Var2, correlation = Freq)

p_cor <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman r") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Predictor Correlation Matrix", x = NULL, y = NULL)

ggsave(path(out_dir_fig, "02b_predictor_correlations.png"), p_cor, width = 10, height = 8, dpi = 300)

cat("\n###############################################################################\n")
cat("# 5. GROUP-LEVEL SAMPLE SIZE CHECKS\n")
cat("###############################################################################\n\n")

# --- 5.1 Observations per Site -----------------------------------------------

cat("--- OBSERVATIONS PER SITE ---\n")

site_counts <- df_complete %>%
  group_by(site_code) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  arrange(desc(n_obs))

cat(sprintf("Total sites: %d\n", nrow(site_counts)))
cat(sprintf("Mean obs per site: %.1f\n", mean(site_counts$n_obs)))
cat(sprintf("Median obs per site: %.1f\n", median(site_counts$n_obs)))
cat(sprintf("Min obs per site: %d\n", min(site_counts$n_obs)))
cat(sprintf("Max obs per site: %d\n", max(site_counts$n_obs)))

cat("\nDistribution of observations per site:\n")
site_bins <- site_counts %>%
  mutate(bin = cut(n_obs, breaks = c(0, 100, 500, 1000, 5000, 10000, Inf),
                   labels = c("1-100", "101-500", "501-1000", "1001-5000", "5001-10000", ">10000"))) %>%
  count(bin)
print(site_bins)

# --- 5.2 Observations per Species --------------------------------------------

cat("\n--- OBSERVATIONS PER SPECIES ---\n")

species_counts <- df_complete %>%
  group_by(species_code) %>%
  summarise(
    n_obs = n(),
    n_sites = n_distinct(site_code),
    .groups = "drop"
  ) %>%
  arrange(desc(n_obs))

cat(sprintf("Total species: %d\n", nrow(species_counts)))
cat(sprintf("Mean obs per species: %.1f\n", mean(species_counts$n_obs)))
cat(sprintf("Median obs per species: %.1f\n", median(species_counts$n_obs)))

cat("\nTop 15 species by observation count:\n")
print(head(species_counts, 15))

cat("\nSpecies with < 100 observations (potentially problematic for random effects):\n")
small_species <- species_counts %>% filter(n_obs < 100)
cat(sprintf("  %d species with < 100 obs\n", nrow(small_species)))

# --- 5.3 Balance Across N Deposition Gradient --------------------------------

cat("\n--- BALANCE ACROSS N DEPOSITION GRADIENT ---\n")

ndep_bins <- df_complete %>%
  mutate(ndep_bin = cut(ndep_total_5yr, breaks = 5)) %>%
  group_by(ndep_bin) %>%
  summarise(
    n_obs = n(),
    n_sites = n_distinct(site_code),
    mean_bai = mean(bai, na.rm = TRUE),
    .groups = "drop"
  )

print(ndep_bins)

# --- 5.4 Balance Across Decades ----------------------------------------------

cat("\n--- BALANCE ACROSS DECADES ---\n")

decade_counts <- df_complete %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade) %>%
  summarise(
    n_obs = n(),
    n_sites = n_distinct(site_code),
    n_species = n_distinct(species_code),
    .groups = "drop"
  )

print(decade_counts)

# --- 5.5 Sample Size Plots ---------------------------------------------------

cat("\nCreating: Sample size distribution plots...\n")

p_site_hist <- ggplot(site_counts, aes(x = n_obs)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  scale_x_log10() +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Observations per Site", 
       x = "Number of observations (log scale)", y = "Number of sites")

p_species_hist <- ggplot(species_counts, aes(x = n_obs)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  scale_x_log10() +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Observations per Species",
       x = "Number of observations (log scale)", y = "Number of species")

p_decade_bar <- ggplot(decade_counts, aes(x = factor(decade), y = n_sites)) +
  geom_col(fill = "coral", alpha = 0.8) +
  geom_text(aes(label = n_sites), vjust = -0.5, size = 3) +
  theme_minimal(base_size = 12) +
  labs(title = "Number of Sites per Decade", x = "Decade", y = "Number of sites")

p_sample <- (p_site_hist + p_species_hist) / p_decade_bar
ggsave(path(out_dir_fig, "02b_sample_size_distributions.png"), p_sample, width = 12, height = 10, dpi = 300)

# --- 6) Save Summary Tables --------------------------------------------------

write_csv(site_counts, path(out_dir_data, "02b_site_counts.csv"))
write_csv(species_counts, path(out_dir_data, "02b_species_counts.csv"))
write_csv(decade_counts, path(out_dir_data, "02b_decade_counts.csv"))

cat("\n###############################################################################\n")
cat("# 02b PRE-LMM DIAGNOSTICS COMPLETE\n")
cat("###############################################################################\n\n")

cat("=== SUMMARY OF KEY FINDINGS ===\n\n")

cat("1. DISTRIBUTION:\n")
cat(sprintf("   - Log-transformation %s normality\n", 
            ifelse(abs(skew_log) < abs(skew_raw), "IMPROVES", "does NOT improve")))

cat("\n2. VARIANCE STRUCTURE:\n")
cat(sprintf("   - ICC (site) = %.3f -> Random effects %s justified\n",
            icc_result$ICC_adjusted,
            ifelse(icc_result$ICC_adjusted > 0.05, "ARE", "may NOT be")))

cat("\n3. BIVARIATE RELATIONSHIPS:\n")
cat(sprintf("   - Strongest correlation: %s (rho = %.3f)\n",
            cor_results$variable[1], cor_results$rho[1]))

cat("\n4. COLLINEARITY:\n")
max_vif <- max(vif_df$VIF)
cat(sprintf("   - Maximum VIF = %.2f -> %s\n", max_vif,
            ifelse(max_vif < 5, "No collinearity issues", 
                   ifelse(max_vif < 10, "Moderate collinearity", "Severe collinearity"))))

cat("\n5. SAMPLE SIZE:\n")
cat(sprintf("   - %d sites, %d species\n", nrow(site_counts), nrow(species_counts)))
cat(sprintf("   - Median obs per site: %.0f\n", median(site_counts$n_obs)))

cat("\nPLOTS SAVED:\n")
plot_files <- dir_ls(out_dir_fig, glob = "*.png")
for (f in plot_files) {
  cat("  - ", path_file(f), "\n")
}

cat("\nDATA FILES SAVED:\n")
data_files <- dir_ls(out_dir_data, glob = "*.csv")
for (f in data_files) {
  cat("  - ", path_file(f), "\n")
}

cat("\n=== 02b COMPLETE ===\n")