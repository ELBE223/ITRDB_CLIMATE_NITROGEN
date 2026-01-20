###############################################################################
# 04a_prepare_panel_filters_and_z.R
#
# Filter panel (year + optional continent), drop sites by criteria,
# compute global z-standardized covariates, write dataset/parquet for LMM.
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, dplyr, fs, arrow, rlang)

CFG <- list(
  # ---- I/O ----
  PATH_IN   = here::here("output", "02_env_merge", "data", "02_panel_with_env.parquet"),
  OUT_DIR   = here::here("output", "04_a_LMM_prepare", "data"),
  OUT_NAME  = "04a_panel_z",
  WRITE_MODE = "dataset",     # "dataset" | "parquet"
  COMPRESSION = "zstd",
  
  # ---- Base filters ----
  YEAR_MIN = 1901,
  YEAR_MAX = NULL,            # NULL for no max
  CONTINENTS_KEEP = NULL,     # NULL = keep all continents
  CONTINENTS_DROP = NULL,     # optional
  CONTINENT_COL_CANDIDATES = c("continent", "continent_raw"),
  
  # ---- Site drop criteria ----
  SITE_ID_COL = "site_id",
  TREE_ID_CANDIDATES = c("tree_uid", "tree_id"),
  
  MIN_TREES_PER_SITE = 10,    # drop sites with < MIN_TREES_PER_SITE unique trees
  
  # >>> NEW: Switch for Time Series Length <<<
  MIN_TIME_SPAN      = 50,    # drop sites with (max_year - min_year + 1) < 50
  
  NDEPO_VAR = "N_depo1",       # column used for site-level N deposition filter
  NDEPO_SITE_STAT = "max",    # "mean" | "max"
  NDEPO_THRESHOLD = 15,        # drop sites with site_stat(NDEPO_VAR) > threshold
  NDEPO_NA_POLICY = "drop",    # "keep" | "drop" (if site_stat is NA)
  
  # ---- Schema ----
  REQUIRED_COLS = c("site_id","tree_id","tree_uid","year","genus","log_BAI_eps"),
  
  # ---- Z standardization ----
  Z_VARS = c(
    "age","year",
    "Temp_1","Nieder_1",
    "Temp_JJA","Temp_MAM","Temp_SON","Temp_DJF",
    "Nieder_JJA","Nieder_MAM","Nieder_SON","Nieder_DJF",
    "N_depo1","N_depo3","N_depo5"
  ),
  Z_SUFFIX = "_z",
  
  SD_METHOD = "sample",       # "sample" | "population"
  SD_POLICY = "floor",        # "na" | "zero" | "floor"
  SD_FLOOR  = 1e-8,
  
  STRICT = TRUE,
  
  # ---- Reporting ----
  SAVE_STATS = TRUE,
  STATS_FILE = "04a_z_stats.csv",
  SAVE_FILTER_REPORT = TRUE,
  FILTER_REPORT_FILE = "04a_filter_report.csv",
  
  # ---- Sanity check ----
  RUN_NDEPO_CHECK = TRUE
)

# ------------------------------ helpers --------------------------------------

pick_first_existing <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

apply_sd_policy <- function(sd, policy, floor = 1e-8) {
  if (policy == "zero")  return(ifelse(is.na(sd) || sd == 0, 1, sd))
  if (policy == "floor") return(ifelse(is.na(sd) || sd < floor, floor, sd))
  sd
}

safe_mean_sd <- function(x, method = "sample") {
  x <- suppressWarnings(as.numeric(x))
  n <- sum(!is.na(x))
  if (is.na(n) || n < 2) return(list(mean = NA_real_, sd = NA_real_))
  mu <- mean(x, na.rm = TRUE)
  if (method == "population") {
    sd <- sqrt(mean((x - mu)^2, na.rm = TRUE))
  } else {
    sd <- stats::sd(x, na.rm = TRUE)
  }
  list(mean = as.numeric(mu), sd = as.numeric(sd))
}

build_z_mutate_exprs <- function(vars, stats, suffix = "_z") {
  out <- list()
  for (v in vars) {
    mu <- stats[[v]]$mean
    sd <- stats[[v]]$sd
    out[[paste0(v, suffix)]] <- rlang::expr((!!rlang::sym(v) - !!mu) / !!sd)
  }
  out
}

safe_snapshot <- function(ds, stage, site_col, tree_col) {
  # Try Arrow-native distinct counts; fallback if needed.
  out <- tryCatch({
    ds %>%
      summarise(
        n_rows  = n(),
        n_sites = n_distinct(.data[[site_col]]),
        n_trees = n_distinct(.data[[tree_col]])
      ) %>% collect()
  }, error = function(e) {
    n_rows <- ds %>% summarise(n_rows = n()) %>% collect()
    n_sites <- ds %>% select(all_of(site_col)) %>% distinct() %>% summarise(n_sites = n()) %>% collect()
    n_trees <- ds %>% select(all_of(tree_col)) %>% distinct() %>% summarise(n_trees = n()) %>% collect()
    data.frame(
      n_rows  = n_rows$n_rows,
      n_sites = n_sites$n_sites,
      n_trees = n_trees$n_trees
    )
  })
  
  data.frame(
    stage = stage,
    n_rows  = as.numeric(out$n_rows[[1]]),
    n_sites = as.numeric(out$n_sites[[1]]),
    n_trees = as.numeric(out$n_trees[[1]])
  )
}

print_snapshot <- function(snap) {
  cat(
    sprintf(
      "%-28s  rows=%s  sites=%s  trees=%s\n",
      snap$stage,
      format(snap$n_rows,  big.mark = ","),
      format(snap$n_sites, big.mark = ","),
      format(snap$n_trees, big.mark = ",")
    )
  )
}

# ------------------------------ main -----------------------------------------

cat("=== 04a Prepare panel (filters + site drops + z) ===\n")
cat("Input: ", CFG$PATH_IN, "\n", sep = "")
stopifnot(fs::file_exists(CFG$PATH_IN))
fs::dir_create(CFG$OUT_DIR)

ds <- arrow::open_dataset(CFG$PATH_IN, format = "parquet")
vars <- names(ds)

continent_col <- pick_first_existing(CFG$CONTINENT_COL_CANDIDATES, vars)
if (is.na(continent_col)) stop("No continent column found (checked: ", paste(CFG$CONTINENT_COL_CANDIDATES, collapse = ", "), ").")

tree_col <- pick_first_existing(CFG$TREE_ID_CANDIDATES, vars)
if (is.na(tree_col)) stop("No tree id column found (checked: ", paste(CFG$TREE_ID_CANDIDATES, collapse = ", "), ").")

req_missing <- setdiff(unique(c(CFG$REQUIRED_COLS, CFG$Z_VARS, CFG$NDEPO_VAR)), vars)
if (CFG$STRICT && length(req_missing) > 0) stop("Missing required columns: ", paste(req_missing, collapse = ", "))

cat("\nContinent counts (raw):\n")
print(
  ds %>%
    group_by(.data[[continent_col]]) %>%
    summarise(n = n()) %>%
    collect() %>%
    arrange(desc(n))
)

report <- list()
snap0 <- safe_snapshot(ds, "00_raw", CFG$SITE_ID_COL, tree_col)
print_snapshot(snap0)
report[[length(report) + 1]] <- snap0

# ---- base filters ----
ds_f <- ds
if (!is.null(CFG$YEAR_MIN)) ds_f <- ds_f %>% filter(year >= CFG$YEAR_MIN)
if (!is.null(CFG$YEAR_MAX)) ds_f <- ds_f %>% filter(year <= CFG$YEAR_MAX)

if (!is.null(CFG$CONTINENTS_KEEP)) {
  ds_f <- ds_f %>% filter(.data[[continent_col]] %in% CFG$CONTINENTS_KEEP)
}
if (!is.null(CFG$CONTINENTS_DROP)) {
  ds_f <- ds_f %>% filter(!(.data[[continent_col]] %in% CFG$CONTINENTS_DROP))
}

snap1 <- safe_snapshot(ds_f, "01_after_year_continent", CFG$SITE_ID_COL, tree_col)
print_snapshot(snap1)
report[[length(report) + 1]] <- snap1

cat("\nContinent counts (filtered):\n")
print(
  ds_f %>%
    group_by(.data[[continent_col]]) %>%
    summarise(n = n()) %>%
    collect() %>%
    arrange(desc(n))
)

n_rows <- ds_f %>% summarise(n = n()) %>% collect() %>% pull(n)
if (n_rows < 2) stop("Too few rows after base filtering.")

# ---- drop sites with too few trees ----
cat("\nDropping sites with <", CFG$MIN_TREES_PER_SITE, " unique trees (tree col: ", tree_col, ")...\n", sep = "")

site_tree_counts <- tryCatch({
  ds_f %>%
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(n_trees = n_distinct(.data[[tree_col]])) %>%
    collect()
}, error = function(e) {
  ds_f %>%
    select(all_of(c(CFG$SITE_ID_COL, tree_col))) %>%
    distinct() %>%
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(n_trees = n(), .groups = "drop") %>%
    collect()
})

sites_keep_trees <- site_tree_counts %>%
  filter(n_trees >= CFG$MIN_TREES_PER_SITE) %>%
  pull(.data[[CFG$SITE_ID_COL]])

cat("Sites before (unique): ", format(nrow(site_tree_counts), big.mark = ","), "\n", sep = "")
cat("Sites dropped (<min trees): ", format(sum(site_tree_counts$n_trees < CFG$MIN_TREES_PER_SITE), big.mark = ","), "\n", sep = "")
cat("Sites kept: ", format(length(sites_keep_trees), big.mark = ","), "\n", sep = "")

ds_f <- ds_f %>% filter(.data[[CFG$SITE_ID_COL]] %in% sites_keep_trees)

snap2 <- safe_snapshot(ds_f, "02_after_min_trees_site", CFG$SITE_ID_COL, tree_col)
print_snapshot(snap2)
report[[length(report) + 1]] <- snap2

n_rows2 <- ds_f %>% summarise(n = n()) %>% collect() %>% pull(n)
if (n_rows2 < 2) stop("Too few rows after min-trees site filter.")

# =============================================================================
# >>> NEW BLOCK: Drop sites by Time Span (min 50 years) <<<
# =============================================================================
if (!is.null(CFG$MIN_TIME_SPAN) && CFG$MIN_TIME_SPAN > 0) {
  cat("\nDropping sites with time span < ", CFG$MIN_TIME_SPAN, " years...\n", sep = "")
  
  # Calculate min and max year per site
  site_span_df <- ds_f %>%
    select(all_of(c(CFG$SITE_ID_COL, "year"))) %>%
    distinct() %>% # unique site-year combos
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(
      y_min = min(year, na.rm = TRUE),
      y_max = max(year, na.rm = TRUE)
    ) %>%
    collect() %>%
    mutate(span = y_max - y_min + 1)
  
  sites_keep_span <- site_span_df %>%
    filter(span >= CFG$MIN_TIME_SPAN) %>%
    pull(.data[[CFG$SITE_ID_COL]])
  
  n_dropped_span <- sum(site_span_df$span < CFG$MIN_TIME_SPAN)
  
  cat("Sites before (unique): ", format(nrow(site_span_df), big.mark = ","), "\n", sep = "")
  cat("Sites dropped (<", CFG$MIN_TIME_SPAN, " yrs): ", format(n_dropped_span, big.mark = ","), "\n", sep = "")
  cat("Sites kept: ", format(length(sites_keep_span), big.mark = ","), "\n", sep = "")
  
  ds_f <- ds_f %>% filter(.data[[CFG$SITE_ID_COL]] %in% sites_keep_span)
  
  # Snapshot for the new step
  snap_span <- safe_snapshot(ds_f, "02b_after_min_timespan", CFG$SITE_ID_COL, tree_col)
  print_snapshot(snap_span)
  report[[length(report) + 1]] <- snap_span
  
  n_rows_span <- ds_f %>% summarise(n = n()) %>% collect() %>% pull(n)
  if (n_rows_span < 2) stop("Too few rows after time span filter.")
}
# =============================================================================

# ---- drop sites by N-depo (deduplicate by site-year to avoid tree weighting) ----
cat("\nDropping sites with ", CFG$NDEPO_SITE_STAT, "(", CFG$NDEPO_VAR, ") > ", CFG$NDEPO_THRESHOLD, "...\n", sep = "")
cat("Collecting only columns for N-depo site aggregation: ", CFG$SITE_ID_COL, ", year, ", CFG$NDEPO_VAR, "\n", sep = "")

ndepo_df <- ds_f %>%
  select(all_of(c(CFG$SITE_ID_COL, "year", CFG$NDEPO_VAR))) %>%
  collect() %>%
  distinct(.data[[CFG$SITE_ID_COL]], year, .keep_all = TRUE)

site_ndepo <- ndepo_df %>%
  group_by(.data[[CFG$SITE_ID_COL]]) %>%
  summarise(
    ndepo_site = if (CFG$NDEPO_SITE_STAT == "max") {
      max(.data[[CFG$NDEPO_VAR]], na.rm = TRUE)
    } else {
      mean(.data[[CFG$NDEPO_VAR]], na.rm = TRUE)
    },
    .groups = "drop"
  )

site_ndepo$ndepo_site[is.nan(site_ndepo$ndepo_site)] <- NA_real_

if (CFG$NDEPO_NA_POLICY == "drop") {
  sites_keep_ndepo <- site_ndepo %>%
    filter(!is.na(ndepo_site) & ndepo_site <= CFG$NDEPO_THRESHOLD) %>%
    pull(.data[[CFG$SITE_ID_COL]])
} else {
  sites_keep_ndepo <- site_ndepo %>%
    filter(is.na(ndepo_site) | ndepo_site <= CFG$NDEPO_THRESHOLD) %>%
    pull(.data[[CFG$SITE_ID_COL]])
}

cat("Sites before (unique): ", format(nrow(site_ndepo), big.mark = ","), "\n", sep = "")
cat("Sites dropped (N-depo > threshold): ", format(sum(site_ndepo$ndepo_site > CFG$NDEPO_THRESHOLD, na.rm = TRUE), big.mark = ","), "\n", sep = "")
cat("Sites with NA site N-depo: ", format(sum(is.na(site_ndepo$ndepo_site)), big.mark = ","), " (policy=", CFG$NDEPO_NA_POLICY, ")\n", sep = "")

cat("\nN-depo site stat summary (before filter):\n")
print(summary(site_ndepo$ndepo_site))

ds_f <- ds_f %>% filter(.data[[CFG$SITE_ID_COL]] %in% sites_keep_ndepo)

snap3 <- safe_snapshot(ds_f, "03_after_ndepo_site", CFG$SITE_ID_COL, tree_col)
print_snapshot(snap3)
report[[length(report) + 1]] <- snap3

n_rows3 <- ds_f %>% summarise(n = n()) %>% collect() %>% pull(n)
if (n_rows3 < 2) stop("Too few rows after N-depo site filter.")

# ---- sanity check for ndepo filter (optional) ----
if (isTRUE(CFG$RUN_NDEPO_CHECK)) {
  ndepo_check <- ds_f %>%
    select(all_of(c(CFG$SITE_ID_COL, "year", CFG$NDEPO_VAR))) %>%
    collect() %>%
    distinct(.data[[CFG$SITE_ID_COL]], year, .keep_all = TRUE) %>%
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(ndepo_site = mean(.data[[CFG$NDEPO_VAR]], na.rm = TRUE), .groups = "drop")
  
  max_nd <- max(ndepo_check$ndepo_site, na.rm = TRUE)
  cat("\nSanity check: max(site mean ", CFG$NDEPO_VAR, ") after filter = ", signif(max_nd, 6), "\n", sep = "")
  if (is.finite(max_nd) && max_nd > CFG$NDEPO_THRESHOLD + 1e-8) {
    warning("Sanity check failed: max N-depo still above threshold.")
  }
}

# ---- compute z-stats in R ----
cat("\nComputing z-stats in R (collecting only Z_VARS)...\n")
df_stats <- ds_f %>% select(any_of(CFG$Z_VARS)) %>% collect()

S <- list()
for (v in CFG$Z_VARS) {
  ms <- safe_mean_sd(df_stats[[v]], CFG$SD_METHOD)
  ms$sd <- apply_sd_policy(ms$sd, CFG$SD_POLICY, CFG$SD_FLOOR)
  S[[v]] <- ms
}

if (CFG$STRICT) {
  bad <- names(S)[!sapply(S, function(x) is.finite(x$mean) && is.finite(x$sd) && x$sd > 0)]
  if (length(bad) > 0) stop("Invalid mean/sd for: ", paste(bad, collapse = ", "))
}

stats_df <- data.frame(
  var  = names(S),
  mean = vapply(S, function(x) x$mean, numeric(1)),
  sd   = vapply(S, function(x) x$sd,   numeric(1))
)

if (isTRUE(CFG$SAVE_STATS)) {
  stats_path <- fs::path(CFG$OUT_DIR, CFG$STATS_FILE)
  utils::write.csv(stats_df, stats_path, row.names = FALSE)
  cat("Saved z-stats:", stats_path, "\n")
}

# ---- add z cols in Arrow query ----
mut_exprs <- build_z_mutate_exprs(CFG$Z_VARS, S, CFG$Z_SUFFIX)
ds_z <- ds_f %>% mutate(!!!mut_exprs)

# ---- write ----
if (CFG$WRITE_MODE == "dataset") {
  out_path <- fs::path(CFG$OUT_DIR, CFG$OUT_NAME)
  if (fs::dir_exists(out_path)) fs::dir_delete(out_path)
  arrow::write_dataset(ds_z, out_path, format = "parquet", compression = CFG$COMPRESSION)
  cat("\nSaved dataset:", out_path, "\n")
} else {
  out_path <- fs::path(CFG$OUT_DIR, paste0(CFG$OUT_NAME, ".parquet"))
  tbl <- ds_z %>% collect()
  arrow::write_parquet(tbl, out_path, compression = CFG$COMPRESSION)
  cat("\nSaved parquet:", out_path, "\n")
}

# ---- filter report ----
report_df <- dplyr::bind_rows(report)
cat("\nFilter report:\n")
print(report_df)

if (isTRUE(CFG$SAVE_FILTER_REPORT)) {
  rep_path <- fs::path(CFG$OUT_DIR, CFG$FILTER_REPORT_FILE)
  utils::write.csv(report_df, rep_path, row.names = FALSE)
  cat("Saved filter report:", rep_path, "\n")
}

cat("\nZ stats preview:\n")
print(stats_df)
cat("\n=== 04a DONE ===\n")
