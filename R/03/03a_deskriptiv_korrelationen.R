###############################################################################
# 03a_deskriptiv_korrelationen.R  
#
# Criteria order:
# 1) Year range
# 2) Min trees per site (unique trees)
# 3) N-depo site filter (site STAT of N_depo1 <= threshold; dedup by site-year)
# 4) Min trees per genus (global)
# 5) Min years per site (>= 50 distinct years)
# 6) OPTIONAL: Drop sites with any NA in selected env vars (e.g., Temp_1, Nieder_1)
#
# Outputs:
#   output/03_a/qc/03a_deskriptiv_report.txt
#   output/03_a/data/03a_site_year.parquet
#   output/03_a/data/03a_cor_growth_env_spearman.csv
#   output/03_a/data/03a_cor_growth_env_by_continent_spearman.csv
#   output/03_a/data/03a_counts_genus_continent.csv
#   output/03_a/data/03a_cor_n_window_redundanz_spearman.csv   (optional)
#   output/03_a/data/03a_filter_report.csv
###############################################################################

# ------------------------------ Packages ---------------------------------------
pkgs <- c("here", "dplyr", "fs", "arrow", "tibble", "rlang")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

options(scipen = 999)

# ------------------------------ Config -----------------------------------------
CFG <- list(
  # Base filters
  YEAR_MIN = 1901L,
  YEAR_MAX = 2010L,
  CONTINENTS_KEEP = NULL,
  CONTINENTS_DROP = NULL,
  CONTINENT_COL_CANDIDATES = c("continent", "continent_raw"),
  
  SITE_ID_COL = "site_id",
  TREE_ID_CANDIDATES = c("tree_uid", "tree_id"),
  
  # Criteria
  MIN_TREES_PER_SITE = 10L,
  
  RUN_NDEPO_SITE_FILTER = TRUE,
  NDEPO_VAR = "N_depo1",
  NDEPO_SITE_STAT = "max",     # "mean" | "max"
  NDEPO_THRESHOLD = 15,
  NDEPO_NA_POLICY = "drop",     # "keep" | "drop" (usually irrelevant if 0 NA)
  
  RUN_GENUS_MIN_TREES = TRUE,
  MIN_TREES_PER_GENUS = 10L,
  
  RUN_MIN_YEARS_PER_SITE = TRUE,
  MIN_YEARS_PER_SITE = 50L,
  
  # OPTIONAL NA-site drop (this is the likely missing step to reach ~1051)
  RUN_DROP_SITES_WITH_NA_ENV = TRUE,
  NA_ENV_VARS = c("Temp_1", "Nieder_1"),   # change if needed
  NA_SITE_POLICY = "any",                  # "any" = drop site if any NA in any checked var
  
  # Expected final sites (optional check)
  EXPECTED_SITES_FINAL = 1051L,
  
  # Correlation filters
  MIN_N        = 5000L,
  MIN_ABS_RHO  = 0.10,
  
  MIN_N_CONT       = 2000L,
  MIN_ABS_RHO_CONT = 0.10,
  
  # Paths
  PATH_FINAL = here::here("output", "02_env_merge", "data", "02_panel_with_env.parquet"),
  
  OUT_DIR = here::here("output", "03_a"),
  OUT_QC  = here::here("output", "03_a", "qc"),
  OUT_DAT = here::here("output", "03_a", "data"),
  
  OUT_TXT = "03a_deskriptiv_report.txt",
  SAVE_FILTER_REPORT = TRUE,
  FILTER_REPORT_FILE = "03a_filter_report.csv"
)

fs::dir_create(CFG$OUT_DIR)
fs::dir_create(CFG$OUT_QC)
fs::dir_create(CFG$OUT_DAT)

stopifnot(fs::file_exists(CFG$PATH_FINAL))
path_txt <- fs::path(CFG$OUT_QC, CFG$OUT_TXT)

# ------------------------------ Helpers ----------------------------------------
pick_existing <- function(cols, candidates) candidates[candidates %in% cols]

pick_first_existing <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

desc_num <- function(x) {
  x_ok <- x[is.finite(x)]
  miss <- sum(!is.finite(x) | is.na(x))
  if (length(x_ok) == 0) {
    return(tibble(
      n = 0L, miss = miss,
      mean = NA_real_, sd = NA_real_,
      min = NA_real_, p05 = NA_real_, p50 = NA_real_, p95 = NA_real_, max = NA_real_
    ))
  }
  tibble(
    n    = length(x_ok),
    miss = miss,
    mean = mean(x_ok),
    sd   = stats::sd(x_ok),
    min  = min(x_ok),
    p05  = as.numeric(stats::quantile(x_ok, 0.05, names = FALSE, type = 7)),
    p50  = as.numeric(stats::quantile(x_ok, 0.50, names = FALSE, type = 7)),
    p95  = as.numeric(stats::quantile(x_ok, 0.95, names = FALSE, type = 7)),
    max  = max(x_ok)
  )
}

cor_spearman_1vN <- function(df, x, ys) {
  stopifnot(x %in% names(df))
  ys <- ys[ys %in% names(df)]
  out <- lapply(ys, function(y) {
    ok <- stats::complete.cases(df[[x]], df[[y]])
    n  <- sum(ok)
    r  <- suppressWarnings(stats::cor(df[[x]][ok], df[[y]][ok], method = "spearman"))
    tibble(var_x = x, var_y = y, rho = as.numeric(r), n = as.integer(n))
  })
  dplyr::bind_rows(out) %>% arrange(desc(abs(rho)))
}

pairwise_spearman_long <- function(df, vars) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) < 2) return(tibble(var_a = character(), var_b = character(), rho = numeric(), n = integer()))
  
  cmat <- suppressWarnings(stats::cor(df[, vars, drop = FALSE], use = "pairwise.complete.obs", method = "spearman"))
  out  <- tibble::as_tibble(as.table(cmat), .name_repair = "minimal")
  names(out) <- c("var_a", "var_b", "rho")
  out <- out %>% filter(var_a < var_b)
  
  nmat <- matrix(NA_integer_, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (i < j) {
        nmat[i, j] <- sum(stats::complete.cases(df[[vars[i]]], df[[vars[j]]]))
      }
    }
  }
  out$n <- as.integer(mapply(function(a, b) nmat[a, b], out$var_a, out$var_b))
  out %>% arrange(desc(abs(rho)))
}

fmt_tbl <- function(df, digits = 3) {
  df %>% mutate(across(where(is.numeric), ~round(.x, digits)))
}

safe_snapshot <- function(ds, stage, site_col, tree_col) {
  out <- tryCatch({
    ds %>%
      summarise(
        n_rows  = n(),
        n_sites = n_distinct(.data[[site_col]]),
        n_trees = n_distinct(.data[[tree_col]])
      ) %>% collect()
  }, error = function(e) {
    n_rows  <- ds %>% summarise(n_rows = n()) %>% collect()
    n_sites <- ds %>% select(all_of(site_col)) %>% distinct() %>% summarise(n_sites = n()) %>% collect()
    n_trees <- ds %>% select(all_of(tree_col)) %>% distinct() %>% summarise(n_trees = n()) %>% collect()
    data.frame(
      n_rows  = n_rows$n_rows,
      n_sites = n_sites$n_sites,
      n_trees = n_trees$n_trees
    )
  })
  
  data.frame(
    stage  = stage,
    n_rows = as.numeric(out$n_rows[[1]]),
    n_sites = as.numeric(out$n_sites[[1]]),
    n_trees = as.numeric(out$n_trees[[1]])
  )
}

print_snapshot <- function(snap) {
  cat(sprintf(
    "%-30s rows=%s  sites=%s  trees=%s\n",
    snap$stage,
    format(snap$n_rows,  big.mark = ","),
    format(snap$n_sites, big.mark = ","),
    format(snap$n_trees, big.mark = ",")
  ))
}

print_loss_vs_prev <- function(prev, curr) {
  cat("Sites lost vs previous: ", prev$n_sites - curr$n_sites, "\n", sep = "")
  cat("Trees lost vs previous: ", prev$n_trees - curr$n_trees, "\n", sep = "")
  cat("Rows  lost vs previous: ", prev$n_rows  - curr$n_rows,  "\n\n", sep = "")
}

# ------------------------------ Load -------------------------------------------
ds   <- arrow::open_dataset(CFG$PATH_FINAL, format = "parquet")
cols <- ds$schema$names

stopifnot(all(c("site_id", "year", "lat", "lon") %in% cols))
stopifnot("log_BAI_eps" %in% cols)
stopifnot("genus" %in% cols)

continent_col <- pick_first_existing(CFG$CONTINENT_COL_CANDIDATES, cols)
if (is.na(continent_col)) stop("No continent column found (checked: ", paste(CFG$CONTINENT_COL_CANDIDATES, collapse = ", "), ").")

tree_uid_col <- pick_first_existing(CFG$TREE_ID_CANDIDATES, cols)
if (is.na(tree_uid_col)) stop("No tree id column found (checked: ", paste(CFG$TREE_ID_CANDIDATES, collapse = ", "), ").")

# Env core selection
env_core <- c()
if ("Temp_1" %in% cols)      env_core <- c(env_core, "Temp_1") else if ("Temp_MAM" %in% cols) env_core <- c(env_core, "Temp_MAM")
if ("Nieder_1" %in% cols)    env_core <- c(env_core, "Nieder_1") else if ("Nieder_JJA" %in% cols) env_core <- c(env_core, "Nieder_JJA")

if (all(c("NHx_1", "NOy_1") %in% cols)) {
  env_core <- c(env_core, "NHx_1", "NOy_1")
} else if ("N_depo1" %in% cols) {
  env_core <- c(env_core, "N_depo1")
}
stopifnot(length(env_core) >= 2)

n_fam <- list(
  N_depo = pick_existing(cols, c("N_depo1", "N_depo3", "N_depo5")),
  NHx    = pick_existing(cols, c("NHx_1", "NHx_3", "NHx_5")),
  NOy    = pick_existing(cols, c("NOy_1", "NOy_3", "NOy_5"))
)

# ------------------------------ Report -----------------------------------------
sink(path_txt, split = TRUE)
cat("=== 03a DESCRIPTIVES + SPEARMAN ===\n")
cat("Input:  ", CFG$PATH_FINAL, "\n", sep = "")
cat("Year filter:  ", CFG$YEAR_MIN, "-", CFG$YEAR_MAX, "\n", sep = "")
cat("Continent col: ", continent_col, "\n", sep = "")
cat("Tree id col:   ", tree_uid_col, "\n", sep = "")
cat("Env (core): ", paste(env_core, collapse = ", "), "\n\n", sep = "")

# ------------------------------ Filters + losses -------------------------------
cat("---- FILTERS / SITE DROPS (with losses) ----\n")

report <- list()
snap0 <- safe_snapshot(ds, "00_raw", CFG$SITE_ID_COL, tree_uid_col)
print_snapshot(snap0)
report[[length(report) + 1]] <- snap0

# 1) Base filter: year (+ optional continent)
ds_filt <- ds %>%
  filter(year >= CFG$YEAR_MIN, year <= CFG$YEAR_MAX)

if (!is.null(CFG$CONTINENTS_KEEP)) {
  ds_filt <- ds_filt %>% filter(.data[[continent_col]] %in% CFG$CONTINENTS_KEEP)
}
if (!is.null(CFG$CONTINENTS_DROP)) {
  ds_filt <- ds_filt %>% filter(!(.data[[continent_col]] %in% CFG$CONTINENTS_DROP))
}

snap1 <- safe_snapshot(ds_filt, "01_after_year_continent", CFG$SITE_ID_COL, tree_uid_col)
print_snapshot(snap1)
report[[length(report) + 1]] <- snap1
print_loss_vs_prev(snap0, snap1)

# 2) Min trees per site
cat("Criterion: MIN_TREES_PER_SITE = ", CFG$MIN_TREES_PER_SITE, "\n", sep = "")

site_tree_counts <- tryCatch({
  ds_filt %>%
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(n_trees = n_distinct(.data[[tree_uid_col]])) %>%
    collect()
}, error = function(e) {
  ds_filt %>%
    select(all_of(c(CFG$SITE_ID_COL, tree_uid_col))) %>%
    distinct() %>%
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(n_trees = n(), .groups = "drop") %>%
    collect()
})

sites_keep_trees <- site_tree_counts %>%
  filter(n_trees >= CFG$MIN_TREES_PER_SITE) %>%
  pull(.data[[CFG$SITE_ID_COL]])

cat("Sites before: ", format(nrow(site_tree_counts), big.mark = ","), "\n", sep = "")
cat("Sites dropped (<min trees): ", format(sum(site_tree_counts$n_trees < CFG$MIN_TREES_PER_SITE), big.mark = ","), "\n", sep = "")
cat("Sites kept: ", format(length(sites_keep_trees), big.mark = ","), "\n", sep = "")

ds_filt <- ds_filt %>% filter(.data[[CFG$SITE_ID_COL]] %in% sites_keep_trees)

snap2 <- safe_snapshot(ds_filt, "02_after_min_trees_site", CFG$SITE_ID_COL, tree_uid_col)
print_snapshot(snap2)
report[[length(report) + 1]] <- snap2
print_loss_vs_prev(snap1, snap2)

# 3) N-depo site filter (dedup by site-year)
if (isTRUE(CFG$RUN_NDEPO_SITE_FILTER) && (CFG$NDEPO_VAR %in% cols)) {
  cat("Criterion: ", CFG$NDEPO_SITE_STAT, "(", CFG$NDEPO_VAR, ") <= ", CFG$NDEPO_THRESHOLD,
      "  (NA policy=", CFG$NDEPO_NA_POLICY, ")\n", sep = "")
  
  ndepo_df <- ds_filt %>%
    select(all_of(c(CFG$SITE_ID_COL, "year", CFG$NDEPO_VAR))) %>%
    collect() %>%
    distinct(.data[[CFG$SITE_ID_COL]], year, .keep_all = TRUE)
  
  # extra diagnostics
  site_ndepo_extra <- ndepo_df %>%
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(
      ndepo_site_mean = mean(.data[[CFG$NDEPO_VAR]], na.rm = TRUE),
      ndepo_site_max  = max(.data[[CFG$NDEPO_VAR]],  na.rm = TRUE),
      n_years         = n(),
      n_years_gt15    = sum(.data[[CFG$NDEPO_VAR]] > CFG$NDEPO_THRESHOLD, na.rm = TRUE),
      .groups = "drop"
    )
  
  site_ndepo_extra$ndepo_site_mean[is.nan(site_ndepo_extra$ndepo_site_mean)] <- NA_real_
  site_ndepo_extra$ndepo_site_max[!is.finite(site_ndepo_extra$ndepo_site_max)] <- NA_real_
  
  cat("Diagnostic: n_years_gt15 summary:\n")
  print(summary(site_ndepo_extra$n_years_gt15))
  cat("Diagnostic: n_years_gt15 quantiles (0,50,90,95,99,100%):\n")
  print(quantile(site_ndepo_extra$n_years_gt15, probs = c(0, .5, .9, .95, .99, 1), na.rm = TRUE))
  cat("\n")
  
  site_ndepo <- site_ndepo_extra %>%
    transmute(
      !!CFG$SITE_ID_COL := .data[[CFG$SITE_ID_COL]],
      ndepo_site = if (CFG$NDEPO_SITE_STAT == "max") ndepo_site_max else ndepo_site_mean
    )
  
  if (CFG$NDEPO_NA_POLICY == "drop") {
    sites_keep_ndepo <- site_ndepo %>%
      filter(!is.na(ndepo_site) & ndepo_site <= CFG$NDEPO_THRESHOLD) %>%
      pull(.data[[CFG$SITE_ID_COL]])
  } else {
    sites_keep_ndepo <- site_ndepo %>%
      filter(is.na(ndepo_site) | ndepo_site <= CFG$NDEPO_THRESHOLD) %>%
      pull(.data[[CFG$SITE_ID_COL]])
  }
  
  cat("Sites before: ", format(nrow(site_ndepo), big.mark = ","), "\n", sep = "")
  cat("Sites dropped (site_stat > threshold): ", format(sum(site_ndepo$ndepo_site > CFG$NDEPO_THRESHOLD, na.rm = TRUE), big.mark = ","), "\n", sep = "")
  cat("Sites with NA site_stat: ", format(sum(is.na(site_ndepo$ndepo_site)), big.mark = ","), "\n\n", sep = "")
  
  ds_filt <- ds_filt %>% filter(.data[[CFG$SITE_ID_COL]] %in% sites_keep_ndepo)
  
  snap3 <- safe_snapshot(ds_filt, "03_after_ndepo_site", CFG$SITE_ID_COL, tree_uid_col)
  print_snapshot(snap3)
  report[[length(report) + 1]] <- snap3
  print_loss_vs_prev(snap2, snap3)
  
} else {
  cat("N-depo site filter skipped.\n\n")
}

# 4) Min trees per genus
if (isTRUE(CFG$RUN_GENUS_MIN_TREES)) {
  cat("Criterion: MIN_TREES_PER_GENUS = ", CFG$MIN_TREES_PER_GENUS, "\n", sep = "")
  
  tree_genus <- ds_filt %>%
    select(all_of(c(tree_uid_col, "genus"))) %>%
    distinct() %>%
    collect()
  
  genus_counts <- tree_genus %>%
    group_by(genus) %>%
    summarise(n_trees = n_distinct(.data[[tree_uid_col]]), .groups = "drop")
  
  genera_keep <- genus_counts %>%
    filter(n_trees >= CFG$MIN_TREES_PER_GENUS) %>%
    pull(genus)
  
  cat("Genera before: ", nrow(genus_counts), "\n", sep = "")
  cat("Genera dropped (<min trees): ", sum(genus_counts$n_trees < CFG$MIN_TREES_PER_GENUS), "\n", sep = "")
  cat("Genera kept: ", length(genera_keep), "\n", sep = "")
  
  ds_filt <- ds_filt %>% filter(genus %in% genera_keep)
  
  snap4 <- safe_snapshot(ds_filt, "04_after_min_trees_genus", CFG$SITE_ID_COL, tree_uid_col)
  print_snapshot(snap4)
  report[[length(report) + 1]] <- snap4
  print_loss_vs_prev(report[[length(report) - 1]], snap4)
}

# 5) Min years per site
if (isTRUE(CFG$RUN_MIN_YEARS_PER_SITE)) {
  cat("Criterion: MIN_YEARS_PER_SITE = ", CFG$MIN_YEARS_PER_SITE, "\n", sep = "")
  
  site_years <- ds_filt %>%
    select(all_of(c(CFG$SITE_ID_COL, "year"))) %>%
    distinct() %>%
    collect() %>%
    group_by(.data[[CFG$SITE_ID_COL]]) %>%
    summarise(n_years = n_distinct(year), .groups = "drop")
  
  sites_keep_years <- site_years %>%
    filter(n_years >= CFG$MIN_YEARS_PER_SITE) %>%
    pull(.data[[CFG$SITE_ID_COL]])
  
  cat("Sites before: ", nrow(site_years), "\n", sep = "")
  cat("Sites dropped (<min years): ", sum(site_years$n_years < CFG$MIN_YEARS_PER_SITE), "\n", sep = "")
  cat("Sites kept: ", length(sites_keep_years), "\n", sep = "")
  
  ds_filt <- ds_filt %>% filter(.data[[CFG$SITE_ID_COL]] %in% sites_keep_years)
  
  snap5 <- safe_snapshot(ds_filt, "05_after_min_years_site", CFG$SITE_ID_COL, tree_uid_col)
  print_snapshot(snap5)
  report[[length(report) + 1]] <- snap5
  print_loss_vs_prev(report[[length(report) - 1]], snap5)
}

# 6) OPTIONAL: Drop sites with NA in env vars
if (isTRUE(CFG$RUN_DROP_SITES_WITH_NA_ENV)) {
  vars_na <- CFG$NA_ENV_VARS[CFG$NA_ENV_VARS %in% cols]
  if (length(vars_na) == 0) {
    cat("NA-site drop skipped (no NA_ENV_VARS found in data).\n\n")
  } else {
    cat("Criterion: DROP_SITES_WITH_NA_ENV (policy=", CFG$NA_SITE_POLICY, ") vars=",
        paste(vars_na, collapse = ", "), "\n", sep = "")
    
    na_df <- ds_filt %>%
      select(all_of(c(CFG$SITE_ID_COL, "year", vars_na))) %>%
      collect() %>%
      distinct(.data[[CFG$SITE_ID_COL]], year, .keep_all = TRUE)
    
    site_na <- na_df %>%
      group_by(.data[[CFG$SITE_ID_COL]]) %>%
      summarise(
        any_na = any(if_any(all_of(vars_na), is.na)),
        .groups = "drop"
      )
    
    n_sites_any_na <- sum(site_na$any_na, na.rm = TRUE)
    cat("Sites with ANY NA in checked vars: ", n_sites_any_na, "\n", sep = "")
    
    if (CFG$NA_SITE_POLICY == "any") {
      sites_keep_na <- site_na %>%
        filter(!any_na) %>%
        pull(.data[[CFG$SITE_ID_COL]])
    } else {
      # fallback: keep all
      sites_keep_na <- site_na %>% pull(.data[[CFG$SITE_ID_COL]])
    }
    
    ds_filt <- ds_filt %>% filter(.data[[CFG$SITE_ID_COL]] %in% sites_keep_na)
    
    snap6 <- safe_snapshot(ds_filt, "06_after_drop_na_sites", CFG$SITE_ID_COL, tree_uid_col)
    print_snapshot(snap6)
    report[[length(report) + 1]] <- snap6
    print_loss_vs_prev(report[[length(report) - 1]], snap6)
  }
}

# Filter report table
report_df <- dplyr::bind_rows(report) %>%
  mutate(
    sites_lost_vs_prev = dplyr::lag(n_sites) - n_sites,
    trees_lost_vs_prev = dplyr::lag(n_trees) - n_trees,
    rows_lost_vs_prev  = dplyr::lag(n_rows)  - n_rows
  )

cat("---- FILTER REPORT (stages) ----\n")
print(report_df)
cat("\n")

if (isTRUE(CFG$SAVE_FILTER_REPORT)) {
  rep_path <- fs::path(CFG$OUT_DAT, CFG$FILTER_REPORT_FILE)
  utils::write.csv(report_df, rep_path, row.names = FALSE)
  cat("Saved filter report: ", rep_path, "\n\n", sep = "")
}

# Expected sites check
final_sites <- tail(report_df$n_sites, 1)
if (!is.null(CFG$EXPECTED_SITES_FINAL) && is.finite(CFG$EXPECTED_SITES_FINAL)) {
  if (as.integer(final_sites) != as.integer(CFG$EXPECTED_SITES_FINAL)) {
    cat("WARNING: Final sites = ", final_sites, " but expected ", CFG$EXPECTED_SITES_FINAL, "\n\n", sep = "")
  } else {
    cat("OK: Final sites match expected = ", CFG$EXPECTED_SITES_FINAL, "\n\n", sep = "")
  }
}

# ------------------------------ Basis (after filters) ---------------------------
n_rows <- ds_filt %>% summarise(n = n()) %>% collect() %>% pull(n)

site_meta <- ds_filt %>%
  select(site_id,
         continent = all_of(continent_col),
         lat, lon,
         any_of("elevation_m")) %>%
  distinct() %>%
  collect()

tree_meta <- ds_filt %>%
  select(!!rlang::sym(tree_uid_col),
         genus,
         continent = all_of(continent_col),
         site_id) %>%
  distinct() %>%
  collect()

site_year_keys <- ds_filt %>%
  select(site_id, year) %>%
  distinct() %>%
  collect()

yr_rng <- ds_filt %>%
  summarise(ymin = min(year, na.rm = TRUE), ymax = max(year, na.rm = TRUE)) %>%
  collect()

cat("---- Basis (AFTER FILTERS) ----\n")
cat("Rows (Panel):          ", n_rows, "\n", sep = "")
cat("Sites:                 ", nrow(site_meta), "\n", sep = "")
cat("Trees (distinct):      ", nrow(tree_meta), "\n", sep = "")
cat("Site-Years (distinct): ", nrow(site_year_keys), "\n", sep = "")
cat("Year-Range (Panel):    ", yr_rng$ymin, "-", yr_rng$ymax, "\n\n", sep = "")

# ------------------------------ Aggregation: Site-Year --------------------------
site_year <- ds_filt %>%
  select(site_id, year, log_BAI_eps, all_of(env_core)) %>%
  group_by(site_id, year) %>%
  summarise(
    n_rows = n(),
    log_BAI_eps_mean = mean(log_BAI_eps, na.rm = TRUE),
    across(all_of(env_core), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  collect() %>%
  left_join(site_meta, by = "site_id") %>%
  mutate(
    year = as.integer(year),
    continent = as.character(continent)
  )

arrow::write_parquet(
  site_year,
  fs::path(CFG$OUT_DAT, "03a_site_year.parquet"),
  compression = "zstd"
)

cat("---- Site-Year ----\n")
cat("Rows (Site-Year): ", nrow(site_year), "\n\n", sep = "")

# ------------------------------ Missingness + Descriptives ----------------------
vars_report <- c("log_BAI_eps_mean", env_core)

cat("---- Missingness (Site-Year) ----\n")
miss <- tibble(
  variable = vars_report,
  miss = vapply(vars_report, function(v) sum(!is.finite(site_year[[v]]) | is.na(site_year[[v]])), integer(1)),
  n = nrow(site_year)
)
print(miss)
cat("\n")

cat("---- Descriptive stats (Site-Year) ----\n")
desc_tbl <- bind_rows(lapply(vars_report, function(v) {
  d <- desc_num(site_year[[v]])
  d$variable <- v
  d
})) %>%
  select(variable, everything())
print(fmt_tbl(desc_tbl, digits = 3))
cat("\n")

# ------------------------------ Correlations: Growth vs Env ---------------------
cor_main <- cor_spearman_1vN(site_year, x = "log_BAI_eps_mean", ys = env_core) %>%
  filter(n >= CFG$MIN_N, abs(rho) >= CFG$MIN_ABS_RHO)

write.csv(
  cor_main,
  fs::path(CFG$OUT_DAT, "03a_cor_growth_env_spearman.csv"),
  row.names = FALSE
)

cat("---- Spearman: Growth vs Environment (filtered) ----\n")
cat("Filter: n >= ", CFG$MIN_N, " and |rho| >= ", CFG$MIN_ABS_RHO, "\n", sep = "")
if (nrow(cor_main) == 0) {
  cat("No hits after filter.\n\n")
} else {
  print(fmt_tbl(head(cor_main, 10), digits = 3))
  cat("\n")
}

cor_by_cont <- site_year %>%
  filter(!is.na(continent), continent != "") %>%
  group_by(continent) %>%
  group_modify(~cor_spearman_1vN(.x, x = "log_BAI_eps_mean", ys = env_core)) %>%
  ungroup() %>%
  filter(n >= CFG$MIN_N_CONT, abs(rho) >= CFG$MIN_ABS_RHO_CONT)

write.csv(
  cor_by_cont,
  fs::path(CFG$OUT_DAT, "03a_cor_growth_env_by_continent_spearman.csv"),
  row.names = FALSE
)

cat("---- Spearman: Growth vs Env by continent (Top 5) ----\n")
cat("Filter: n >= ", CFG$MIN_N_CONT, " and |rho| >= ", CFG$MIN_ABS_RHO_CONT, "\n", sep = "")
if (nrow(cor_by_cont) == 0) {
  cat("No hits after filter.\n\n")
} else {
  top_cont <- cor_by_cont %>%
    group_by(continent) %>%
    arrange(desc(abs(rho)), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()
  print(fmt_tbl(top_cont, digits = 3))
  cat("\n")
}

# ------------------------------ Optional: N-window redundancy -------------------
all_n_windows <- unique(unlist(n_fam))
if (length(all_n_windows) >= 2) {
  n_red <- pairwise_spearman_long(site_year, all_n_windows) %>%
    filter(n >= CFG$MIN_N)
  
  if (nrow(n_red) > 0) {
    write.csv(
      n_red,
      fs::path(CFG$OUT_DAT, "03a_cor_n_window_redundanz_spearman.csv"),
      row.names = FALSE
    )
    cat("---- N-window redundancy (Spearman, Top 10 |rho|) ----\n")
    print(fmt_tbl(head(n_red, 10), digits = 3))
    cat("\n")
  }
}

# ------------------------------ Counts Genus x Continent ------------------------
counts_genus_cont <- tree_meta %>%
  mutate(genus = as.character(genus), continent = as.character(continent)) %>%
  group_by(continent, genus) %>%
  summarise(
    n_trees = n(),
    n_sites = n_distinct(site_id),
    .groups = "drop"
  ) %>%
  arrange(continent, desc(n_trees))

write.csv(
  counts_genus_cont,
  fs::path(CFG$OUT_DAT, "03a_counts_genus_continent.csv"),
  row.names = FALSE
)

cat("---- Genus x Continent (Top per continent) ----\n")
top_gen <- counts_genus_cont %>%
  group_by(continent) %>%
  slice_head(n = 3) %>%
  ungroup()
print(top_gen)
cat("\n")

cat("=== 03a DONE ===\n")
cat("Report: ", path_txt, "\n", sep = "")
cat("Data:   ", CFG$OUT_DAT, "\n", sep = "")

sink()
message("03a fertig: ", path_txt)
