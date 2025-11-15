############################################################
# 03_compute_BAI_all_sites_global.R
# Compute BAI for all ITRDB sites with ring-width data
# Uses:
#   - metadata/itrdb_global_site_metadata.csv   (from 02)
#   - classic *.rwl files in rwl/<region>/      (from 01)
############################################################

## --- 0) Packages ---------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(here, dplR)

## --- 1) Paths ------------------------------------------------------

base_dir       <- here::here("itrdb_global_example")
rwl_root       <- file.path(base_dir, "rwl")
meta_dir       <- file.path(base_dir, "metadata")
site_meta_path <- file.path(meta_dir, "itrdb_global_site_metadata.csv")
bai_root       <- file.path(base_dir, "bai")

if (!file.exists(site_meta_path)) {
  stop("Site metadata not found: ", site_meta_path,
       "\nRun 02_build_BAI_base_global.R first.")
}

dir.create(bai_root, recursive = TRUE, showWarnings = FALSE)

## --- 2) Read site metadata ----------------------------------------

meta_df <- read.csv(site_meta_path, stringsAsFactors = FALSE)

# Required columns (adapted to 'region' instead of 'region_path')
req_cols <- c(
  "region", "site_code", "classic_rwl_path", "has_classic_rwl",
  "species_code", "species_name", "lat", "lon", "elevation_m",
  "first_year", "last_year"
)

missing_cols <- setdiff(req_cols, names(meta_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ",
       paste(missing_cols, collapse = ", "))
}

# Keep only sites where we have a classic *.rwl locally
meta_df <- meta_df[meta_df$has_classic_rwl, , drop = FALSE]

if (nrow(meta_df) == 0L) {
  stop("No sites with classic RWL files available. Check scripts 01 and 02.")
}

# Optional: restrict to subset of sites for testing
site_subset <- NULL  # e.g. c("AUS116", "CAN001", "AK001")
if (!is.null(site_subset)) {
  meta_df <- meta_df[meta_df$site_code %in% site_subset, , drop = FALSE]
}

cat("Number of sites to process:", nrow(meta_df), "\n\n")

## --- 3) Prepare log file ------------------------------------------

log_path <- file.path(bai_root, "bai_global_log.txt")
log_con  <- file(log_path, open = "wt")

writeLines(paste0("Global BAI computation log - ", Sys.time()), log_con)
writeLines(paste("Number of sites:", nrow(meta_df)), log_con)
writeLines("", log_con)

# List to collect per-site summary stats
summary_list <- list()

## --- 4) Loop over sites and compute BAI ---------------------------

for (i in seq_len(nrow(meta_df))) {
  row <- meta_df[i, ]
  
  msg <- paste0("[", i, "/", nrow(meta_df), "] ",
                row$site_code, " - ", row$site_name,
                " (", row$region, ")")
  cat(msg, "\n")
  writeLines(msg, log_con)
  
  rwl_path <- row$classic_rwl_path
  
  # Safety: file exists?
  if (!file.exists(rwl_path)) {
    warn <- paste("  Missing RWL file:", rwl_path)
    cat(warn, "\n")
    writeLines(warn, log_con)
    next
  }
  
  ## --- Read RWL ----------------------------------------------------
  rwl <- try(read.rwl(rwl_path), silent = TRUE)
  if (inherits(rwl, "try-error") || !is.data.frame(rwl) ||
      nrow(rwl) == 0L || ncol(rwl) == 0L) {
    warn <- paste("  Could not read RWL for", row$site_code)
    cat(warn, "\n")
    writeLines(warn, log_con)
    next
  }
  
  ## --- Compute BAI -------------------------------------------------
  bai <- try(bai.in(rwl), silent = TRUE)
  if (inherits(bai, "try-error") || !is.data.frame(bai) ||
      nrow(bai) == 0L || ncol(bai) == 0L) {
    warn <- paste("  BAI computation failed for", row$site_code)
    cat(warn, "\n")
    writeLines(warn, log_con)
    next
  }
  
  # Check for all-NA BAI
  bai_mat <- as.matrix(bai)
  if (all(is.na(bai_mat))) {
    warn <- paste("  BAI all NA for", row$site_code, "- skipped")
    cat(warn, "\n")
    writeLines(warn, log_con)
    next
  }
  
  # Years from rownames
  years <- suppressWarnings(as.integer(rownames(bai)))
  if (all(is.na(years))) {
    warn <- paste("  Could not parse years for", row$site_code)
    cat(warn, "\n")
    writeLines(warn, log_con)
    next
  }
  
  bai_df <- data.frame(
    year = years,
    bai,
    check.names = FALSE
  )
  
  ## --- Save per-site BAI table ------------------------------------
  
  # Create region-specific BAI folder
  # 'region' may contain slashes (e.g. "northamerica/usa") -> nested folders
  bai_dir_region <- file.path(bai_root, row$region)
  dir.create(bai_dir_region, recursive = TRUE, showWarnings = FALSE)
  
  bai_file <- file.path(
    bai_dir_region,
    paste0(row$site_code, "_BAI.csv")
  )
  
  write.csv(bai_df, bai_file, row.names = FALSE)
  
  ok_msg <- paste("  BAI written to:", bai_file)
  cat(ok_msg, "\n")
  writeLines(ok_msg, log_con)
  
  ## --- Per-site summary statistics -------------------------------
  
  bai_vals <- as.numeric(bai_mat)
  bai_vals <- bai_vals[!is.na(bai_vals)]
  if (!length(bai_vals)) {
    warn <- paste("  No non-NA BAI values for", row$site_code,
                  " -> summary skipped")
    cat(warn, "\n")
    writeLines(warn, log_con)
    next
  }
  
  summary_list[[length(summary_list) + 1L]] <- data.frame(
    region       = row$region,
    site_code    = row$site_code,
    site_name    = row$site_name,
    species_code = row$species_code,
    species_name = row$species_name,
    lat          = row$lat,
    lon          = row$lon,
    elevation_m  = row$elevation_m,
    n_series     = ncol(bai),
    first_year   = min(bai_df$year, na.rm = TRUE),
    last_year    = max(bai_df$year, na.rm = TRUE),
    bai_min      = min(bai_vals, na.rm = TRUE),
    bai_q1       = as.numeric(quantile(bai_vals, 0.25, na.rm = TRUE)),
    bai_median   = median(bai_vals, na.rm = TRUE),
    bai_mean     = mean(bai_vals, na.rm = TRUE),
    bai_q3       = as.numeric(quantile(bai_vals, 0.75, na.rm = TRUE)),
    bai_max      = max(bai_vals, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

close(log_con)

## --- 5) Combine and save summary table ----------------------------

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  
  summary_path <- file.path(meta_dir, "itrdb_global_BAI_site_summary.csv")
  write.csv(summary_df, summary_path, row.names = FALSE)
  
  cat("\nGlobal BAI site summary written to:\n  ",
      normalizePath(summary_path), "\n")
  cat("Columns include:\n",
      "  region, site_code, site_name, species_code, species_name,\n",
      "  lat, lon, elevation_m, n_series, first_year, last_year,\n",
      "  bai_min, bai_q1, bai_median, bai_mean, bai_q3, bai_max\n\n")
} else {
  cat("\nNo sites produced valid BAI summaries. Check the log file:\n  ",
      normalizePath(log_path), "\n")
}

############################################################
# Result:
#  - Per-site BAI time series:
#      bai/<region>/<SITE_CODE>_BAI.csv
#      e.g. bai/northamerica/usa/AZ001_BAI.csv
#  - Global summary:
#      metadata/itrdb_global_BAI_site_summary.csv
#  - Log file:
#      bai/bai_global_log.txt
############################################################
