###############################################################################
# 01d_QC_BAI_global_IMPROVED.R
#
# Description:
#   QC for high-latitude sites (>= 50°N).
#   Compares recomputed BAI values against stored values to verify accuracy.
###############################################################################

# --- 1) Packages & Setup -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  stringr,
  furrr,
  dplR,
  ggplot2,
  rnaturalearth,
  rnaturalearthdata,
  sf,
  tidyr,
  forcats
)

num_cores <- parallel::detectCores() - 1
plan(multisession, workers = num_cores)

TOLERANCE <- 1e-4

# --- 2) Directories & Paths --------------------------------------------------

step_01a <- "01a_download_ITRDB_RWL_global"
step_01b <- "01b_build_BAI_base_global"
step_01c <- "01c_compute_BAI_all_sites_global"

meta_path <- here::here("output", step_01b, "tables", "01b_ITRDB_site_metadata.csv")
bai_root  <- here::here("output", step_01c, "data", "bai")
rwl_root  <- here::here("output", step_01a, "data", "rwl")

base_out_dir <- here::here("output", "01d_QC_BAI_global")
dir_tables   <- path(base_out_dir, "tables")
dir_figures  <- path(base_out_dir, "figures")

dir_create(c(dir_tables, dir_figures))

cat("=== 01d QC BAI Global ===\n\n")
cat("Metadata: ", meta_path, "\n")
cat("BAI dir:  ", bai_root,  "\n")
cat("RWL dir:  ", rwl_root, "\n")
cat("Output:   ", base_out_dir, "\n\n")

# --- 3) Load & Filter Sites --------------------------------------------------

if (!file_exists(meta_path)) stop("Metadata not found. Did 01b run?")

raw_meta <- read_csv(meta_path, show_col_types = FALSE)
cat("Total sites in metadata: ", nrow(raw_meta), "\n")

# Clean species name function
clean_species_name <- function(x) {
  result <- character(length(x))
  
  for (i in seq_along(x)) {
    xi <- x[i]
    
    if (is.na(xi) || xi == "" || trimws(xi) == "") {
      result[i] <- NA_character_
      next
    }
    
    xi <- str_to_sentence(trimws(xi))
    words <- str_split(xi, "\\s+")[[1]]
    
    if (length(words) >= 2) {
      temp <- paste(words[1], words[2])
    } else if (length(words) == 1) {
      temp <- words[1]
    } else {
      result[i] <- NA_character_
      next
    }
    
    temp <- str_remove_all(temp, "[^A-Za-z ]")
    temp <- trimws(temp)
    
    if (temp == "" || nchar(temp) < 3) {
      result[i] <- NA_character_
    } else {
      result[i] <- temp
    }
  }
  
  return(result)
}

# Extract genus with fallback
extract_genus <- function(sci_name, species_code) {
  result <- character(length(sci_name))
  
  for (i in seq_along(sci_name)) {
    sci <- sci_name[i]
    code <- species_code[i]
    
    if (!is.na(sci) && sci != "") {
      genus <- word(sci, 1)
      if (!is.na(genus) && nchar(genus) >= 3) {
        result[i] <- genus
        next
      }
    }
    
    if (!is.na(code) && nchar(code) >= 4) {
      code_upper <- toupper(substr(code, 1, 2))
      genus_map <- c(
        "PI" = "Pinus",
        "PC" = "Picea", 
        "AB" = "Abies",
        "LA" = "Larix",
        "TS" = "Tsuga",
        "QU" = "Quercus",
        "BE" = "Betula",
        "JU" = "Juniperus",
        "PS" = "Pseudotsuga"
      )
      
      if (code_upper %in% names(genus_map)) {
        result[i] <- genus_map[code_upper]
        next
      }
    }
    
    result[i] <- "Unknown"
  }
  
  return(result)
}

# Filter and clean
filtered_meta <- raw_meta %>%
  filter(!is.na(lat), lat >= 50) %>%
  filter(has_rwl_tucson_local == TRUE) %>%
  mutate(
    scientific_name = clean_species_name(species_name),
    genus           = extract_genus(scientific_name, species_code),
    path_rwl_raw = file.path(rwl_root, region, rwl_tucson_file),
    path_bai_csv = file.path(bai_root, region, paste0(site_code, "_BAI.csv"))
  )

cat("Sites >= 50°N with Tucson RWL: ", nrow(filtered_meta), "\n\n")
if (nrow(filtered_meta) == 0) stop("No suitable sites above 50°N.")

# --- Species Summary ---------------------------------------------------------

cat("=== SPECIES SUMMARY ===\n")

species_summary <- filtered_meta %>%
  count(genus, scientific_name, species_code, sort = TRUE) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

write_csv(species_summary, path(dir_tables, "01d_species_list.csv"))

cat("Top 10 genera:\n")
species_summary %>%
  group_by(genus) %>%
  summarise(n_sites = sum(n), .groups = "drop") %>%
  arrange(desc(n_sites)) %>%
  head(10) %>%
  mutate(pct = round(100 * n_sites / nrow(filtered_meta), 1)) %>%
  print(n = 10)

cat("\nSpecies list saved.\n")

n_unknown <- sum(filtered_meta$genus == "Unknown")
if (n_unknown > 0) {
  cat(sprintf("\nWarning: %d sites (%.1f%%) have unknown genus\n", 
              n_unknown, 100 * n_unknown / nrow(filtered_meta)))
}

# --- Duplicate Check ---------------------------------------------------------

dup_site_codes <- raw_meta %>%
  count(site_code, sort = TRUE) %>%
  filter(n > 1)

write_csv(dup_site_codes, path(dir_tables, "01d_duplicate_site_codes.csv"))

if (nrow(dup_site_codes) > 0) {
  cat(sprintf("\nWarning: %d duplicate site_codes found\n", nrow(dup_site_codes)))
} else {
  cat("\nNo duplicate site_codes found.\n")
}

# --- 4) QC: Recompute & Compare BAI ------------------------------------------

cat("\n=== STARTING QC: Recompute BAI vs Stored ===\n")

check_site_qc <- function(site_code, path_rwl, path_bai) {
  
  make_result <- function(status, max_diff = NA_real_, detail = "") {
    tibble(
      site_code   = site_code, 
      qc_status   = status, 
      max_diff    = max_diff,
      qc_detail   = detail
    )
  }
  
  if (!file_exists(path_rwl)) {
    return(make_result("Missing RWL", detail = "RWL file not found"))
  }
  if (!file_exists(path_bai)) {
    return(make_result("Missing BAI", detail = "BAI file not found"))
  }
  
  rwl <- tryCatch(
    suppressWarnings(dplR::read.rwl(path_rwl)), 
    error = function(e) NULL
  )
  
  bai_stored <- tryCatch(
    read_csv(path_bai, show_col_types = FALSE), 
    error = function(e) NULL
  )
  
  if (is.null(rwl) || nrow(rwl) == 0 || ncol(rwl) == 0) {
    return(make_result("Read Error", detail = "Could not read RWL"))
  }
  
  if (is.null(bai_stored) || nrow(bai_stored) == 0) {
    return(make_result("Read Error", detail = "Could not read stored BAI"))
  }
  
  bai_mm2 <- tryCatch(
    suppressWarnings(dplR::bai.in(rwl)), 
    error = function(e) NULL
  )
  
  if (is.null(bai_mm2) || all(is.na(as.matrix(bai_mm2)))) {
    return(make_result("Calc Failed", detail = "BAI calculation failed"))
  }
  
  bai_cm2_new <- bai_mm2 / 100
  bai_cm2_new[bai_cm2_new > 100] <- NA
  
  years_new <- suppressWarnings(as.integer(rownames(bai_cm2_new)))
  
  if (any(is.na(years_new))) {
    return(make_result("Invalid Years", detail = "Could not parse year rownames"))
  }
  
  if (!"year" %in% names(bai_stored)) {
    return(make_result("Invalid Format", detail = "No year column in stored BAI"))
  }
  
  years_stored <- bai_stored$year
  
  common_years <- intersect(years_new, years_stored)
  
  if (length(common_years) == 0) {
    return(make_result("No Overlap", detail = "No overlapping years"))
  }
  
  series_new    <- colnames(bai_cm2_new)
  series_stored <- setdiff(colnames(bai_stored), "year")
  common_series <- intersect(series_new, series_stored)
  
  if (length(common_series) == 0) {
    return(make_result("No Overlap", detail = "No overlapping series"))
  }
  
  idx_new    <- match(common_years, years_new)
  idx_stored <- match(common_years, years_stored)
  
  mat_new <- as.matrix(bai_cm2_new[idx_new, common_series, drop = FALSE])
  mat_stored <- as.matrix(bai_stored[idx_stored, common_series, drop = FALSE])
  
  diff_mat <- abs(mat_new - mat_stored)
  
  if (all(is.na(diff_mat))) {
    return(make_result("All NA", detail = "All values are NA"))
  }
  
  max_diff <- max(diff_mat, na.rm = TRUE)
  
  if (max_diff < TOLERANCE) {
    status_str <- "PASS"
    detail_str <- sprintf("%d years, %d series", length(common_years), length(common_series))
  } else {
    status_str <- "FAIL"
    detail_str <- sprintf("Max diff: %.6f cm²", max_diff)
  }
  
  make_result(status_str, max_diff = max_diff, detail = detail_str)
}

qc_results <- filtered_meta %>%
  select(site_code, path_rwl = path_rwl_raw, path_bai = path_bai_csv) %>%
  future_pmap_dfr(
    check_site_qc,
    .progress = TRUE,
    .options  = furrr_options(seed = TRUE)
  )

final_qc_df <- filtered_meta %>%
  left_join(qc_results, by = "site_code")

write_csv(final_qc_df, path(dir_tables, "01d_QC_results.csv"))

# Summary statistics
cat("\n=== QC RESULTS SUMMARY ===\n")
qc_summary <- final_qc_df %>%
  count(qc_status, sort = TRUE) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

print(qc_summary)

cat("\nBy QC category:\n")
cat(sprintf("  ✓ PASS:           %4d (%.1f%%)\n", 
            sum(final_qc_df$qc_status == "PASS"), 
            100 * mean(final_qc_df$qc_status == "PASS")))
cat(sprintf("  ✗ FAIL:           %4d (%.1f%%)\n", 
            sum(final_qc_df$qc_status == "FAIL"), 
            100 * mean(final_qc_df$qc_status == "FAIL")))
cat(sprintf("  ⚠ Issues:         %4d (%.1f%%)\n", 
            sum(!final_qc_df$qc_status %in% c("PASS", "FAIL")), 
            100 * mean(!final_qc_df$qc_status %in% c("PASS", "FAIL"))))

failures <- final_qc_df %>%
  filter(qc_status != "PASS") %>%
  select(site_code, site_name, region, genus, species_code, 
         qc_status, qc_detail, max_diff)

write_csv(failures, path(dir_tables, "01d_QC_failures_detail.csv"))

if (nrow(failures) > 0) {
  cat(sprintf("\nDetailed failure report saved (%d sites)\n", nrow(failures)))
}

# --- 5) Map ------------------------------------------------------------------

cat("\n=== GENERATING MAP ===\n")

world <- ne_countries(scale = "medium", returnclass = "sf")

sites_sf <- final_qc_df %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

sites_sf <- sites_sf %>%
  mutate(
    qc_category = case_when(
      qc_status == "PASS"        ~ "PASS",
      qc_status == "FAIL"        ~ "FAIL",
      qc_status == "Missing BAI" ~ "No BAI",
      qc_status == "Missing RWL" ~ "No RWL",
      TRUE                       ~ "Other Issue"
    ),
    qc_category = factor(qc_category, 
                         levels = c("PASS", "FAIL", "No BAI", "No RWL", "Other Issue"))
  )

top_genera <- sites_sf %>%
  st_drop_geometry() %>%
  count(genus, sort = TRUE) %>%
  filter(genus != "Unknown") %>%
  slice(1:8) %>%
  pull(genus)

sites_sf <- sites_sf %>%
  mutate(
    genus_display = case_when(
      genus == "Unknown"          ~ "Unknown",
      genus %in% top_genera       ~ genus,
      TRUE                        ~ "Other"
    ),
    genus_display = factor(genus_display, 
                           levels = c(sort(top_genera), "Other", "Unknown"))
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
  "Other"      = "#7f7f7f",
  "Unknown"    = "#d62728"
)

shape_map <- c(
  "PASS"        = 16,
  "FAIL"        = 17,
  "No BAI"      = 4,
  "No RWL"      = 3,
  "Other Issue" = 1
)

size_map <- c(
  "PASS"        = 1.8,
  "FAIL"        = 2.5,
  "No BAI"      = 2.0,
  "No RWL"      = 2.0,
  "Other Issue" = 1.5
)

alpha_map <- c(
  "PASS"        = 0.8,
  "FAIL"        = 1.0,
  "No BAI"      = 0.7,
  "No RWL"      = 0.7,
  "Other Issue" = 0.6
)

n_pass <- sum(sites_sf$qc_category == "PASS")
n_fail <- sum(sites_sf$qc_category == "FAIL")
n_other <- nrow(sites_sf) - n_pass - n_fail

p_map <- ggplot() +
  geom_sf(data = world, fill = "#f0f0f0", color = "#d0d0d0", linewidth = 0.2) +
  geom_sf(
    data = sites_sf,
    aes(color = genus_display, 
        shape = qc_category, 
        size = qc_category,
        alpha = qc_category)
  ) +
  scale_color_manual(
    values = genus_colors, 
    name = "Genus",
    drop = FALSE
  ) +
  scale_shape_manual(
    values = shape_map, 
    name = "QC Status",
    drop = FALSE
  ) +
  scale_size_manual(
    values = size_map, 
    name = "QC Status",
    drop = FALSE
  ) +
  scale_alpha_manual(
    values = alpha_map,
    guide = "none"
  ) +
  coord_sf(ylim = c(45, 90), expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 60)) +
  scale_y_continuous(breaks = seq(50, 80, 10)) +
  theme_minimal() +
  theme(
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9),
    panel.grid.major  = element_line(color = "grey90", linetype = "dashed"),
    panel.background  = element_rect(fill = "aliceblue", color = NA),
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(size = 10, color = "grey30"),
    axis.text         = element_text(color = "grey40", size = 9),
    axis.title        = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2),
    size = guide_legend(order = 2)
  ) +
  labs(
    title    = "ITRDB High-Latitude Sites (≥50°N) - Quality Control",
    subtitle = sprintf(
      "Total: %d sites | PASS: %d (%.1f%%) | FAIL: %d (%.1f%%) | Other: %d (%.1f%%)",
      nrow(sites_sf),
      n_pass, 100 * n_pass / nrow(sites_sf),
      n_fail, 100 * n_fail / nrow(sites_sf),
      n_other, 100 * n_other / nrow(sites_sf)
    ),
    x = "Longitude",
    y = "Latitude"
  )

map_path <- path(dir_figures, "01d_map_50N_plus_improved.png")
ggsave(
  filename = map_path,
  plot     = p_map,
  width    = 14,
  height   = 8,
  bg       = "white",
  dpi      = 300
)

cat("Map saved to:", map_path, "\n")

# --- 6) Additional Summary Tables --------------------------------------------

region_summary <- final_qc_df %>%
  count(region, qc_status) %>%
  pivot_wider(names_from = qc_status, values_from = n, values_fill = 0)

write_csv(region_summary, path(dir_tables, "01d_QC_by_region.csv"))

genus_summary <- final_qc_df %>%
  count(genus, qc_status) %>%
  pivot_wider(names_from = qc_status, values_from = n, values_fill = 0) %>%
  arrange(desc(PASS))

write_csv(genus_summary, path(dir_tables, "01d_QC_by_genus.csv"))

cat("\n=== OUTPUT FILES ===\n")
cat("  QC results:       ", path(dir_tables, "01d_QC_results.csv"), "\n")
cat("  Failure details:  ", path(dir_tables, "01d_QC_failures_detail.csv"), "\n")
cat("  Species list:     ", path(dir_tables, "01d_species_list.csv"), "\n")
cat("  By region:        ", path(dir_tables, "01d_QC_by_region.csv"), "\n")
cat("  By genus:         ", path(dir_tables, "01d_QC_by_genus.csv"), "\n")
cat("  Map:              ", map_path, "\n")

cat("\n=== 01d QC BAI Global COMPLETE ===\n")
cat(sprintf("Successfully QC'd %d/%d sites (%.1f%%)\n", 
            n_pass, nrow(final_qc_df), 100 * n_pass / nrow(final_qc_df)))