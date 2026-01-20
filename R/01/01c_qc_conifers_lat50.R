###############################################################################
# 01c_qc_conifers_lat50.R
#
# Goal:
#   1) Use 01a metadata (LAT>=50 already) and classify Conifer vs non-Conifer.
#   2) Filter for Conifers only.
#   3) Summaries by Continent and Genus.
#   4) Filter: keep only Genera with > 15 sites (GLOBAL threshold).
#   5) Save key tables + print to console.
#   6) FINAL GEO QC (NO DOWNLOADS): sanity-check continent vs site_id prefix and
#      coarse lon bounding boxes; flag "Asia but looks like NorthAmerica".
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, dplyr, readr, stringr, fs, tibble)

# ----------------------------- CONFIG ----------------------------------------

CFG <- list(
  LAT_MIN_KEEP = 50,
  MIN_SITES_GENUS_GLOBAL = 15L,  # keep genera with > 15 sites (global)
  N_DEBUG_OTHER = 50L,
  N_DEBUG_GEO   = 20L
)

# ----------------------------- PATHS -----------------------------------------

path_in <- here::here("output", "01a_analyse_ITRDB", "tables", "01a_ITRDB_site_metadata_lat50.csv")
if (!fs::file_exists(path_in)) stop("Input file not found: ", path_in)

out_dir <- here::here("output", "01c_qc_descriptives_lat50_conifers")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

path_out_global_types <- fs::path(out_dir, "01c_summary_global_types.csv")
path_out_cont_summary <- fs::path(out_dir, "01c_summary_conifers_by_continent.csv")
path_out_cont_table   <- fs::path(out_dir, "01c_conifers_by_continent_genus_filtered.csv")
path_out_genus_table  <- fs::path(out_dir, "01c_conifers_global_genus_filtered.csv")
path_out_other_debug  <- fs::path(out_dir, "01c_other_unknown_genera_examples.csv")
path_out_geo_flags    <- fs::path(out_dir, "01c_geo_qc_flags.csv")

cat("=== 01c QC + Descriptives (LAT>=50, Conifers) ===\n")
cat("Input:  ", path_in, "\n", sep = "")
cat("Output: ", out_dir, "\n\n", sep = "")

# ----------------------------- LOAD ------------------------------------------

df_raw <- readr::read_csv(path_in, show_col_types = FALSE)

req_cols <- c("site_id", "continent", "lat", "lon", "species_name")
missing_cols <- setdiff(req_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------- PREP ------------------------------------------

conifer_genera_ext <- c(
  "Abies","Cedrus","Chamaecyparis","Cryptomeria","Cupressus","Juniperus","Larix",
  "Metasequoia","Picea","Pinus","Pseudotsuga","Sequoia","Sequoiadendron",
  "Taxodium","Taxus","Thuja","Tsuga","Podocarpus","Araucaria","Agathis"
)

broadleaf_genera_hint <- c(
  "Alnus","Betula","Fagus","Fraxinus","Populus","Quercus","Salix","Tilia","Acer","Ulmus"
)

df_prep <- df_raw %>%
  mutate(
    genus_from_name = str_to_title(str_extract(species_name, "^[A-Za-z]+")),
    genus = if ("genus" %in% names(.)) coalesce(.data$genus, genus_from_name) else genus_from_name,
    
    leaf_type_norm = if ("leaf_type" %in% names(.)) str_to_lower(.data$leaf_type) else NA_character_,
    
    Type = case_when(
      leaf_type_norm == "conifer" ~ "Conifer",
      is.na(leaf_type_norm) & genus %in% conifer_genera_ext ~ "Conifer",
      genus %in% broadleaf_genera_hint ~ "Broadleaf",
      TRUE ~ "Other"
    ),
    
    leaf_type_safe = if ("leaf_type" %in% names(.)) as.character(.data$leaf_type) else NA_character_
  )

cat("Rows loaded (LAT>=50 table): ", nrow(df_prep), "\n", sep = "")
cat("Missing species_name:        ", sum(is.na(df_prep$species_name) | df_prep$species_name == "", na.rm = TRUE), "\n", sep = "")
cat("Missing genus (after prep):  ", sum(is.na(df_prep$genus) | df_prep$genus == "", na.rm = TRUE), "\n\n", sep = "")

# ----------------------------- FINAL GEO QC (NO DOWNLOADS) --------------------

to360 <- function(x) ifelse(is.na(x), NA_real_, ifelse(x < 0, x + 360, x))

df_qc <- df_prep %>%
  mutate(
    lon360 = to360(lon),
    prefix = stringr::str_extract(site_id, "^[a-z]+"),
    continent_from_prefix = dplyr::case_when(
      prefix == "northamerica" ~ "NorthAmerica",
      prefix == "southamerica" ~ "SouthAmerica",
      prefix == "europe"       ~ "Europe",
      prefix == "asia"         ~ "Asia",
      prefix == "africa"       ~ "Africa",
      prefix == "australia"    ~ "Australia",
      TRUE ~ NA_character_
    ),
    
    # Prefix mismatch (allow asia__rus* -> Europe as valid after 01a fix)
    flag_prefix_mismatch = !is.na(continent_from_prefix) &
      continent != continent_from_prefix &
      !(continent == "Europe" & stringr::str_detect(site_id, "^asia__rus")),
    
    # Coarse lon sanity checks using lon360:
    # Europe approx: lon360 in [330..360] U [0..75]
    # Asia approx:   lon360 in [20..200]
    # NorthAmerica:  lon360 in [170..330]
    flag_geo_bbox = dplyr::case_when(
      continent == "Europe" &
        !(lon360 >= 330 | lon360 <= 75) ~ TRUE,
      continent == "Asia" &
        !(lon360 >= 20 & lon360 <= 200) ~ TRUE,
      continent == "NorthAmerica" &
        !(lon360 >= 170 & lon360 <= 330) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Special: labelled Asia but looks like NorthAmerica (lon360 230..330)
    flag_asia_looks_na = (continent == "Asia") &
      !is.na(lon360) & (lon360 >= 230 & lon360 <= 330)
  )

cat("=== FINAL GEO QC (no downloads) ===\n")
cat("Prefix mismatches (excluding asia__rus -> Europe): ", sum(df_qc$flag_prefix_mismatch, na.rm = TRUE), "\n", sep = "")
cat("Geo bbox suspicious (very coarse):                ", sum(df_qc$flag_geo_bbox, na.rm = TRUE), "\n", sep = "")
cat("Asia but looks like NorthAmerica (lon360 230..330):", sum(df_qc$flag_asia_looks_na, na.rm = TRUE), "\n\n", sep = "")

if (any(df_qc$flag_prefix_mismatch, na.rm = TRUE)) {
  cat("Top prefix-mismatch examples (up to ", CFG$N_DEBUG_GEO, "):\n", sep = "")
  print(
    df_qc %>%
      filter(flag_prefix_mismatch) %>%
      select(site_id, continent, continent_from_prefix, lat, lon, lon360, species_name) %>%
      slice_head(n = CFG$N_DEBUG_GEO),
    n = CFG$N_DEBUG_GEO
  )
  cat("\n")
}

if (any(df_qc$flag_geo_bbox, na.rm = TRUE)) {
  cat("Top geo-bbox suspicious examples (up to ", CFG$N_DEBUG_GEO, "):\n", sep = "")
  print(
    df_qc %>%
      filter(flag_geo_bbox) %>%
      select(site_id, continent, lat, lon, lon360, species_name) %>%
      slice_head(n = CFG$N_DEBUG_GEO),
    n = CFG$N_DEBUG_GEO
  )
  cat("\n")
}

if (any(df_qc$flag_asia_looks_na, na.rm = TRUE)) {
  cat("Asia-but-looks-NorthAmerica examples (up to ", CFG$N_DEBUG_GEO, "):\n", sep = "")
  print(
    df_qc %>%
      filter(flag_asia_looks_na) %>%
      select(site_id, continent, lat, lon, lon360, species_name) %>%
      slice_head(n = CFG$N_DEBUG_GEO),
    n = CFG$N_DEBUG_GEO
  )
  cat("\n")
}

df_qc %>%
  filter(flag_prefix_mismatch | flag_geo_bbox | flag_asia_looks_na) %>%
  select(
    site_id, continent, continent_from_prefix, lat, lon, lon360, species_name,
    flag_prefix_mismatch, flag_geo_bbox, flag_asia_looks_na
  ) %>%
  readr::write_csv(path_out_geo_flags)

cat("Saved geo QC flags: ", path_out_geo_flags, "\n\n", sep = "")

# ----------------------------- 1) GLOBAL SUMMARY -----------------------------

cat("=== 1) GLOBAL SUMMARY: ALL TYPES ===\n")
summary_global <- df_prep %>%
  count(Type, name = "N_Sites") %>%
  mutate(Percent = round(N_Sites / sum(N_Sites) * 100, 1)) %>%
  arrange(desc(N_Sites))

print(summary_global)

readr::write_csv(summary_global, path_out_global_types)
cat("\nSaved: ", path_out_global_types, "\n\n", sep = "")

# ----------------------------- 2) CONIFERS ONLY ------------------------------

df_conifers <- df_prep %>% filter(Type == "Conifer")

cat("=== 2) FILTER: CONIFERS ONLY ===\n")
cat("Conifer sites: ", nrow(df_conifers), "\n\n", sep = "")

# ----------------------------- 3) CONIFERS BY CONTINENT ----------------------

cat("=== 3) CONIFERS: DISTRIBUTION BY CONTINENT ===\n")
summary_continents <- df_conifers %>%
  count(continent, name = "N_Sites") %>%
  mutate(
    Total_Conifers = sum(N_Sites),
    Percent_Global = round((N_Sites / Total_Conifers) * 100, 1)
  ) %>%
  select(Continent = continent, N_Sites, Percent_Global) %>%
  arrange(desc(N_Sites))

print(summary_continents)

readr::write_csv(summary_continents, path_out_cont_summary)
cat("\nSaved: ", path_out_cont_summary, "\n\n", sep = "")

# ----------------------------- 4) GENUS FILTER (>15 GLOBAL) -------------------

cat(
  "=== 4) GENUS COUNTS (CONIFERS) + FILTER (GLOBAL > ",
  CFG$MIN_SITES_GENUS_GLOBAL,
  " SITES) ===\n",
  sep = ""
)

genus_counts_global <- df_conifers %>%
  count(genus, name = "N_Sites") %>%
  arrange(desc(N_Sites))

cat("Top genera (before filter):\n")
print(head(genus_counts_global, 15))

keep_genera <- genus_counts_global %>%
  filter(N_Sites > CFG$MIN_SITES_GENUS_GLOBAL) %>%
  pull(genus)

cat("\nGenera kept (global): ", length(keep_genera), "\n", sep = "")
cat("Kept genera: ", paste(keep_genera, collapse = ", "), "\n\n", sep = "")

# ----------------------------- 5) TABLE A: CONTINENT x GENUS ------------------

cat("=== 5) TABLE A: GENUS PER CONTINENT (only kept genera) ===\n")

final_table_cont <- df_conifers %>%
  filter(genus %in% keep_genera) %>%
  count(continent, genus, name = "N_Sites") %>%
  group_by(continent) %>%
  mutate(
    Total_in_Continent = sum(N_Sites),
    Percent_in_Continent = round((N_Sites / Total_in_Continent) * 100, 1)
  ) %>%
  ungroup() %>%
  select(
    Continent = continent,
    Genus = genus,
    N_Sites,
    Percent = Percent_in_Continent
  ) %>%
  arrange(Continent, desc(N_Sites))

print(final_table_cont, n = Inf)

readr::write_csv(final_table_cont, path_out_cont_table)
cat("\nSaved: ", path_out_cont_table, "\n\n", sep = "")

# ----------------------------- 6) TABLE B: GLOBAL GENUS -----------------------

cat("=== 6) TABLE B: GLOBAL GENUS (only kept genera) ===\n")

final_table_genus <- genus_counts_global %>%
  mutate(
    Total_Conifers = sum(N_Sites),
    Percent_Global = round((N_Sites / Total_Conifers) * 100, 1)
  ) %>%
  filter(genus %in% keep_genera) %>%
  select(
    Genus = genus,
    N_Sites,
    Percent_Global
  ) %>%
  arrange(desc(N_Sites))

print(final_table_genus, n = Inf)

readr::write_csv(final_table_genus, path_out_genus_table)
cat("\nSaved: ", path_out_genus_table, "\n\n", sep = "")

# ----------------------------- 7) DEBUG: OTHER/UNKNOWN ------------------------

cat("=== 7) DEBUG: 'Other' / unknown genus examples ===\n")

other_debug <- df_prep %>%
  filter(Type == "Other" | is.na(genus) | genus == "") %>%
  transmute(
    site_id,
    continent,
    lat, lon,
    species_name,
    genus,
    leaf_type = leaf_type_safe,
    Type
  ) %>%
  slice_head(n = CFG$N_DEBUG_OTHER)

cat("Other/unknown rows (showing up to ", CFG$N_DEBUG_OTHER, "): ", nrow(other_debug), "\n", sep = "")
print(other_debug, n = min(nrow(other_debug), 25))

readr::write_csv(other_debug, path_out_other_debug)
cat("\nSaved: ", path_out_other_debug, "\n\n", sep = "")

# ----------------------------- DONE ------------------------------------------

cat("=== 01c DONE ===\n")
cat("Outputs:\n")
cat(" - ", path_out_geo_flags, "\n", sep = "")
cat(" - ", path_out_global_types, "\n", sep = "")
cat(" - ", path_out_cont_summary, "\n", sep = "")
cat(" - ", path_out_cont_table, "\n", sep = "")
cat(" - ", path_out_genus_table, "\n", sep = "")
cat(" - ", path_out_other_debug, "\n", sep = "")
