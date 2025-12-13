###############################################################################
# 02c_site_maps.R
#
# Description:
#   Creates spatial visualizations (maps) of site distributions.
#   Generates 4 maps (Before/After filtering).
#
#   1. Map 1: ALL SITES - Genus
#   2. Map 2: ALL SITES - N Deposition
#   3. Map 3: FILTERED SITES (Conifers Only, No "Others", <= 15 kg N) - Genus
#   4. Map 4: FILTERED SITES (Conifers Only, No "Others", <= 15 kg N) - N Dep
#
#   * Note on N-Deposition: The variable 'ndep_total' represents the sum
#     of Reduced (NHx) and Oxidized (NOy) Nitrogen.
#
# Inputs:
#   - output/01d_merge_climate_and_ndep/data/01d_Analysis_Data_Merged.rds
#   - input/itrdb_species_with_leaf_type_and_genus.csv
#
# Outputs:
#   - output/02c_site_maps/figures/*.png
#   - Console logs of data counts.
###############################################################################

# --- 1) Packages -------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  here,
  fs,
  dplyr,
  readr,
  tidyr,      # Required for drop_na
  ggplot2,
  sf,
  rnaturalearth,
  rnaturalearthdata,
  viridis,
  scales
)

# --- 2) Paths & Setup --------------------------------------------------------

step_01d <- "01d_merge_climate_and_ndep"
current  <- "02c_site_maps"

# Input Files
path_data_merged <- here::here("output", step_01d, "data", "01d_Analysis_Data_Merged.rds")
path_spp_meta    <- here::here("input", "itrdb_species_with_leaf_type_and_genus.csv")

# Output Directories
base_out_dir <- here::here("output", current)
out_dir_fig  <- fs::path(base_out_dir, "figures")
fs::dir_create(out_dir_fig)

# Graphics Theme
theme_map <- theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#F0F8FF", color = NA), # Light blue water
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    plot.caption = element_text(size = 10, face = "italic", hjust = 1) # Right-aligned caption
  )

cat("###############################################################################\n")
cat("# 02c SITE MAPS (STRICT FILTERING & CAPTIONS)\n")
cat("###############################################################################\n\n")

# --- 3) Load Data ------------------------------------------------------------

if (!file.exists(path_data_merged)) stop("Merged data not found!")
df <- readRDS(path_data_merged)

if (!file.exists(path_spp_meta)) stop("Species metadata CSV not found!")
spp_meta <- read_csv(path_spp_meta, show_col_types = FALSE)

cat("=== DATA LOADED ===\n")
cat(sprintf("Main Data:    %d rows\n", nrow(df)))

# --- 4) Prepare Site Data (Aggregated) ---------------------------------------

cat("\n--- Aggregating sites (Optimized Speed) ---\n")
# NOTE: ndep_total is the sum of NHx (Reduced) and NOy (Oxidized)

# 1. Aggregate main data to site level first
site_agg <- df %>%
  group_by(site_code, lat, lon, species_code) %>%
  summarise(
    mean_ndep_g_m2 = mean(ndep_total_5yr, na.rm = TRUE),
    n_trees = n_distinct(tree_id),
    .groups = "drop"
  )

# 2. Join metadata
site_info_all <- site_agg %>%
  left_join(spp_meta, by = "species_code") %>%
  mutate(
    mean_ndep_kg_ha = mean_ndep_g_m2 * 10 # Convert to kg/ha
  )

# 3. Define Main Groups vs Others
top_genera <- c("Abies", "Picea", "Pinus", "Larix", "Tsuga", "Pseudotsuga", "Quercus", "Fagus")

site_info_all <- site_info_all %>%
  mutate(
    genus_plot = ifelse(genus %in% top_genera, genus, "Other"),
    genus_plot = factor(genus_plot, levels = c(top_genera, "Other"))
  )

# --- 5) CONSOLE REPORT: BEFORE FILTERING -------------------------------------

cat("\n###############################################################################\n")
cat("# REPORT: ALL SITES (BEFORE FILTERING)\n")
cat("###############################################################################\n")

cat(sprintf("Total Sites: %d\n", nrow(site_info_all)))
cat("--- Counts by Genus (All) ---\n")
print(sort(table(site_info_all$genus), decreasing = TRUE))

# --- 6) MAPS 1 & 2: ALL SITES ------------------------------------------------

cat("\nGenerating maps for ALL sites...\n")

world <- ne_countries(scale = "medium", returnclass = "sf")
sites_sf_all <- st_as_sf(site_info_all %>% drop_na(lat, lon), coords = c("lon", "lat"), crs = 4326)
count_all <- nrow(sites_sf_all)

# MAP 1: By Genus
p1 <- ggplot() +
  geom_sf(data = world, fill = "grey92", color = "grey80", linewidth = 0.2) +
  geom_sf(data = sites_sf_all, aes(color = genus_plot), size = 1.2, alpha = 0.8) +
  scale_color_brewer(palette = "Paired", name = "Genus") +
  coord_sf(ylim = c(45, 85), expand = FALSE) +
  labs(
    title = "Map 1: All Sites by Genus", 
    subtitle = "Before Filtering",
    caption = sprintf("Total Sites: N = %d", count_all)
  ) +
  theme_map +
  guides(color = guide_legend(override.aes = list(size = 3), nrow = 2))

ggsave(fs::path(out_dir_fig, "02c_Map1_All_Genus.png"), p1, width = 12, height = 6)
cat("  -> Saved: 02c_Map1_All_Genus.png\n")

# MAP 2: By N-Deposition
p2 <- ggplot() +
  geom_sf(data = world, fill = "grey92", color = "grey80", linewidth = 0.2) +
  geom_sf(data = sites_sf_all, aes(color = mean_ndep_kg_ha), size = 1.2, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "N Dep (kg/ha)", end = 0.9) +
  coord_sf(ylim = c(45, 85), expand = FALSE) +
  labs(
    title = "Map 2: All Sites by N-Deposition", 
    subtitle = "Before Filtering (Total N = NHx + NOy)",
    caption = sprintf("Total Sites: N = %d", count_all)
  ) +
  theme_map

ggsave(fs::path(out_dir_fig, "02c_Map2_All_NDep.png"), p2, width = 12, height = 6)
cat("  -> Saved: 02c_Map2_All_NDep.png\n")

# --- 7) STRICT FILTERING -----------------------------------------------------

NDEP_LIMIT <- 15

cat("\n###############################################################################\n")
cat("# APPLYING STRICT FILTERS\n")
cat("###############################################################################\n")
cat("1. KEEP ONLY: leaf_type == 'conifer'\n")
cat("2. REMOVE:    genus_plot == 'Other' (e.g., Chamaecyparis, Thuja)\n")
cat(sprintf("3. KEEP ONLY: N Deposition <= %s kg/ha\n", NDEP_LIMIT))

site_info_filt <- site_info_all %>%
  filter(
    leaf_type == "conifer",           # 1. Conifers only
    genus_plot != "Other",            # 2. No "Others"
    mean_ndep_kg_ha <= NDEP_LIMIT     # 3. Low N-Dep
  ) %>%
  mutate(
    # Remove unused factor levels
    genus_plot = droplevels(genus_plot)
  )

# Calculate Stats
n_start     <- nrow(site_info_all)
n_end       <- nrow(site_info_filt)
n_removed   <- n_start - n_end

# Breakdown
n_broadleaf <- nrow(site_info_all %>% filter(leaf_type != "conifer"))
n_others    <- nrow(site_info_all %>% filter(leaf_type == "conifer", genus_plot == "Other"))
n_high_n    <- nrow(site_info_all %>% filter(leaf_type == "conifer", genus_plot != "Other", mean_ndep_kg_ha > NDEP_LIMIT))

# --- 8) CONSOLE REPORT: AFTER FILTERING --------------------------------------

cat("\n###############################################################################\n")
cat("# REPORT: FILTERED SITES (FINAL)\n")
cat("###############################################################################\n")

cat(sprintf("Original Sites:  %d\n", n_start))
cat(sprintf("Remaining Sites: %d\n", n_end))
cat(sprintf("Total Removed:   %d\n", n_removed))
cat("   - Broadleafs removed: ", n_broadleaf, "\n")
cat("   - 'Other' Conifers removed:", n_others, "\n")
cat("   - High N-Dep (>15) removed:", n_high_n, "\n\n")

cat("--- Remaining Genera ---\n")
print(sort(table(site_info_filt$genus), decreasing = TRUE))

cat("\n--- Remaining N Deposition Stats ---\n")
print(summary(site_info_filt$mean_ndep_kg_ha))

# --- 9) MAPS 3 & 4: FILTERED SITES -------------------------------------------

cat("\nGenerating maps for FILTERED sites...\n")

sites_sf_filt <- st_as_sf(site_info_filt %>% drop_na(lat, lon), coords = c("lon", "lat"), crs = 4326)
count_filt <- nrow(sites_sf_filt)

# MAP 3: Filtered - By Genus
p3 <- ggplot() +
  geom_sf(data = world, fill = "grey92", color = "grey80", linewidth = 0.2) +
  geom_sf(data = sites_sf_filt, aes(color = genus_plot), size = 1.2, alpha = 0.8) +
  scale_color_brewer(palette = "Paired", name = "Genus") +
  coord_sf(ylim = c(45, 85), expand = FALSE) +
  labs(
    title = "Map 3: Filtered Sites by Genus", 
    subtitle = "Major Conifers Only, Low N-Deposition",
    caption = sprintf("Included Sites: N = %d", count_filt)
  ) +
  theme_map +
  guides(color = guide_legend(override.aes = list(size = 3), nrow = 1))

ggsave(fs::path(out_dir_fig, "02c_Map3_Filtered_Genus.png"), p3, width = 12, height = 6)
cat("  -> Saved: 02c_Map3_Filtered_Genus.png\n")

# MAP 4: Filtered - By N-Deposition
p4 <- ggplot() +
  geom_sf(data = world, fill = "grey92", color = "grey80", linewidth = 0.2) +
  geom_sf(data = sites_sf_filt, aes(color = mean_ndep_kg_ha), size = 1.2, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "N Dep (kg/ha)", end = 0.9) +
  coord_sf(ylim = c(45, 85), expand = FALSE) +
  labs(
    title = "Map 4: Filtered Sites by N-Deposition", 
    subtitle = "Major Conifers Only, Low N-Deposition (Total N = NHx + NOy)",
    caption = sprintf("Included Sites: N = %d", count_filt)
  ) +
  theme_map

ggsave(fs::path(out_dir_fig, "02c_Map4_Filtered_NDep.png"), p4, width = 12, height = 6)
cat("  -> Saved: 02c_Map4_Filtered_NDep.png\n")

cat("\n###############################################################################\n")
cat("# 02c COMPLETE\n")
cat("###############################################################################\n")