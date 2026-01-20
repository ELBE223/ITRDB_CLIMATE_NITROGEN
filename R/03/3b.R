###############################################################################
# 03b_deskriptiv_maps_box_bar.R
#
# Uses final sites from 03a (output/03_a/data/03a_site_year.parquet).
# Restricts ALL plots to those final sites (and to Asia/Europe/NorthAmerica).
# N deposition map + N boxplot are capped to <= 15 kg N/ha/yr (legend + plotting).
# Prints:
# - Final site counts + percentages
# - Top 6 genera (trees + sites + %)
# - N totals + shares by dominant genus (site-year sums + site-mean sums)
###############################################################################

# ------------------------------ Packages --------------------------------------
pkgs <- c(
  "here", "dplyr", "fs", "arrow", "tibble", "rlang",
  "ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "scales"
)
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

options(scipen = 999)

# ------------------------------ Config ----------------------------------------
CFG <- list(
  YEAR_MIN = 1901L,
  YEAR_MAX = 2010L,
  
  N_CAP = 15,
  Y_CAP_Q = 0.995,
  
  PATH_FINAL = here::here("output", "02_env_merge", "data", "02_panel_with_env.parquet"),
  PATH_03A_SITE_YEAR = here::here("output", "03_a", "data", "03a_site_year.parquet"),
  
  OUT_DIR   = here::here("output", "03_b"),
  OUT_PLOTS = here::here("output", "03_b", "plots"),
  OUT_DATA  = here::here("output", "03_b", "data")
)

fs::dir_create(CFG$OUT_DIR)
fs::dir_create(CFG$OUT_PLOTS)
fs::dir_create(CFG$OUT_DATA)

stopifnot(fs::file_exists(CFG$PATH_FINAL))
stopifnot(fs::file_exists(CFG$PATH_03A_SITE_YEAR))

# ------------------------------ Helpers ---------------------------------------
continent_de <- function(x) {
  dplyr::case_when(
    x %in% c("Asia") ~ "Asien",
    x %in% c("Europe") ~ "Europa",
    x %in% c("NorthAmerica", "North America", "N. America") ~ "Nordamerika",
    TRUE ~ as.character(x)
  )
}

save_png <- function(p, filename, w = 10, h = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = fs::path(CFG$OUT_PLOTS, filename),
    plot = p, width = w, height = h, units = "in", dpi = dpi
  )
}

theme_map <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      plot.title   = element_text(face = "bold")
    )
}

theme_orig_style <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "plain")
    )
}

fmt_num <- function(x, digits = 2) {
  if (!is.finite(x) || is.na(x)) return(NA_character_)
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

# ------------------------------ Load ------------------------------------------
ds   <- arrow::open_dataset(CFG$PATH_FINAL, format = "parquet")
cols <- ds$schema$names
stopifnot(all(c("site_id", "year", "continent", "lat", "lon", "genus") %in% cols))

tree_uid_col <- if ("tree_uid" %in% cols) "tree_uid" else if ("tree_id" %in% cols) "tree_id" else NA_character_
stopifnot(!is.na(tree_uid_col))

# ------------------------------ Final sites (03a) ------------------------------
site_year_03a <- arrow::read_parquet(CFG$PATH_03A_SITE_YEAR)
stopifnot("site_id" %in% names(site_year_03a))
final_sites <- unique(site_year_03a$site_id)

# ------------------------------ Base year filter (for % reporting) -------------
ds_year <- ds %>%
  filter(year >= CFG$YEAR_MIN, year <= CFG$YEAR_MAX)

site_meta_year_all <- ds_year %>%
  select(site_id, continent) %>%
  distinct() %>%
  collect() %>%
  mutate(
    continent = as.character(continent),
    continent_de = continent_de(continent)
  )

site_meta_year_3cont <- site_meta_year_all %>%
  filter(continent_de %in% c("Asien", "Europa", "Nordamerika"))

n_sites_year_all   <- nrow(site_meta_year_all)
n_sites_year_3cont <- nrow(site_meta_year_3cont)

n_sites_final_all <- length(final_sites)

# ------------------------------ Apply final site filter ------------------------
ds_filt <- ds_year %>%
  filter(site_id %in% final_sites)

# ------------------------------ Site meta (3 continents) -----------------------
site_meta <- ds_filt %>%
  select(site_id, continent, lat, lon) %>%
  distinct() %>%
  collect() %>%
  mutate(
    continent = as.character(continent),
    continent_de = continent_de(continent)
  ) %>%
  filter(continent_de %in% c("Asien", "Europa", "Nordamerika"))

final_sites_3cont <- unique(site_meta$site_id)
n_sites_final_3cont <- length(final_sites_3cont)

pct_final_all_vs_year_all     <- 100 * n_sites_final_all / n_sites_year_all
pct_final_3cont_vs_year_3cont <- 100 * n_sites_final_3cont / n_sites_year_3cont

# Restrict dataset further to the plotted continents (keeps everything consistent)
ds_filt <- ds_filt %>% filter(site_id %in% final_sites_3cont)

# ------------------------------ N deposition selection -------------------------
use_n_depo1 <- "N_depo1" %in% cols
use_n_sum   <- all(c("NHx_1", "NOy_1") %in% cols)
stopifnot(use_n_depo1 | use_n_sum)

n_expr <- if (use_n_depo1) {
  rlang::sym("N_depo1")
} else {
  rlang::expr(NHx_1 + NOy_1)
}
n_label_src <- if (use_n_depo1) "N_depo1" else "NHx_1 + NOy_1"

# ------------------------------ Tree meta -------------------------------------
tree_meta <- ds_filt %>%
  select(!!rlang::sym(tree_uid_col), site_id, continent, genus) %>%
  distinct() %>%
  collect() %>%
  mutate(
    continent = as.character(continent),
    genus = as.character(genus),
    continent_de = continent_de(continent)
  ) %>%
  filter(continent_de %in% c("Asien", "Europa", "Nordamerika"))

# Top 6 genera by tree count (within FINAL sites)
top6_tbl <- tree_meta %>%
  count(genus, sort = TRUE, name = "n_trees") %>%
  slice_head(n = 6)

top6_genera <- top6_tbl$genus

classic_order <- c("Abies", "Larix", "Picea", "Pinus", "Pseudotsuga", "Tsuga")
if (all(classic_order %in% top6_genera)) top6_genera <- classic_order

# Dominant genus per site (by tree counts)
dom_genus_site <- tree_meta %>%
  count(site_id, genus, name = "n_trees") %>%
  group_by(site_id) %>%
  slice_max(order_by = n_trees, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(genus6 = if_else(genus %in% top6_genera, genus, "Andere")) %>%
  left_join(site_meta %>% select(site_id, continent_de, lat, lon), by = "site_id") %>%
  mutate(
    continent_de = factor(continent_de, levels = c("Asien", "Europa", "Nordamerika")),
    genus6 = factor(genus6, levels = c(top6_genera, "Andere"))
  )

# ------------------------------ Site-year N + site mean N ----------------------
site_year_n <- ds_filt %>%
  mutate(N_depo_total = !!n_expr) %>%
  select(site_id, year, N_depo_total) %>%
  group_by(site_id, year) %>%
  summarise(N_depo_total = mean(N_depo_total, na.rm = TRUE), .groups = "drop") %>%
  collect()

site_n_mean <- site_year_n %>%
  group_by(site_id) %>%
  summarise(N_depo_mean = mean(N_depo_total, na.rm = TRUE), .groups = "drop")

site_meta_n <- site_meta %>%
  left_join(site_n_mean, by = "site_id") %>%
  mutate(continent_de = factor(continent_de, levels = c("Asien", "Europa", "Nordamerika")))

# ------------------------------ Tree mean BAI (cm²/yr) -------------------------
EPS_BAI <- 1e-4

bai_cm2_col <- dplyr::case_when(
  "BAI_cm2" %in% cols ~ "BAI_cm2",
  "bai_cm2" %in% cols ~ "bai_cm2",
  "BAI"     %in% cols ~ "BAI",
  TRUE ~ NA_character_
)

if (!is.na(bai_cm2_col)) {
  tree_bai <- ds_filt %>%
    select(!!rlang::sym(tree_uid_col), site_id, continent, genus, !!rlang::sym(bai_cm2_col)) %>%
    group_by(!!rlang::sym(tree_uid_col), site_id, continent, genus) %>%
    summarise(BAI_cm2 = mean(!!rlang::sym(bai_cm2_col), na.rm = TRUE), .groups = "drop") %>%
    collect()
} else {
  stopifnot("log_BAI_eps" %in% cols)
  tree_bai <- ds_filt %>%
    select(!!rlang::sym(tree_uid_col), site_id, continent, genus, log_BAI_eps) %>%
    group_by(!!rlang::sym(tree_uid_col), site_id, continent, genus) %>%
    summarise(log_BAI_eps_mean = mean(log_BAI_eps, na.rm = TRUE), .groups = "drop") %>%
    collect() %>%
    mutate(BAI_cm2 = pmax(0, exp(log_BAI_eps_mean) - EPS_BAI))
}

tree_bai <- tree_bai %>%
  mutate(
    continent = as.character(continent),
    genus = as.character(genus),
    continent_de = continent_de(continent),
    genus6 = if_else(genus %in% top6_genera, genus, NA_character_)
  ) %>%
  filter(continent_de %in% c("Asien", "Europa", "Nordamerika")) %>%
  filter(!is.na(genus6)) %>%
  left_join(site_n_mean, by = "site_id") %>%
  filter(is.finite(BAI_cm2), !is.na(BAI_cm2)) %>%
  mutate(
    continent_de = factor(continent_de, levels = c("Asien", "Europa", "Nordamerika")),
    genus6 = factor(genus6, levels = top6_genera)
  )

ylim_N   <- CFG$N_CAP
ylim_BAI <- as.numeric(stats::quantile(tree_bai$BAI_cm2, probs = CFG$Y_CAP_Q, na.rm = TRUE))
if (!is.finite(ylim_BAI)) ylim_BAI <- max(tree_bai$BAI_cm2, na.rm = TRUE)

# ------------------------------ Prints: sites + % ------------------------------
message("---- 03b SUMMARY (final sites from 03a) ----")
message("Sites after year filter (all continents): ", format(n_sites_year_all, big.mark = ","))
message("Sites final from 03a (all continents):    ", format(n_sites_final_all, big.mark = ","),
        " (", fmt_num(pct_final_all_vs_year_all, 1), "%)")
message("Sites after year filter (3 continents):   ", format(n_sites_year_3cont, big.mark = ","))
message("Sites final (3 continents plotted):       ", format(n_sites_final_3cont, big.mark = ","),
        " (", fmt_num(pct_final_3cont_vs_year_3cont, 1), "%)")
message("N source: ", n_label_src, " | Plot cap: ", CFG$N_CAP, " kg N/ha/yr")

# ------------------------------ Prints: Top-6 with % (trees + sites) ----------
tot_trees_all <- dplyr::n_distinct(tree_meta[[tree_uid_col]])
tot_sites_all <- dplyr::n_distinct(tree_meta$site_id)

top6_print <- tree_meta %>%
  count(genus, sort = TRUE, name = "n_trees") %>%
  slice_head(n = 6) %>%
  left_join(
    tree_meta %>% distinct(site_id, genus) %>% count(genus, name = "n_sites"),
    by = "genus"
  ) %>%
  mutate(
    pct_trees = 100 * n_trees / tot_trees_all,
    pct_sites = 100 * n_sites / tot_sites_all
  )

message("Top 6 genera (within final sites; incl. % of trees and % of sites):")
print(top6_print)

# ------------------------------ Prints: N totals + shares by dominant genus ----
site_dom <- dom_genus_site %>% select(site_id, dom_genus = genus6)

# N summary overall
n_site_mean <- site_n_mean$N_depo_mean
n_siteyear  <- site_year_n$N_depo_total

message("N (site mean)  mean=", fmt_num(mean(n_site_mean, na.rm = TRUE), 2),
        " median=", fmt_num(stats::median(n_site_mean, na.rm = TRUE), 2),
        " min=", fmt_num(min(n_site_mean, na.rm = TRUE), 2),
        " max=", fmt_num(max(n_site_mean, na.rm = TRUE), 2),
        " sum(site means)=", fmt_num(sum(n_site_mean, na.rm = TRUE), 2))

message("N (site-year)  mean=", fmt_num(mean(n_siteyear, na.rm = TRUE), 2),
        " sum(site-years)=", fmt_num(sum(n_siteyear, na.rm = TRUE), 2))

# Share of N_total across dominant genus (site-year sums)
n_by_dom_siteyear <- site_year_n %>%
  left_join(site_dom, by = "site_id") %>%
  filter(!is.na(dom_genus), dom_genus != "Andere") %>%
  group_by(dom_genus) %>%
  summarise(
    N_sum_siteyear = sum(N_depo_total, na.rm = TRUE),
    N_mean_siteyear = mean(N_depo_total, na.rm = TRUE),
    n_siteyears = n(),
    n_sites = n_distinct(site_id),
    .groups = "drop"
  ) %>%
  mutate(
    pct_N_siteyear = 100 * N_sum_siteyear / sum(N_sum_siteyear, na.rm = TRUE)
  ) %>%
  arrange(desc(N_sum_siteyear))

message("N shares by dominant genus (site-year sums):")
print(n_by_dom_siteyear)

# Share of N across dominant genus (sum of site means)
n_by_dom_sitemean <- site_n_mean %>%
  left_join(site_dom, by = "site_id") %>%
  filter(!is.na(dom_genus), dom_genus != "Andere") %>%
  group_by(dom_genus) %>%
  summarise(
    N_sum_sitemean = sum(N_depo_mean, na.rm = TRUE),
    N_mean_sitemean = mean(N_depo_mean, na.rm = TRUE),
    n_sites = n(),
    .groups = "drop"
  ) %>%
  mutate(
    pct_N_sitemean = 100 * N_sum_sitemean / sum(N_sum_sitemean, na.rm = TRUE)
  ) %>%
  arrange(desc(N_sum_sitemean))

message("N shares by dominant genus (sum of site means):")
print(n_by_dom_sitemean)

# ------------------------------ Optional exports ------------------------------
arrow::write_parquet(site_meta_n, fs::path(CFG$OUT_DATA, "03b_site_meta_n.parquet"), compression = "zstd")
arrow::write_parquet(dom_genus_site, fs::path(CFG$OUT_DATA, "03b_dom_genus_site.parquet"), compression = "zstd")
arrow::write_parquet(tree_bai, fs::path(CFG$OUT_DATA, "03b_tree_bai.parquet"), compression = "zstd")

# ------------------------------ World basemap ---------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(4326)

# ------------------------------ MAP 01: continents ----------------------------
map01 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = site_meta %>% filter(is.finite(lat), is.finite(lon)),
    aes(x = lon, y = lat, color = continent_de),
    size = 1.4, alpha = 0.85
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(30, 90), expand = FALSE) +
  labs(
    title = "",
    color = "Kontinent"
  ) +
  theme_map()

save_png(map01, "03b_map_01_continents.png", w = 11, h = 6.5)

# ------------------------------ MAP 02: genus ---------------------------------
map02 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = dom_genus_site %>% filter(is.finite(lat), is.finite(lon), genus6 != "Andere"),
    aes(x = lon, y = lat, color = genus6),
    size = 1.4, alpha = 0.85
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(30, 90), expand = FALSE) +
  labs(
    title = "",
    color = "Baumgattung"
  ) +
  theme_map()

save_png(map02, "03b_map_02_genus.png", w = 11, h = 6.5)

# ------------------------------ MAP 03: N deposition (capped) ------------------
map03 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = site_meta_n %>% filter(is.finite(lat), is.finite(lon), is.finite(N_depo_mean)),
    aes(x = lon, y = lat, color = N_depo_mean),
    size = 1.6, alpha = 0.9
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(30, 90), expand = FALSE) +
  scale_color_viridis_c(
    option = "C",
    limits = c(0, CFG$N_CAP),
    oob = scales::squish,
    labels = scales::label_number(accuracy = 0.1)
  ) +
  labs(
    title = paste0(""),
    color = "kg N/ha/Jahr"
  ) +
  theme_map()

save_png(map03, "03b_map_03_ndep.png", w = 11, h = 6.5)

# ------------------------------ BOXPLOT N (fixed y to cap) ---------------------
boxN <- tree_bai %>%
  filter(is.finite(N_depo_mean), !is.na(N_depo_mean)) %>%
  ggplot(aes(x = genus6, y = N_depo_mean)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~continent_de, nrow = 1) +
  coord_cartesian(ylim = c(0, ylim_N)) +
  labs(
    title = "",
    x = "Baumgattung",
    y = paste0("Stickstoffdeposition (kg N/ha/Jahr)")
  ) +
  theme_orig_style()

save_png(boxN, "03b_box_N_orig_style.png", w = 12, h = 4.8)

# ------------------------------ BOXPLOT BAI (fixed y) --------------------------
boxBAI <- tree_bai %>%
  ggplot(aes(x = genus6, y = BAI_cm2)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~continent_de, nrow = 1) +
  coord_cartesian(ylim = c(0, ylim_BAI)) +
  labs(
    title = "",
    x = "Baumgattung",
    y = "BAI (cm²/Jahr)"
  ) +
  theme_orig_style()

save_png(boxBAI, "03b_box_BAI_orig_style.png", w = 12, h = 4.8)

# ------------------------------ BAR: sites per dominant genus ------------------
bar_sites <- dom_genus_site %>%
  filter(genus6 != "Andere") %>%
  count(genus6, name = "n_sites") %>%
  mutate(genus6 = factor(genus6, levels = top6_genera)) %>%
  ggplot(aes(x = genus6, y = n_sites)) +
  geom_col() +
  labs(
    title = "",
    x = "Baumgattung",
    y = "Anzahl der Standorte"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_png(bar_sites, "03b_bar_sites.png", w = 9, h = 5)

message("03b fertig. Output: ", CFG$OUT_DIR)
message("Plots:  ", CFG$OUT_PLOTS)
message("Data:   ", CFG$OUT_DATA)
