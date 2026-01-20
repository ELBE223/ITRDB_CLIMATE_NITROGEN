###############################################################################
# 05a_collect_outputs.R
#
# Step 1: Load CSV outputs from 04c / 04d / 04e and print ALL created values.
# Writes one combined long CSV for quick inspection.
#
# Fix: tidyr::pivot_longer() uses vctrs and does NOT implicitly cast
# numeric <-> character. We therefore transform all pivoted values to character.
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, fs, readr, dplyr, tidyr, stringr)

options(width = 180)

CFG <- list(
  H1_DIR  = here::here("output", "04_c_LMM_M0_M1", "data"),
  H2_DIR  = here::here("output", "04_d_LMM_H2_JJA_DJF", "data"),
  E1_DIR  = here::here("output", "04_e_LMM_E1_genus_slopes", "data"),
  OUT_DIR = here::here("output", "05_collect_results", "data"),
  
  H1_FILES = c(
    "04c_M0_M1_model_compare.csv",
    "04c_M0_M1_fixed_effects.csv"
  ),
  H2_FILES = c(
    "04d_H2_JJA_DJF_model_compare_all.csv",
    "04d_H2_JJA_DJF_fixed_effects_all.csv",
    "04d_H2_JJA_model_compare.csv",
    "04d_H2_JJA_fixed_effects.csv",
    "04d_H2_DJF_model_compare.csv",
    "04d_H2_DJF_fixed_effects.csv"
  ),
  E1_FILES = c(
    "04e_E1_model_compare.csv",
    "04e_E1_fixed_effects.csv",
    "04e_E1_genus_slopes.csv"
  ),
  
  OUT_LONG = "05a_all_outputs_long.csv"
)

fs::dir_create(CFG$OUT_DIR)

read_csv_checked <- function(path) {
  if (!fs::file_exists(path)) {
    stop("Missing file: ", path, "\nRun the corresponding 04* script first.")
  }
  readr::read_csv(path, show_col_types = FALSE)
}

print_df_info <- function(df, title, n_head = 10, n_max_full = 200) {
  cat("\n--- ", title, " ---\n", sep = "")
  cat("Rows:", format(nrow(df), big.mark = ","), " | Cols:", ncol(df), "\n")
  cat("Columns:", paste(names(df), collapse = ", "), "\n")
  
  if (nrow(df) == 0) return(invisible(NULL))
  
  if (nrow(df) <= n_max_full) {
    print(df)
  } else {
    cat("\n(Showing head and tail because table is large)\n")
    print(utils::head(df, n_head))
    cat("\n...")
    print(utils::tail(df, n_head))
  }
}

cat("=== 05a: Collecting CSV outputs ===\n")

# ------------------------------ H1 (04c) -------------------------------------
h1_paths <- fs::path(CFG$H1_DIR, CFG$H1_FILES)
H1 <- lapply(h1_paths, read_csv_checked)
names(H1) <- CFG$H1_FILES

print_df_info(H1[["04c_M0_M1_model_compare.csv"]], "H1: 04c_M0_M1_model_compare")
print_df_info(H1[["04c_M0_M1_fixed_effects.csv"]], "H1: 04c_M0_M1_fixed_effects")

h1_key <- H1[["04c_M0_M1_fixed_effects.csv"]] %>%
  filter(model == "M1", term == "N_depo1_z")
cat("\nH1 key term (M1: N_depo1_z):\n")
print(h1_key)

# ------------------------------ H2 (04d) -------------------------------------
h2_paths <- fs::path(CFG$H2_DIR, CFG$H2_FILES)
H2 <- lapply(h2_paths, read_csv_checked)
names(H2) <- CFG$H2_FILES

print_df_info(H2[["04d_H2_JJA_DJF_model_compare_all.csv"]], "H2: model_compare_all")
print_df_info(H2[["04d_H2_JJA_DJF_fixed_effects_all.csv"]], "H2: fixed_effects_all")

h2_fe_all <- H2[["04d_H2_JJA_DJF_fixed_effects_all.csv"]]

h2_int <- h2_fe_all %>%
  filter(model == "M2") %>%
  filter(str_detect(term, "N_depo1_z") & str_detect(term, ":"))

cat("\nH2 interaction terms (M2, all seasons; pattern: contains 'N_depo1_z' and ':'):\n")
if (nrow(h2_int) == 0) {
  cat("(No rows matched. Term naming may differ. Printing all M2 terms containing 'N_depo1_z' instead.)\n")
  print(h2_fe_all %>% filter(model == "M2") %>% filter(str_detect(term, "N_depo1_z")))
} else {
  print(h2_int)
}

# ------------------------------ E1 (04e) -------------------------------------
e1_paths <- fs::path(CFG$E1_DIR, CFG$E1_FILES)
E1 <- lapply(e1_paths, read_csv_checked)
names(E1) <- CFG$E1_FILES

print_df_info(E1[["04e_E1_model_compare.csv"]], "E1: 04e_E1_model_compare")
print_df_info(E1[["04e_E1_fixed_effects.csv"]], "E1: 04e_E1_fixed_effects")
print_df_info(E1[["04e_E1_genus_slopes.csv"]], "E1: 04e_E1_genus_slopes")

cat("\nE1 slope range (min/max estimate):\n")
if ("estimate" %in% names(E1[["04e_E1_genus_slopes.csv"]])) {
  print(range(E1[["04e_E1_genus_slopes.csv"]]$estimate, na.rm = TRUE))
}

# ------------------------------ combined long export --------------------------
to_long <- function(df, source_file) {
  df %>%
    mutate(.source_file = source_file, .row_id = row_number()) %>%
    pivot_longer(
      cols = -c(.source_file, .row_id),
      names_to = "field",
      values_to = "value",
      values_transform = list(value = as.character)
    )
}

all_long <- bind_rows(
  bind_rows(lapply(names(H1), function(nm) to_long(H1[[nm]], nm))),
  bind_rows(lapply(names(H2), function(nm) to_long(H2[[nm]], nm))),
  bind_rows(lapply(names(E1), function(nm) to_long(E1[[nm]], nm)))
)

out_path <- fs::path(CFG$OUT_DIR, CFG$OUT_LONG)
readr::write_csv(all_long, out_path)

cat("\nWrote combined long CSV to:\n", out_path, "\n", sep = "")
cat("\nDone.\n")
