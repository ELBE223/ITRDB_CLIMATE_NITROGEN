###############################################################################
# 02a_extract_env_site_monthly.R
#
# Fast monthly extraction at sites:
# - CRU TS: tmp, pre
# - ISIMIP3a: nhx, noy  -> converted to kg N ha-1 mon-1
#
# Output (parquet):
#   output/02_env_extract/data/
#     02a_site_monthly_climate.parquet
#     02a_site_monthly_ndep.parquet
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, dplyr, stringr, fs, tibble,
  arrow, terra, R.utils, ncdf4, progress
)

terra::terraOptions(progress = 0)

CFG <- list(
  PATH_PANEL_01D = here::here("output", "01d_prepare_LMM_panel",
                              "01d_LMM_panel_conifers_lat50_genus15.parquet"),
  DIR_00B = here::here("output", "00b_download_ISIMIP_CRU", "data", "climate"),
  
  YEAR_MIN = 1901L,
  YEAR_MAX = 2010L,
  
  CHUNK_MONTHS = 120L,
  
  OUT_DIR  = here::here("output", "02_env_extract", "data"),
  OUT_CLIM = "02a_site_monthly_climate.parquet",
  OUT_NDEP = "02a_site_monthly_ndep.parquet"
)

fs::dir_create(CFG$OUT_DIR)
out_clim <- fs::path(CFG$OUT_DIR, CFG$OUT_CLIM)
out_ndep <- fs::path(CFG$OUT_DIR, CFG$OUT_NDEP)

cat("=== 02a Extract monthly env at sites ===\n")
cat("01d panel: ", CFG$PATH_PANEL_01D, "\n", sep = "")
cat("00b dir:   ", CFG$DIR_00B, "\n", sep = "")
cat("Years:     ", CFG$YEAR_MIN, "-", CFG$YEAR_MAX, "\n", sep = "")
cat("Chunk:     ", CFG$CHUNK_MONTHS, " months\n\n", sep = "")

stopifnot(fs::file_exists(CFG$PATH_PANEL_01D))
stopifnot(fs::dir_exists(CFG$DIR_00B))

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ----------------------------- sites -----------------------------------------

panel_ds <- arrow::open_dataset(CFG$PATH_PANEL_01D, format = "parquet")
sites <- panel_ds %>%
  select(site_id, lat, lon) %>%
  distinct() %>%
  collect()

cat("Sites found: ", nrow(sites), "\n\n", sep = "")

pts <- terra::vect(
  sites %>% transmute(x = lon, y = lat),
  geom = c("x", "y"),
  crs = "EPSG:4326"
)
pts$site_id <- sites$site_id

# ----------------------------- helpers ---------------------------------------

find_one <- function(dir, pattern) {
  f <- fs::dir_ls(dir, recurse = TRUE, type = "file", regexp = pattern)
  if (length(f) == 0) stop("File not found for pattern: ", pattern, " in ", dir)
  f <- f[order(stringr::str_detect(f, "\\.gz$"), f)]  # prefer .nc over .nc.gz
  f[1]
}

ensure_unzipped <- function(path) {
  if (!stringr::str_detect(path, "\\.gz$")) return(path)
  out <- stringr::str_remove(path, "\\.gz$")
  if (fs::file_exists(out)) return(out)
  cat("Unzipping: ", fs::path_file(path), "\n", sep = "")
  R.utils::gunzip(path, destname = out, remove = FALSE, overwrite = FALSE)
  out
}

rast_var <- function(nc_path, varname) {
  r <- try(terra::rast(paste0("NETCDF:", nc_path, ":", varname)), silent = TRUE)
  if (!inherits(r, "try-error") && terra::nlyr(r) > 0) return(r)
  
  s <- terra::sds(nc_path)
  nm <- names(s)
  i <- which(tolower(nm) == tolower(varname))
  if (length(i) == 0) i <- which(stringr::str_detect(tolower(nm), paste0(":", tolower(varname), "$")))
  if (length(i) == 0) i <- which(stringr::str_detect(tolower(nm), paste0("\\b", tolower(varname), "\\b")))
  if (length(i) == 0) stop("Variable '", varname, "' not found in ", nc_path, "\nFound: ", paste(nm, collapse = ", "))
  
  terra::rast(s[[i[1]]])
}

align_lon_to_raster <- function(pts, r) {
  e <- terra::ext(r)
  xmin <- e[1]; xmax <- e[2]
  
  if (xmin >= 0 && xmax > 180) {
    xy <- terra::crds(pts)
    xy[, 1] <- ifelse(xy[, 1] < 0, xy[, 1] + 360, xy[, 1])
    p2 <- terra::vect(xy, type = "points", crs = terra::crs(pts))
    p2$site_id <- pts$site_id
    return(p2)
  }
  
  if (xmin < 0 && xmax <= 180) {
    xy <- terra::crds(pts)
    xy[, 1] <- ifelse(xy[, 1] > 180, xy[, 1] - 360, xy[, 1])
    p2 <- terra::vect(xy, type = "points", crs = terra::crs(pts))
    p2$site_id <- pts$site_id
    return(p2)
  }
  
  pts
}

# robust year-month from terra time; fallback: monthly sequence from start_date
get_time_ym <- function(r, start_date = NULL) {
  tt <- terra::time(r)
  
  if (length(tt) == terra::nlyr(r) && !all(is.na(tt)) &&
      inherits(tt, c("POSIXct", "POSIXt", "Date"))) {
    return(tibble::tibble(
      idx = seq_len(terra::nlyr(r)),
      year = as.integer(format(tt, "%Y")),
      month = as.integer(format(tt, "%m"))
    ))
  }
  
  if (!is.null(start_date)) {
    d0 <- as.Date(start_date)
    if (is.na(d0)) stop("start_date not convertible to Date: ", start_date)
    
    dd <- seq.Date(d0, by = "month", length.out = terra::nlyr(r))
    return(tibble::tibble(
      idx = seq_len(terra::nlyr(r)),
      year = as.integer(format(dd, "%Y")),
      month = as.integer(format(dd, "%m"))
    ))
  }
  
  stop("No valid time() in raster and no start_date fallback provided.")
}

chunk_indices <- function(idxs, chunk_size) {
  split(idxs, ceiling(seq_along(idxs) / chunk_size))
}

extract_two_vars_monthly <- function(r1, r2, pts, ym, v1, v2, pb_label) {
  stopifnot(terra::nlyr(r1) == terra::nlyr(r2))
  
  pts2 <- align_lon_to_raster(pts, r1)
  site_ids <- pts2$site_id
  n_sites <- length(site_ids)
  
  keep <- ym[ym$year >= CFG$YEAR_MIN & ym$year <= CFG$YEAR_MAX, , drop = FALSE]
  keep <- keep[order(keep$idx), , drop = FALSE]
  
  idxs <- keep$idx
  chunks <- chunk_indices(idxs, CFG$CHUNK_MONTHS)
  
  pb <- progress::progress_bar$new(
    format = paste0(pb_label, " [:bar] :percent | :current/:total | ETA :eta"),
    total = length(chunks), clear = FALSE, width = 70
  )
  
  out_list <- vector("list", length(chunks))
  
  keep_idx   <- keep$idx
  keep_year  <- keep$year
  keep_month <- keep$month
  
  for (i in seq_along(chunks)) {
    pb$tick()
    
    idx_chunk <- chunks[[i]]
    k <- length(idx_chunk)
    
    pos <- match(idx_chunk, keep_idx)
    yr <- keep_year[pos]
    mo <- keep_month[pos]
    
    ex1 <- terra::extract(r1[[idx_chunk]], pts2, ID = FALSE)
    ex2 <- terra::extract(r2[[idx_chunk]], pts2, ID = FALSE)
    
    m1 <- as.matrix(ex1)
    m2 <- as.matrix(ex2)
    
    out_list[[i]] <- tibble::tibble(
      site_id = rep(site_ids, each = k),
      year    = rep(yr, times = n_sites),
      month   = rep(mo, times = n_sites),
      !!v1 := as.vector(t(m1)),
      !!v2 := as.vector(t(m2))
    )
  }
  
  dplyr::bind_rows(out_list) %>% dplyr::arrange(site_id, year, month)
}

read_units <- function(nc_path, var) {
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  as.character(ncdf4::ncatt_get(nc, var, "units")$value)
}

ndep_factor_to_kg_ha <- function(units_in) {
  u <- tolower(trimws(units_in %||% ""))
  u <- stringr::str_replace_all(u, "\\s+", " ")
  if (u == "" || is.na(u)) stop("Ndep units missing.")
  
  if (stringr::str_detect(u, "s-1") || stringr::str_detect(u, "/s") || stringr::str_detect(u, "sec")) {
    stop("Flux per second units not supported: '", units_in, "'.")
  }
  
  # 1 g m-2 = 10 kg ha-1
  if (stringr::str_detect(u, "mg") && stringr::str_detect(u, "m-2")) return(list(f = 0.01,   out = "kgN ha-1 mon-1"))
  if (stringr::str_detect(u, "g")  && stringr::str_detect(u, "m-2")) return(list(f = 10,     out = "kgN ha-1 mon-1"))
  if (stringr::str_detect(u, "kg") && stringr::str_detect(u, "m-2")) return(list(f = 10000,  out = "kgN ha-1 mon-1"))
  
  if (stringr::str_detect(u, "mg") && stringr::str_detect(u, "ha-1")) return(list(f = 1e-6,  out = "kgN ha-1 mon-1"))
  if (stringr::str_detect(u, "g")  && stringr::str_detect(u, "ha-1")) return(list(f = 1e-3,  out = "kgN ha-1 mon-1"))
  if (stringr::str_detect(u, "kg") && stringr::str_detect(u, "ha-1")) return(list(f = 1,     out = "kgN ha-1 mon-1"))
  
  stop("Unsupported Ndep units: '", units_in, "'")
}

# ----------------------------- CRU -------------------------------------------

cru_root <- fs::path(CFG$DIR_00B, "CRU_TS_v4.09")
tmp_path <- ensure_unzipped(find_one(cru_root, "cru_ts4\\.09\\..*\\.tmp\\.dat\\.nc(\\.gz)?$"))
pre_path <- ensure_unzipped(find_one(cru_root, "cru_ts4\\.09\\..*\\.pre\\.dat\\.nc(\\.gz)?$"))

cat("CRU tmp: ", tmp_path, "\n", sep = "")
cat("CRU pre: ", pre_path, "\n\n", sep = "")

r_tmp <- rast_var(tmp_path, "tmp")
r_pre <- rast_var(pre_path, "pre")
stopifnot(terra::nlyr(r_tmp) == terra::nlyr(r_pre))

ym_cru <- get_time_ym(r_tmp)
clim <- extract_two_vars_monthly(r_tmp, r_pre, pts, ym_cru, "tmp", "pre", "CRU")

exp_rows <- nrow(sites) * (CFG$YEAR_MAX - CFG$YEAR_MIN + 1L) * 12L
cat("\nClimate rows: ", nrow(clim), " (expected ~", exp_rows, ")\n", sep = "")
cat("Unique site-year-month: ", nrow(dplyr::distinct(clim, site_id, year, month)), "\n\n", sep = "")

arrow::write_parquet(clim, out_clim, compression = "zstd")
cat("Saved: ", out_clim, "\n\n", sep = "")

rm(clim); gc()

# ----------------------------- ISIMIP Ndep -----------------------------------

isimip_root <- fs::path(CFG$DIR_00B, "ISIMIP3a", "ndep", "histsoc")

nhx_1850 <- find_one(isimip_root, "ndep-nhx_histsoc_monthly_1850_1900\\.nc$")
nhx_1901 <- find_one(isimip_root, "ndep-nhx_histsoc_monthly_1901_2021\\.nc$")
noy_1850 <- find_one(isimip_root, "ndep-noy_histsoc_monthly_1850_1900\\.nc$")
noy_1901 <- find_one(isimip_root, "ndep-noy_histsoc_monthly_1901_2021\\.nc$")

u_nhx <- read_units(nhx_1901, "nhx")
u_noy <- read_units(noy_1901, "noy")

cat("Units nhx: '", u_nhx, "'\n", sep = "")
cat("Units noy: '", u_noy, "'\n\n", sep = "")

conv_nhx <- ndep_factor_to_kg_ha(u_nhx)
conv_noy <- ndep_factor_to_kg_ha(u_noy)

cat("Convert nhx factor: ", conv_nhx$f, " -> ", conv_nhx$out, "\n", sep = "")
cat("Convert noy factor: ", conv_noy$f, " -> ", conv_noy$out, "\n\n", sep = "")

r_nhx <- c(rast_var(nhx_1850, "nhx"), rast_var(nhx_1901, "nhx"))
r_noy <- c(rast_var(noy_1850, "noy"), rast_var(noy_1901, "noy"))
stopifnot(terra::nlyr(r_nhx) == terra::nlyr(r_noy))

# ISIMIP time fallback: monthly sequence starting 1850-01
ym_nd <- get_time_ym(r_nhx, start_date = "1850-01-01")

nd_raw <- extract_two_vars_monthly(r_nhx, r_noy, pts, ym_nd, "nhx", "noy", "NDEP")

ndep <- nd_raw %>%
  mutate(
    nhx_kg_ha = nhx * conv_nhx$f,
    noy_kg_ha = noy * conv_noy$f,
    n_total_kg_ha = nhx_kg_ha + noy_kg_ha,
    units_in_nhx = u_nhx,
    units_in_noy = u_noy,
    units_out = conv_nhx$out
  ) %>%
  select(site_id, year, month, nhx_kg_ha, noy_kg_ha, n_total_kg_ha,
         units_in_nhx, units_in_noy, units_out)

cat("\nNdep rows: ", nrow(ndep), " (expected ~", exp_rows, ")\n", sep = "")
cat("Unique site-year-month: ", nrow(distinct(ndep, site_id, year, month)), "\n\n", sep = "")

arrow::write_parquet(ndep, out_ndep, compression = "zstd")
cat("Saved: ", out_ndep, "\n\n", sep = "")

cat("=== 02a DONE ===\n")
