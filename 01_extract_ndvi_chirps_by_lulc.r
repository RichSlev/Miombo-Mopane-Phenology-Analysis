---
title: "01_extract_ndvi_chirps_by_lulc"
author: "Richard Slevin & Eamon O'Cathain"
date: "2025-01-28"
output: html_document
---

# ============================================================
# Script 1: Extract class-level NDVI and CHIRPS time series
# Purpose:  Compute per-date mean NDVI and rainfall per LULC class
# Inputs:   rasters/NDVI.tif (multi-layer)
#           rasters/CHIRPS_pentad.tif (multi-layer)
#           rasters/LULC_MCC.tif
#           Rasters/detailled_LULC_Kariba_area.tif
#           Shapefile/MCC.shp
#           Shapefile/Zone_etude_Kariba_lake.shp
# Outputs:  results/timeseries/ndvi_mcc_means.csv
#           results/timeseries/chirps_mcc_means.csv
#           results/timeseries/ndvi_kariba_means.csv
#           results/timeseries/chirps_kariba_means.csv
#           (+ equivalent .rds files)
# Notes:    CSV schema: date (YYYY-MM-DD) + one column per LULC class
#           NDVI is exported as scaled integer (approx -2000..10000)
#           CHIRPS is exported as mean rainfall per pentad (raster units)
# ============================================================


## Purpose

This notebook:
1. Loads LULC, NDVI, and CHIRPS rasters plus study-area vectors
2. Aligns all rasters to the NDVI grid (CRS, extent, resolution)
3. Crops and masks to the study area
4. Extracts class-level mean time series for NDVI and CHIRPS for:
   - Kariba LULC classes
   - MCC LULC classes (including detailed Miombo and Mopane sub-classes)
5. Writes clean CSV outputs for downstream plotting and phenology analysis

---

## Setup

```{r, setup}
rm(list = ls())

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(lubridate)

terraOptions(progress = 1)


```



```{r paths}

# ---------------------------
# Relative paths (repo root)
# ---------------------------
dir_rasters   <- "./rasters"
dir_rasters2  <- "./Rasters"      # keep both because repo uses both
dir_vectors   <- "./Shapefile"
dir_results   <- "./results/timeseries"

# ---------------------------
# Output directory
# ---------------------------
dir.create(dir_results, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Input file paths
# ---------------------------
fp_lulc_mcc    <- file.path(dir_rasters,  "LULC_MCC.tif")
fp_lulc_kariba <- file.path(dir_rasters2, "detailled_LULC_Kariba_area.tif")
fp_ndvi        <- file.path(dir_rasters,  "NDVI.tif")
fp_chirps      <- file.path(dir_rasters,  "CHIRPS_pentad.tif")

fp_study_area      <- file.path(dir_vectors, "Zone_etude_Kariba_lake.shp")
fp_mcc             <- file.path(dir_vectors, "MCC.shp")
fp_protected_areas <- file.path(dir_vectors, "protected_areas.shp")

# ---------------------------
# Quick existence check
# ---------------------------
inputs <- c(
  fp_lulc_mcc, fp_lulc_kariba, fp_ndvi, fp_chirps,
  fp_study_area, fp_mcc, fp_protected_areas
)

missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing input files:\n- ", paste(missing_inputs, collapse = "\n- "))
}


```


```{r load-data}

# ---------------------------
# Load rasters
# ---------------------------
lulc_mcc    <- rast(fp_lulc_mcc)
lulc_kariba <- rast(fp_lulc_kariba)
ndvi        <- rast(fp_ndvi)
chirps      <- rast(fp_chirps)

# ---------------------------
# Load vectors
# ---------------------------
study_area      <- vect(fp_study_area)
mcc             <- vect(fp_mcc)
protected_areas <- vect(fp_protected_areas)

# ---------------------------
# Basic sanity checks
# ---------------------------
stopifnot(inherits(lulc_mcc, "SpatRaster"))
stopifnot(inherits(lulc_kariba, "SpatRaster"))
stopifnot(inherits(ndvi, "SpatRaster"))
stopifnot(inherits(chirps, "SpatRaster"))

stopifnot(inherits(study_area, "SpatVector"))
stopifnot(inherits(mcc, "SpatVector"))
stopifnot(inherits(protected_areas, "SpatVector"))



```




```{r}

# ---------------------------
# Reproject vectors to NDVI CRS
# ---------------------------
study_area      <- project(study_area, ndvi)
mcc             <- project(mcc, ndvi)
protected_areas <- project(protected_areas, ndvi)

# ---------------------------
# Sanity check: CRS alignment
# ---------------------------
stopifnot(same.crs(study_area, ndvi))
stopifnot(same.crs(mcc, ndvi))
stopifnot(same.crs(protected_areas, ndvi))

```


```{r align rasters}

# --------------------------------------------------
# Helper to align rasters to NDVI grid, then crop/mask
# --------------------------------------------------
align_crop_mask <- function(x, ref, mask_vect, method = c("near", "bilinear")) {
  method <- match.arg(method)

  # Reproject to reference CRS if needed
  if (!same.crs(x, ref)) {
    x <- project(x, ref, method = method)
  }

  # Resample to reference grid if resolution or extent differ
  if (!isTRUE(all.equal(res(x), res(ref))) ||
      !isTRUE(all.equal(ext(x), ext(ref)))) {
    x <- resample(x, ref, method = method)
  }

  # Crop and mask to study area
  x <- crop(x, mask_vect)
  x <- mask(x, mask_vect)

  return(x)
}


# --------------------------------------------------
# Align, crop, and mask rasters to the NDVI grid
# --------------------------------------------------

# NDVI is the reference grid (continuous)
ndvi_masked <- align_crop_mask(
  x = ndvi,
  ref = ndvi,
  mask_vect = study_area,
  method = "bilinear"
)

# CHIRPS is continuous (rainfall)
chirps_masked <- align_crop_mask(
  x = chirps,
  ref = ndvi,
  mask_vect = study_area,
  method = "bilinear"
)

# LULC rasters are categorical
lulc_mcc_masked <- align_crop_mask(
  x = lulc_mcc,
  ref = ndvi,
  mask_vect = study_area,
  method = "near"
)

lulc_kariba_masked <- align_crop_mask(
  x = lulc_kariba,
  ref = ndvi,
  mask_vect = study_area,
  method = "near"
)

# --------------------------------------------------
# Sanity checks
# --------------------------------------------------
stopifnot(same.crs(ndvi_masked, chirps_masked))
stopifnot(same.crs(ndvi_masked, lulc_mcc_masked))
stopifnot(same.crs(ndvi_masked, lulc_kariba_masked))

stopifnot(all.equal(res(ndvi_masked), res(chirps_masked)))
stopifnot(all.equal(res(ndvi_masked), res(lulc_mcc_masked)))
stopifnot(all.equal(res(ndvi_masked), res(lulc_kariba_masked)))



```


```{r classes-kariba}

# --------------------------------------------------
# Define land-cover indices: Kariba LULC classes
# --------------------------------------------------

miombo_kariba   <- which(lulc_kariba_masked[] == 4)
mopane_kariba   <- which(lulc_kariba_masked[] == 5)
riparian_kariba <- which(lulc_kariba_masked[] == 7)
crop_kariba     <- which(lulc_kariba_masked[] == 2)
all_kariba      <- which(!is.na(values(lulc_kariba_masked)))

list_kariba <- list(
  miombo_kariba   = miombo_kariba,
  mopane_kariba   = mopane_kariba,
  riparian_kariba = riparian_kariba,
  crop_kariba     = crop_kariba,
  all_kariba      = all_kariba
)

# --------------------------------------------------
# Sanity checks
# --------------------------------------------------
stopifnot(length(all_kariba) > 0)


```


```{r classes-mcc}

# --------------------------------------------------
# Define land-cover indices: MCC LULC classes
# --------------------------------------------------

miombo_mcc_all <- which(lulc_mcc_masked[] %in% c(41, 42, 43))
miombo_mcc_41  <- which(lulc_mcc_masked[] == 41)
miombo_mcc_42  <- which(lulc_mcc_masked[] == 42)
miombo_mcc_43  <- which(lulc_mcc_masked[] == 43)

mopane_mcc_all <- which(lulc_mcc_masked[] %in% c(51, 52, 53))
mopane_mcc_51  <- which(lulc_mcc_masked[] == 51)
mopane_mcc_52  <- which(lulc_mcc_masked[] == 52)
mopane_mcc_53  <- which(lulc_mcc_masked[] == 53)

riparian_mcc <- which(lulc_mcc_masked[] == 7)
crop_mcc     <- which(lulc_mcc_masked[] == 2)
all_mcc      <- which(!is.na(values(lulc_mcc_masked)))

list_mcc <- list(
  miombo_mcc_all = miombo_mcc_all,
  miombo_mcc_41  = miombo_mcc_41,
  miombo_mcc_42  = miombo_mcc_42,
  miombo_mcc_43  = miombo_mcc_43,
  mopane_mcc_all = mopane_mcc_all,
  mopane_mcc_51  = mopane_mcc_51,
  mopane_mcc_52  = mopane_mcc_52,
  mopane_mcc_53  = mopane_mcc_53,
  riparian_mcc   = riparian_mcc,
  crop_mcc       = crop_mcc,
  all_mcc        = all_mcc
)

# --------------------------------------------------
# Sanity checks
# --------------------------------------------------
stopifnot(length(all_mcc) > 0)


```


```{r dates-chirps}
# --------------------------------------------------
# Parse CHIRPS dates from raster layer names
# Expected format: YYYYMMDD_precipitation
# --------------------------------------------------

chirps_layer_names <- names(chirps_masked)

# Remove the literal suffix "_precipitation"
chirps_dates_chr <- sub("_precipitation$", "", chirps_layer_names)

# Convert to Date
chirps_dates <- as.Date(chirps_dates_chr, format = "%Y%m%d")

# --------------------------------------------------
# Sanity checks
# --------------------------------------------------
if (any(is.na(chirps_dates))) {
  stop("CHIRPS date parsing failed: NA values produced. Check layer name format.")
}

stopifnot(length(chirps_dates) == nlyr(chirps_masked))


```


```{r extract-timeseries}

# --------------------------------------------------
# Helper to extract class-level mean time series
# --------------------------------------------------
extract_means <- function(index_list, raster_stack, dates) {

  # Pull raster values once (matrix: ncell x nlyr)
  v <- terra::values(raster_stack, mat = TRUE)

  # Initialise output dataframe
  out <- data.frame(date = dates)

  for (nm in names(index_list)) {
    idx <- index_list[[nm]]

    if (length(idx) == 0) {
      out[[nm]] <- NA_real_
    } else {
      # Subset rows (pixels) for this class, then average across pixels
      out[[nm]] <- colMeans(v[idx, , drop = FALSE], na.rm = TRUE)
    }
  }

  return(out)
}


# --------------------------------------------------
# Extract class-mean time series for NDVI and CHIRPS
# --------------------------------------------------

# NDVI
ndvi_mcc_df    <- extract_means(list_mcc, ndvi_masked, ndvi_dates)
ndvi_kariba_df <- extract_means(list_kariba, ndvi_masked, ndvi_dates)

# CHIRPS
chirps_mcc_df    <- extract_means(list_mcc, chirps_masked, chirps_dates)
chirps_kariba_df <- extract_means(list_kariba, chirps_masked, chirps_dates)


```


```{r validate-sort}

# --------------------------------------------------
# Order, validate, and sanity-check extracted time series
# --------------------------------------------------

# Ensure time series are ordered by date
ndvi_mcc_df    <- ndvi_mcc_df    |> arrange(date)
ndvi_kariba_df <- ndvi_kariba_df |> arrange(date)

chirps_mcc_df    <- chirps_mcc_df    |> arrange(date)
chirps_kariba_df <- chirps_kariba_df |> arrange(date)

# --------------------------------------------------
# Sanity checks: monotonic dates and expected lengths
# --------------------------------------------------
stopifnot(all(diff(ndvi_mcc_df$date)    > 0))
stopifnot(all(diff(ndvi_kariba_df$date) > 0))
stopifnot(all(diff(chirps_mcc_df$date)  > 0))
stopifnot(all(diff(chirps_kariba_df$date) > 0))

stopifnot(nrow(ndvi_mcc_df)    == length(ndvi_dates))
stopifnot(nrow(ndvi_kariba_df) == length(ndvi_dates))
stopifnot(nrow(chirps_mcc_df)  == length(chirps_dates))
stopifnot(nrow(chirps_kariba_df) == length(chirps_dates))


```


```{r write-outputs}

# --------------------------------------------------
# Write outputs for downstream scripts (Script 2, Script 3)
# --------------------------------------------------

# Output directory (already created in chunk 2, but ensure it exists)
out_dir <- "./results/timeseries"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# CSV outputs
write_csv(ndvi_mcc_df,     file.path(out_dir, "ndvi_mcc_means.csv"))
write_csv(chirps_mcc_df,   file.path(out_dir, "chirps_mcc_means.csv"))
write_csv(ndvi_kariba_df,  file.path(out_dir, "ndvi_kariba_means.csv"))
write_csv(chirps_kariba_df,file.path(out_dir, "chirps_kariba_means.csv"))

# RDS outputs (preserve types, faster reload)
saveRDS(ndvi_mcc_df,      file.path(out_dir, "ndvi_mcc_means.rds"))
saveRDS(chirps_mcc_df,    file.path(out_dir, "chirps_mcc_means.rds"))
saveRDS(ndvi_kariba_df,   file.path(out_dir, "ndvi_kariba_means.rds"))
saveRDS(chirps_kariba_df, file.path(out_dir, "chirps_kariba_means.rds"))

message("All time-series outputs written to: ", normalizePath(out_dir))


```



---

## Output schema (for Script 2)

### Files produced
All outputs are written to:

- `results/timeseries/ndvi_mcc_means.csv`
- `results/timeseries/chirps_mcc_means.csv`
- `results/timeseries/ndvi_kariba_means.csv`
- `results/timeseries/chirps_kariba_means.csv`  
(and equivalent `.rds` files)

### Shared structure (all CSVs)
- **Rows:** one row per raster band (time step)
- **Column `date`:**
  - Type: `YYYY-MM-DD`
  - Parsed with `as.Date(date)`
- **Remaining columns:** one column per LULC class
  - Numeric (double)
  - Value = spatial mean across all pixels of that class
  - `NA` if class has

---

