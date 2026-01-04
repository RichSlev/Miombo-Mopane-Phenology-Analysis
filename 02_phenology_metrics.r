---
title: "02_phenology_metrics"
author: "Richard Slevin & Eamon O'Cathain"
date: "2025-01-28"
output: html_document
---


# ============================================================
# Script 2: Phenology metrics 
# Inputs:  results/timeseries/*_means.csv (NDVI and CHIRPS)
# Outputs: results/phenology/phenology_metrics_*.csv (+ plots)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  library(purrr)
  library(zoo)
})

# ----------------------------
# Paths
# ----------------------------
data_dir <- file.path("results", "timeseries")
out_dir  <- file.path("results", "phenology")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------------------------
# Helper: read Script 1 outputs
# ----------------------------
read_class_timeseries <- function(path, variable = c("ndvi", "chirps")) {
  variable <- match.arg(variable)

  df <- readr::read_csv(path, show_col_types = FALSE)

  if (!("date" %in% names(df))) {
    stop("Expected a 'date' column in: ", path)
  }

  df <- df %>%
    mutate(date = as.Date(date)) %>%
    arrange(date)

  # Identify site from filename: ndvi_kariba_means.csv -> kariba
  site <- basename(path) %>%
    str_remove("\\.csv$") %>%
    str_remove(paste0("^", variable, "_")) %>%
    str_remove("_means$")

  attr(df, "site") <- site
  attr(df, "variable") <- variable
  df
}

# ----------------------------
# NDVI rescaling
# Script 1 notes: NDVI values are scaled integers approx -2000..10000
# Default: scale factor 1/10000
# ----------------------------
rescale_ndvi <- function(df, scale_factor = 10000) {
  class_cols <- setdiff(names(df), "date")

  df %>%
    mutate(across(all_of(class_cols), ~ .x / scale_factor))
}

# ----------------------------
# Pivot to long format
# Returns: date, class, value, site, variable
# ----------------------------
to_long <- function(df) {
  site     <- attr(df, "site")
  variable <- attr(df, "variable")

  df %>%
    pivot_longer(cols = -date, names_to = "class", values_to = "value") %>%
    mutate(
      site = site,
      variable = variable
    )
}

# ----------------------------
# Regularize time axis by site-variable-class
# Inserts missing dates (pentads) explicitly as NA
# If cadence is pentad, expected spacing is 5 days, but we infer from data.
# ----------------------------
regularize_dates <- function(df_long) {
  df_long %>%
    group_by(site, variable, class) %>%
    arrange(date, .by_group = TRUE) %>%
    group_modify(function(.x, .g) {
      dates <- .x$date
      if (length(dates) < 2) return(.x)

      # infer step as the modal difference (days)
      dd <- diff(dates) %>% as.numeric()
      step <- suppressWarnings(as.numeric(stats::median(dd, na.rm = TRUE)))
      if (!is.finite(step) || step <= 0) step <- 5

      full_dates <- seq(min(dates, na.rm = TRUE), max(dates, na.rm = TRUE), by = step)
      tibble(date = full_dates) %>%
        left_join(.x %>% select(date, value), by = "date") %>%
        mutate(
          site = .g$site,
          variable = .g$variable,
          class = .g$class
        )
    }) %>%
    ungroup()
}

# ----------------------------
# Smoothing
# Rolling mean over k timesteps. Pentad data: k=3 is 15 days, k=5 is 25 days.
# Keep raw values intact.
# ----------------------------
add_smoothing <- function(df_long, k = 5) {
  df_long %>%
    group_by(site, variable, class) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      value_smooth = zoo::rollapply(value, width = k, FUN = mean, align = "center", fill = NA, na.rm = TRUE)
    ) %>%
    ungroup()
}

# ----------------------------
# Phenology extraction (per site-class-year) for NDVI only
# Method: amplitude-fraction threshold on smoothed NDVI
# Definitions:
#  base = 10th percentile (smoothed) within year
#  peak = max (smoothed) within year
#  amp  = peak - base
#  thresh = base + frac * amp
#  SOS = first date value_smooth >= thresh
#  EOS = last date value_smooth >= thresh
# Peak date = date at max value_smooth
# Guardrails:
#  - require enough non-NA points
#  - require amp > min_amp
# ----------------------------
phenology_metrics_year <- function(df_ndvi, frac = 0.2, min_points = 20, min_amp = 0.05) {
  df_ndvi %>%
    mutate(year = lubridate::year(date)) %>%
    group_by(site, class, year) %>%
    group_modify(function(.x, .g) {
      x <- .x %>% arrange(date)

      v <- x$value_smooth
      ok <- is.finite(v)

      if (sum(ok) < min_points) {
        return(tibble(
          n_points = sum(ok),
          base_ndvi = NA_real_,
          peak_ndvi = NA_real_,
          amp_ndvi  = NA_real_,
          thresh_ndvi = NA_real_,
          sos_date = as.Date(NA),
          eos_date = as.Date(NA),
          peak_date = as.Date(NA),
          sos_doy = NA_integer_,
          eos_doy = NA_integer_,
          peak_doy = NA_integer_,
          los_days = NA_real_
        ))
      }

      v_ok <- v[ok]
      base <- as.numeric(stats::quantile(v_ok, probs = 0.10, na.rm = TRUE, type = 7))
      peak <- max(v_ok, na.rm = TRUE)
      amp  <- peak - base

      if (!is.finite(amp) || amp < min_amp) {
        return(tibble(
          n_points = sum(ok),
          base_ndvi = base,
          peak_ndvi = peak,
          amp_ndvi  = amp,
          thresh_ndvi = NA_real_,
          sos_date = as.Date(NA),
          eos_date = as.Date(NA),
          peak_date = as.Date(NA),
          sos_doy = NA_integer_,
          eos_doy = NA_integer_,
          peak_doy = NA_integer_,
          los_days = NA_real_
        ))
      }

      thresh <- base + frac * amp
      above  <- ok & (v >= thresh)

      if (!any(above, na.rm = TRUE)) {
        return(tibble(
          n_points = sum(ok),
          base_ndvi = base,
          peak_ndvi = peak,
          amp_ndvi  = amp,
          thresh_ndvi = thresh,
          sos_date = as.Date(NA),
          eos_date = as.Date(NA),
          peak_date = as.Date(NA),
          sos_doy = NA_integer_,
          eos_doy = NA_integer_,
          peak_doy = NA_integer_,
          los_days = NA_real_
        ))
      }

      sos_idx <- which(above)[1]
      eos_idx <- tail(which(above), 1)

      peak_idx <- which.max(ifelse(ok, v, -Inf))

      sos_date  <- x$date[sos_idx]
      eos_date  <- x$date[eos_idx]
      peak_date <- x$date[peak_idx]

      tibble(
        n_points = sum(ok),
        base_ndvi = base,
        peak_ndvi = peak,
        amp_ndvi  = amp,
        thresh_ndvi = thresh,
        sos_date = sos_date,
        eos_date = eos_date,
        peak_date = peak_date,
        sos_doy = lubridate::yday(sos_date),
        eos_doy = lubridate::yday(eos_date),
        peak_doy = lubridate::yday(peak_date),
        los_days = as.numeric(eos_date - sos_date)
      )
    }) %>%
    ungroup()
}

# ----------------------------
# Interannual variability summary (per site-class)
# ----------------------------
summarise_interannual <- function(metrics_year) {
  metrics_year %>%
    group_by(site, class) %>%
    summarise(
      years_n = sum(is.finite(amp_ndvi)),
      amp_mean = mean(amp_ndvi, na.rm = TRUE),
      amp_sd   = sd(amp_ndvi, na.rm = TRUE),
      peak_ndvi_mean = mean(peak_ndvi, na.rm = TRUE),
      peak_ndvi_sd   = sd(peak_ndvi, na.rm = TRUE),
      sos_doy_mean = mean(sos_doy, na.rm = TRUE),
      sos_doy_sd   = sd(sos_doy, na.rm = TRUE),
      eos_doy_mean = mean(eos_doy, na.rm = TRUE),
      eos_doy_sd   = sd(eos_doy, na.rm = TRUE),
      los_mean = mean(los_days, na.rm = TRUE),
      los_sd   = sd(los_days, na.rm = TRUE),
      .groups = "drop"
    )
}

# ----------------------------
# Lagged correlation NDVI vs CHIRPS (per site-class)
# Uses smoothed NDVI and raw CHIRPS (or could smooth CHIRPS similarly)
# Lags in timesteps (pentads). Default 0..8 pentads (0..40 days).
# Returns best lag and r.
# ----------------------------
lagged_correlation <- function(ndvi_long, chirps_long, max_lag = 8) {
  # Join by date, site, class
  df <- ndvi_long %>%
    filter(variable == "ndvi") %>%
    select(site, class, date, ndvi = value_smooth) %>%
    inner_join(
      chirps_long %>%
        filter(variable == "chirps") %>%
        select(site, class, date, rain = value),
      by = c("site", "class", "date")
    ) %>%
    group_by(site, class) %>%
    group_modify(function(.x, .g) {
      x <- .x %>% arrange(date)

      # compute correlations for lags: rain leads NDVI by lag
      res <- map_dfr(0:max_lag, function(lg) {
        rain_lag <- dplyr::lag(x$rain, n = lg)
        ok <- is.finite(x$ndvi) & is.finite(rain_lag)
        r <- if (sum(ok) >= 10) stats::cor(x$ndvi[ok], rain_lag[ok]) else NA_real_
        tibble(lag = lg, r = r)
      })

      best <- res %>%
        filter(is.finite(r)) %>%
        arrange(desc(abs(r))) %>%
        slice(1)

      tibble(
        best_lag = ifelse(nrow(best) == 0, NA_integer_, best$lag),
        best_r   = ifelse(nrow(best) == 0, NA_real_, best$r)
      ) %>%
        mutate(curve = list(res))
    }) %>%
    ungroup()

  df
}

# ----------------------------
# Plots
# ----------------------------
plot_timeseries <- function(df_long, variable_label, out_path) {
  p <- ggplot(df_long, aes(x = date, y = value, colour = class)) +
    geom_line(linewidth = 0.5, na.rm = TRUE) +
    facet_wrap(~ site, scales = "free_x") +
    labs(x = NULL, y = variable_label, colour = "Class") +
    theme_minimal(base_size = 12)

  ggsave(out_path, p, width = 11, height = 6, dpi = 300)
}

plot_ndvi_with_markers <- function(ndvi_long, metrics_year, out_path) {
  markers <- metrics_year %>%
    select(site, class, year, sos_date, eos_date, peak_date)

  df <- ndvi_long %>%
    filter(variable == "ndvi") %>%
    mutate(year = lubridate::year(date))

  p <- ggplot(df, aes(x = date, y = value_smooth, colour = class)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    geom_vline(
      data = markers,
      aes(xintercept = as.numeric(sos_date), colour = class),
      linetype = "dashed",
      alpha = 0.35,
      show.legend = FALSE
    ) +
    geom_vline(
      data = markers,
      aes(xintercept = as.numeric(eos_date), colour = class),
      linetype = "dashed",
      alpha = 0.35,
      show.legend = FALSE
    ) +
    geom_vline(
      data = markers,
      aes(xintercept = as.numeric(peak_date), colour = class),
      linetype = "solid",
      alpha = 0.25,
      show.legend = FALSE
    ) +
    facet_grid(site ~ year, scales = "free_x") +
    labs(x = NULL, y = "NDVI (scaled to 0..1)", colour = "Class") +
    theme_minimal(base_size = 11)

  ggsave(out_path, p, width = 14, height = 8, dpi = 300)
}

plot_lag_curve <- function(lag_df, out_path) {
  curve_df <- lag_df %>%
    select(site, class, curve) %>%
    tidyr::unnest(curve)

  p <- ggplot(curve_df, aes(x = lag, y = r, colour = class)) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    facet_wrap(~ site) +
    labs(x = "Rainfall lead (pentads)", y = "Correlation (NDVI vs lagged rain)", colour = "Class") +
    theme_minimal(base_size = 12)

  ggsave(out_path, p, width = 11, height = 6, dpi = 300)
}

# ============================================================
# Main pipeline
# ============================================================

# Input files (Script 1 outputs)
files <- list(
  ndvi_kariba  = file.path(data_dir, "ndvi_kariba_means.csv"),
  chirps_kariba = file.path(data_dir, "chirps_kariba_means.csv"),
  ndvi_mcc     = file.path(data_dir, "ndvi_mcc_means.csv"),
  chirps_mcc    = file.path(data_dir, "chirps_mcc_means.csv")
)

missing <- names(files)[!file.exists(unlist(files))]
if (length(missing) > 0) {
  stop("Missing expected Script 1 outputs in ", data_dir, ": ", paste(missing, collapse = ", "))
}

# Load
ndvi_kariba  <- read_class_timeseries(files$ndvi_kariba,  variable = "ndvi")  %>% rescale_ndvi()
chirps_kariba <- read_class_timeseries(files$chirps_kariba, variable = "chirps")

ndvi_mcc     <- read_class_timeseries(files$ndvi_mcc,     variable = "ndvi") %>% rescale_ndvi()
chirps_mcc    <- read_class_timeseries(files$chirps_mcc,    variable = "chirps")

# Long format
ndvi_long   <- bind_rows(to_long(ndvi_kariba), to_long(ndvi_mcc))   %>% regularize_dates() %>% add_smoothing(k = 5)
chirps_long <- bind_rows(to_long(chirps_kariba), to_long(chirps_mcc)) %>% regularize_dates()

# Phenology metrics per year
metrics_year <- phenology_metrics_year(
  df_ndvi = ndvi_long %>% filter(variable == "ndvi"),
  frac = 0.2,
  min_points = 20,
  min_amp = 0.05
)

# Interannual summaries
metrics_interannual <- summarise_interannual(metrics_year)

# Lagged coupling (optional but included)
lag_df <- lagged_correlation(
  ndvi_long = ndvi_long,
  chirps_long = chirps_long,
  max_lag = 8
)

# Write outputs
readr::write_csv(metrics_year, file.path(out_dir, "phenology_metrics_by_year.csv"))
readr::write_csv(metrics_interannual, file.path(out_dir, "phenology_metrics_interannual_summary.csv"))
readr::write_csv(lag_df %>% select(-curve), file.path(out_dir, "ndvi_rain_lag_summary.csv"))

saveRDS(metrics_year, file.path(out_dir, "phenology_metrics_by_year.rds"))
saveRDS(metrics_interannual, file.path(out_dir, "phenology_metrics_interannual_summary.rds"))
saveRDS(lag_df, file.path(out_dir, "ndvi_rain_lag_full.rds"))

# Plots
plot_timeseries(
  df_long = ndvi_long %>% filter(variable == "ndvi"),
  variable_label = "NDVI (scaled to 0..1)",
  out_path = file.path(out_dir, "ts_ndvi_by_class.png")
)

plot_timeseries(
  df_long = chirps_long %>% filter(variable == "chirps"),
  variable_label = "CHIRPS rainfall (pentad mean, raster units)",
  out_path = file.path(out_dir, "ts_chirps_by_class.png")
)

plot_ndvi_with_markers(
  ndvi_long = ndvi_long,
  metrics_year = metrics_year,
  out_path = file.path(out_dir, "ndvi_with_sos_eos_peak.png")
)

plot_lag_curve(
  lag_df = lag_df,
  out_path = file.path(out_dir, "ndvi_rain_lag_curves.png")
)

message("Done. Outputs written to: ", out_dir)
