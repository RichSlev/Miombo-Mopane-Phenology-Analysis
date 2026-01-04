---
title: "03_cross-site_synthesis"
author: "Richard Slevin & Eamon O'Cathain"
date: "2025-01-28"
output: html_document
---

# ============================================================
# Script 3: Cross-site synthesis, interpretation, and final figures
# Purpose:  Compare phenology metrics across sites and classes, quantify
#           interannual anomalies, summarize rainfall coupling diagnostics,
#           and generate paper-ready tables and figures.
# Inputs:   results/phenology/phenology_metrics_by_year.csv
#           results/phenology/phenology_metrics_interannual_summary.csv
#           results/phenology/ndvi_rain_lag_summary.csv
# Outputs:  results/reporting/tables/phenology_summary_table.csv
#           results/reporting/tables/site_comparison_table.csv
#           results/reporting/tables/anomaly_table.csv
#           results/reporting/figures/Fig_2_phenology_metrics.png
#           results/reporting/figures/Fig_3_rainfall_coupling.png
# Notes:    Script 3 assumes Script 2 has already been run.
#           No spatial processing. No time-series preprocessing.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  library(purrr)
})

# ----------------------------
# Paths
# ----------------------------
in_dir  <- file.path("results", "phenology")
out_dir <- file.path("results", "reporting")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Inputs (from Script 2)
# ----------------------------
paths <- list(
  metrics_by_year   = file.path(in_dir, "phenology_metrics_by_year.csv"),
  metrics_summary   = file.path(in_dir, "phenology_metrics_interannual_summary.csv"),
  lag_summary       = file.path(in_dir, "ndvi_rain_lag_summary.csv")
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) {
  stop("Missing Script 2 outputs in ", in_dir, ": ", paste(missing, collapse = ", "))
}

metrics_year <- readr::read_csv(paths$metrics_by_year, show_col_types = FALSE) %>%
  mutate(
    sos_date  = as.Date(sos_date),
    eos_date  = as.Date(eos_date),
    peak_date = as.Date(peak_date)
  )

metrics_summary <- readr::read_csv(paths$metrics_summary, show_col_types = FALSE)

lag_df <- readr::read_csv(paths$lag_summary, show_col_types = FALSE)

# ----------------------------
# Optional: class label mapping
# If Script 1 used numeric codes, you can map them to names here.
# Leave as identity if you do not want relabeling.
# ----------------------------
class_map <- tibble::tibble(
  class_raw = character(),
  class_label = character()
)

apply_class_labels <- function(df) {
  if (nrow(class_map) == 0) return(df)
  df %>%
    left_join(class_map, by = c("class" = "class_raw")) %>%
    mutate(class = ifelse(is.na(class_label), class, class_label)) %>%
    select(-class_label)
}

metrics_year    <- apply_class_labels(metrics_year)
metrics_summary <- apply_class_labels(metrics_summary)
lag_df          <- apply_class_labels(lag_df)

# ----------------------------
# Utilities
# ----------------------------
safe_mean <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd   <- function(x) if (sum(is.finite(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)

fmt_mean_sd <- function(mu, sd) {
  if (!is.finite(mu)) return(NA_character_)
  if (!is.finite(sd)) return(sprintf("%.2f", mu))
  sprintf("%.2f ± %.2f", mu, sd)
}

# ============================================================
# 1) Summary table for writing (site x class)
# ============================================================
phenology_summary_table <- metrics_year %>%
  group_by(site, class) %>%
  summarise(
    years_n = sum(is.finite(amp_ndvi)),
    amp_mean = safe_mean(amp_ndvi),
    amp_sd   = safe_sd(amp_ndvi),
    peak_ndvi_mean = safe_mean(peak_ndvi),
    peak_ndvi_sd   = safe_sd(peak_ndvi),
    sos_doy_mean = safe_mean(sos_doy),
    sos_doy_sd   = safe_sd(sos_doy),
    eos_doy_mean = safe_mean(eos_doy),
    eos_doy_sd   = safe_sd(eos_doy),
    los_mean = safe_mean(los_days),
    los_sd   = safe_sd(los_days),
    .groups = "drop"
  ) %>%
  mutate(
    amp_mean_sd       = map2_chr(amp_mean, amp_sd, fmt_mean_sd),
    peak_ndvi_mean_sd = map2_chr(peak_ndvi_mean, peak_ndvi_sd, fmt_mean_sd),
    sos_doy_mean_sd   = map2_chr(sos_doy_mean, sos_doy_sd, fmt_mean_sd),
    eos_doy_mean_sd   = map2_chr(eos_doy_mean, eos_doy_sd, fmt_mean_sd),
    los_mean_sd       = map2_chr(los_mean, los_sd, fmt_mean_sd)
  ) %>%
  select(
    site, class, years_n,
    amp_mean_sd, peak_ndvi_mean_sd, sos_doy_mean_sd, eos_doy_mean_sd, los_mean_sd
  ) %>%
  arrange(site, class)

readr::write_csv(phenology_summary_table, file.path(tab_dir, "phenology_summary_table.csv"))

# ============================================================
# 2) Site comparison table (effect sizes per class)
# This expects both sites exist for a class; otherwise returns NA rows.
# Comparison: (site_B - site_A) where site_A is alphabetically first
# ============================================================
site_comparison_table <- metrics_year %>%
  group_by(site, class) %>%
  summarise(
    amp_mean = safe_mean(amp_ndvi),
    sos_mean = safe_mean(sos_doy),
    eos_mean = safe_mean(eos_doy),
    los_mean = safe_mean(los_days),
    peak_mean = safe_mean(peak_ndvi),
    .groups = "drop"
  ) %>%
  group_by(class) %>%
  group_modify(function(.x, .g) {
    if (nrow(.x) < 2) {
      return(tibble(
        site_A = NA_character_,
        site_B = NA_character_,
        amp_diff = NA_real_,
        sos_diff = NA_real_,
        eos_diff = NA_real_,
        los_diff = NA_real_,
        peak_diff = NA_real_
      ))
    }
    .x <- .x %>% arrange(site)
    A <- .x[1, ]
    B <- .x[2, ]

    tibble(
      site_A = A$site,
      site_B = B$site,
      amp_diff  = B$amp_mean  - A$amp_mean,
      sos_diff  = B$sos_mean  - A$sos_mean,
      eos_diff  = B$eos_mean  - A$eos_mean,
      los_diff  = B$los_mean  - A$los_mean,
      peak_diff = B$peak_mean - A$peak_mean
    )
  }) %>%
  ungroup() %>%
  mutate(class = rep(unique(metrics_year$class), each = 1)) %>%
  select(class, everything()) %>%
  arrange(class)

readr::write_csv(site_comparison_table, file.path(tab_dir, "site_comparison_table.csv"))

# ============================================================
# 3) Interannual anomaly table (site x class x year)
# Anomaly = value - site-class mean
# ============================================================
site_class_means <- metrics_year %>%
  group_by(site, class) %>%
  summarise(
    amp_mu = safe_mean(amp_ndvi),
    sos_mu = safe_mean(sos_doy),
    eos_mu = safe_mean(eos_doy),
    los_mu = safe_mean(los_days),
    peak_mu = safe_mean(peak_ndvi),
    .groups = "drop"
  )

anomaly_table <- metrics_year %>%
  left_join(site_class_means, by = c("site", "class")) %>%
  mutate(
    amp_anom  = amp_ndvi  - amp_mu,
    sos_anom  = sos_doy   - sos_mu,
    eos_anom  = eos_doy   - eos_mu,
    los_anom  = los_days  - los_mu,
    peak_anom = peak_ndvi - peak_mu
  ) %>%
  select(site, class, year, amp_anom, sos_anom, eos_anom, los_anom, peak_anom) %>%
  arrange(site, class, year)

readr::write_csv(anomaly_table, file.path(tab_dir, "anomaly_table.csv"))

# ============================================================
# 4) Rainfall coupling summary for writing
# Adds best lag in days assuming pentad cadence (5 days)
# If your cadence differs, update lag_days_per_step.
# ============================================================
lag_days_per_step <- 5

rain_coupling_table <- lag_df %>%
  mutate(
    best_lag_days = best_lag * lag_days_per_step
  ) %>%
  select(site, class, best_lag, best_lag_days, best_r) %>%
  arrange(site, class)

readr::write_csv(rain_coupling_table, file.path(tab_dir, "rainfall_coupling_table.csv"))

# ============================================================
# Figures (publication-oriented)
# Fig 2: Distributions of key phenology metrics by site and class
# Fig 3: Rainfall coupling best lag and strength by site and class
# ============================================================

# ----------------------------
# Fig 2: Phenology metrics summary
# ----------------------------
fig2_df <- metrics_year %>%
  filter(is.finite(amp_ndvi)) %>%
  mutate(
    class = as.factor(class),
    site  = as.factor(site)
  )

p_fig2_amp <- ggplot(fig2_df, aes(x = class, y = amp_ndvi, fill = site)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.4, na.rm = TRUE) +
  labs(x = "LULC class", y = "Seasonal amplitude (NDVI)", fill = "Site") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fig2_sos <- ggplot(fig2_df, aes(x = class, y = sos_doy, fill = site)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.4, na.rm = TRUE) +
  labs(x = "LULC class", y = "Start of season (DOY)", fill = "Site") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fig2_los <- ggplot(fig2_df, aes(x = class, y = los_days, fill = site)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.4, na.rm = TRUE) +
  labs(x = "LULC class", y = "Length of season (days)", fill = "Site") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save as separate panels (robust for manuscripts and layout tools)
ggsave(file.path(fig_dir, "Fig_2a_amplitude.png"), p_fig2_amp, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_2b_SOS.png"),       p_fig2_sos, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_2c_LOS.png"),       p_fig2_los, width = 10, height = 5.5, dpi = 300)

# Also save a single "Fig 2" placeholder image for workflows that expect one file
# Here: amplitude panel as the canonical Fig 2 output
ggsave(file.path(fig_dir, "Fig_2_phenology_metrics.png"), p_fig2_amp, width = 10, height = 5.5, dpi = 300)

# ----------------------------
# Fig 3: Rainfall coupling diagnostics
# ----------------------------
fig3_df <- rain_coupling_table %>%
  filter(is.finite(best_r)) %>%
  mutate(
    class = as.factor(class),
    site  = as.factor(site)
  )

p_fig3_lag <- ggplot(fig3_df, aes(x = class, y = best_lag_days, fill = site)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, na.rm = TRUE) +
  labs(x = "LULC class", y = "Best rainfall lead (days)", fill = "Site") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_fig3_r <- ggplot(fig3_df, aes(x = class, y = best_r, fill = site)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, na.rm = TRUE) +
  labs(x = "LULC class", y = "Correlation at best lag", fill = "Site") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "Fig_3a_best_lag_days.png"), p_fig3_lag, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "Fig_3b_best_r.png"),        p_fig3_r,   width = 10, height = 5.5, dpi = 300)

# Single "Fig 3" placeholder output
ggsave(file.path(fig_dir, "Fig_3_rainfall_coupling.png"), p_fig3_r, width = 10, height = 5.5, dpi = 300)

# ============================================================
# Final message
# ============================================================
message("Done. Reporting outputs written to: ", out_dir)
message("Tables: ", tab_dir)
message("Figures: ", fig_dir)
