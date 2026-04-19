#!/usr/bin/env Rscript
# run_unified_analysis.R — Full analysis of the unified simulation
#
# 1. Aggregate scenario-level means (including AC1)
# 2. Re-derive GRASS thresholds from the larger grid
# 3. Apply thresholds and evaluate classification accuracy
# 4. Analyze by regime (prevalence-dominant vs bias-dominant)
# 5. Generate figures

source("R/00_packages.R")
library(data.table)
library(ggplot2)

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
cat("Loading unified simulation results...\n")
results <- readRDS("output/unified_sim/sim_results/all_results.rds")
grid    <- as.data.table(readRDS("output/unified_sim/parameter_grid.rds"))

# Add regime to grid
grid[, regime := ifelse(PI_theoretical > BI_theoretical, "prevalence_dominant",
                 ifelse(BI_theoretical > PI_theoretical, "bias_dominant", "balanced"))]

# ------------------------------------------------------------------
# 2. Aggregate scenario-level means
# ------------------------------------------------------------------
cat("Aggregating...\n")
scenario_means <- results[, .(
  kappa_mean = mean(kappa, na.rm = TRUE),
  PABAK_mean = mean(PABAK),
  AC1_mean   = mean(AC1, na.rm = TRUE),
  P0_mean    = mean(P0),
  pos_agree_mean = mean(pos_agreement, na.rm = TRUE),
  neg_agree_mean = mean(neg_agreement, na.rm = TRUE),
  prev_idx_mean  = mean(prevalence_index),
  bias_idx_mean  = mean(bias_index),
  kappa_mcse = sd(kappa, na.rm = TRUE) / sqrt(sum(!is.na(kappa))),
  kappa_na_pct = 100 * sum(is.na(kappa)) / .N
), by = scenario_id]

scenario_means <- merge(scenario_means, grid, by = "scenario_id")

# ------------------------------------------------------------------
# 3. Re-derive GRASS thresholds via ROC optimization
# ------------------------------------------------------------------
cat("Deriving thresholds...\n")
scenario_means[, is_high := quality_band == "high"]

# For each prevalence stratum, find threshold maximizing Youden's J
derive_threshold <- function(metric_vals, labels) {
  if (length(unique(labels)) < 2) return(list(threshold = NA, J = NA))
  candidates <- sort(unique(metric_vals))
  best_j <- -Inf
  best_t <- NA
  for (t in candidates) {
    pred <- metric_vals >= t
    se <- sum(pred & labels) / max(sum(labels), 1)
    sp <- sum(!pred & !labels) / max(sum(!labels), 1)
    j <- se + sp - 1
    if (j > best_j) { best_j <- j; best_t <- t }
  }
  list(threshold = best_t, J = best_j)
}

prevalences <- sort(unique(scenario_means$prevalence))
thresh_list <- list()
for (p in prevalences) {
  sub <- scenario_means[prevalence == p]
  k <- derive_threshold(sub$kappa_mean, sub$is_high)
  pb <- derive_threshold(sub$PABAK_mean, sub$is_high)
  ac <- derive_threshold(sub$AC1_mean, sub$is_high)
  thresh_list[[as.character(p)]] <- data.table(
    prevalence = p,
    kappa_thresh = k$threshold, J_kappa = k$J,
    PABAK_thresh = pb$threshold, J_PABAK = pb$J,
    AC1_thresh = ac$threshold, J_AC1 = ac$J
  )
}
thresholds <- rbindlist(thresh_list)

cat("\nRe-derived GRASS thresholds:\n")
print(thresholds, digits = 3)

# ------------------------------------------------------------------
# 4. Apply thresholds and classify
# ------------------------------------------------------------------
scenario_means <- merge(scenario_means, thresholds[, .(prevalence, kappa_thresh, PABAK_thresh, AC1_thresh)],
                        by = "prevalence")

scenario_means[, ground_truth := quality_band == "high"]
scenario_means[, kappa_classifies_high := kappa_mean >= kappa_thresh]
scenario_means[, PABAK_classifies_high := PABAK_mean >= PABAK_thresh]
scenario_means[, AC1_classifies_high := AC1_mean >= AC1_thresh]
scenario_means[, kappa_correct := kappa_classifies_high == ground_truth]
scenario_means[, PABAK_correct := PABAK_classifies_high == ground_truth]
scenario_means[, AC1_correct := AC1_classifies_high == ground_truth]

# ------------------------------------------------------------------
# 5. Summary statistics
# ------------------------------------------------------------------
cat("\n=== Overall Classification Accuracy ===\n")
cat(sprintf("  kappa: %.1f%%\n", 100 * mean(scenario_means$kappa_correct)))
cat(sprintf("  PABAK: %.1f%%\n", 100 * mean(scenario_means$PABAK_correct)))
cat(sprintf("  AC1:   %.1f%%\n", 100 * mean(scenario_means$AC1_correct)))

cat("\n=== By Regime ===\n")
regime_acc <- scenario_means[, .(
  n = .N,
  kappa = mean(kappa_correct),
  PABAK = mean(PABAK_correct),
  AC1 = mean(AC1_correct)
), by = regime]
print(regime_acc)

cat("\n=== By Regime and Quality Band ===\n")
regime_qual <- scenario_means[, .(
  n = .N,
  kappa = mean(kappa_correct),
  PABAK = mean(PABAK_correct),
  AC1 = mean(AC1_correct)
), by = .(regime, quality_band)]
print(regime_qual[order(regime, quality_band)])

cat("\n=== Metric Agreement ===\n")
scenario_means[, all_agree := (kappa_classifies_high == PABAK_classifies_high) &
                               (PABAK_classifies_high == AC1_classifies_high)]
discord <- scenario_means[, .(
  n = .N,
  agree = sum(all_agree),
  discord = sum(!all_agree),
  pct_agree = 100 * mean(all_agree)
), by = regime]
print(discord)

cat("\n=== Grid Composition (for transparency) ===\n")
cat(sprintf("Total scenarios: %d\n", nrow(scenario_means)))
cat(sprintf("Prevalence-dominant: %d (%.1f%%)\n",
    sum(scenario_means$regime == "prevalence_dominant"),
    100 * mean(scenario_means$regime == "prevalence_dominant")))
cat(sprintf("Bias-dominant: %d (%.1f%%)\n",
    sum(scenario_means$regime == "bias_dominant"),
    100 * mean(scenario_means$regime == "bias_dominant")))
cat(sprintf("Balanced: %d (%.1f%%)\n",
    sum(scenario_means$regime == "balanced"),
    100 * mean(scenario_means$regime == "balanced")))

# Landis-Koch label accuracy (for the "fixed labels" finding)
lk_label <- function(x) {
  cut(x, breaks = c(-Inf, 0, 0.20, 0.40, 0.60, 0.80, Inf),
      labels = c("poor", "slight", "fair", "moderate", "substantial", "almost_perfect"),
      right = TRUE)
}

lk_to_band <- function(lk) {
  ifelse(lk %in% c("substantial", "almost_perfect"), "high",
  ifelse(lk %in% c("moderate", "fair"), "medium", "low"))
}

scenario_means[, lk_kappa := lk_label(kappa_mean)]
scenario_means[, lk_kappa_band := lk_to_band(lk_kappa)]
scenario_means[, lk_correct := lk_kappa_band == quality_band]

cat(sprintf("\n=== Landis-Koch Classification Accuracy ===\n"))
cat(sprintf("Overall: %.1f%%\n", 100 * mean(scenario_means$lk_correct)))
lk_by_regime <- scenario_means[, .(
  accuracy = mean(lk_correct),
  n = .N
), by = regime]
print(lk_by_regime)

# ------------------------------------------------------------------
# 6. Figures
# ------------------------------------------------------------------
cat("\nGenerating figures...\n")

# Clean profile names
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("(^|\\s)(\\w)", "\\1\\U\\2", x, perl = TRUE)
  x
}

profile_order <- unique(scenario_means[, .(profile_name, quality_band, min_accuracy)])
profile_order <- profile_order[order(
  -factor(quality_band, levels = c("high", "medium", "low")),
  -min_accuracy, profile_name
)]
profile_order[, display_name := paste0(clean_name(profile_name), "  [",
                                        toupper(substr(quality_band, 1, 1)), "]")]

# Tri-color accuracy heatmap
long <- rbindlist(list(
  scenario_means[, .(scenario_id, prevalence, profile_name, quality_band, N, regime,
                     metric = "kappa", correct = kappa_correct)],
  scenario_means[, .(scenario_id, prevalence, profile_name, quality_band, N, regime,
                     metric = "PABAK", correct = PABAK_correct)],
  scenario_means[, .(scenario_id, prevalence, profile_name, quality_band, N, regime,
                     metric = "AC1", correct = AC1_correct)]
))

agg <- long[, .(prop_correct = mean(correct)),
            by = .(prevalence, profile_name, quality_band, metric)]
agg[, display_name := paste0(clean_name(profile_name), "  [",
                              toupper(substr(quality_band, 1, 1)), "]")]
agg[, display_name := factor(display_name, levels = profile_order$display_name)]
agg[, metric := factor(metric, levels = c("kappa", "PABAK", "AC1"))]

band_breaks <- c(
  sum(profile_order$quality_band == "high") + 0.5,
  sum(profile_order$quality_band %in% c("high", "medium")) + 0.5
)

fig_accuracy <- ggplot(agg, aes(
    x = as.numeric(factor(prevalence)) + (as.numeric(metric) - 2) * 0.28,
    y = display_name,
    fill = prop_correct)) +
  geom_tile(width = 0.26, height = 0.85, color = "grey85", linewidth = 0.15) +
  geom_hline(yintercept = band_breaks, color = "grey30", linewidth = 0.6, linetype = "dashed") +
  scale_fill_gradient2(
    low = "#C1272D", mid = "#FFF3B0", high = "#006837",
    midpoint = 0.5, limits = c(0, 1),
    name = "Classification\naccuracy"
  ) +
  scale_x_continuous(
    breaks = seq_along(sort(unique(agg$prevalence))),
    labels = sort(unique(agg$prevalence)),
    name = "Prevalence"
  ) +
  scale_y_discrete(name = NULL) +
  labs(
    title = "Metric classification accuracy across the unified simulation grid",
    subtitle = expression(paste(
      "Each cell: ", kappa, " | PABAK | AC1 (left to right). Aggregated across sample sizes. ",
      "[H] = High, [M] = Medium, [L] = Low quality."))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6, family = "mono"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

# Best metric heatmap
best <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy = mean(AC1_correct)
), by = .(prevalence, profile_name, quality_band)]

best[, best_metric := ifelse(
  kappa_accuracy >= PABAK_accuracy & kappa_accuracy >= AC1_accuracy, "kappa",
  ifelse(PABAK_accuracy >= AC1_accuracy, "PABAK", "AC1")
)]
best[, all_equal := (kappa_accuracy == PABAK_accuracy) & (PABAK_accuracy == AC1_accuracy)]
best[all_equal == TRUE, best_metric := "All equal"]

best[, display_name := paste0(clean_name(profile_name), "  [",
                               toupper(substr(quality_band, 1, 1)), "]")]
best[, display_name := factor(display_name, levels = profile_order$display_name)]

fig_best <- ggplot(best, aes(x = factor(prevalence), y = display_name, fill = best_metric)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_hline(yintercept = band_breaks, color = "grey30", linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(
    values = c("kappa" = "#E41A1C", "PABAK" = "#377EB8",
               "AC1" = "#4DAF4A", "All equal" = "#D9D9D9"),
    name = "Best metric"
  ) +
  labs(
    x = "Prevalence", y = NULL,
    title = "Which metric best recovers ground-truth rater quality?",
    subtitle = "Grey = all three metrics perform equally. [H] = High, [M] = Medium, [L] = Low quality."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6, family = "mono"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

# Save
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/fig_unified_accuracy.png", fig_accuracy, width = 16, height = 14, dpi = 300)
ggsave("output/figures/fig_unified_accuracy.pdf", fig_accuracy, width = 16, height = 14)
ggsave("output/figures/fig_unified_best.png", fig_best, width = 14, height = 14, dpi = 300)
ggsave("output/figures/fig_unified_best.pdf", fig_best, width = 14, height = 14)

# Save all results
saveRDS(list(
  scenario_means = scenario_means,
  thresholds = thresholds,
  regime_acc = regime_acc,
  regime_qual = regime_qual,
  discord = discord
), "output/unified_sim/unified_results.rds")

cat("\nAll outputs saved.\n")
