#!/usr/bin/env Rscript
# run_balanced_analysis.R — Analyze the balanced bias stress test results
#
# Applies GRASS thresholds and compares performance across
# prevalence-dominant vs bias-dominant regimes.
#
# Produces: output/balanced_sim/balanced_results.rds
#           output/figures/fig_balanced_accuracy.png
#           output/figures/fig_balanced_by_regime.png

source("R/00_packages.R")
library(data.table)
library(ggplot2)

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
cat("Loading balanced simulation results...\n")
results <- readRDS("output/balanced_sim/sim_results/all_results.rds")
grid    <- as.data.table(readRDS("output/balanced_sim/parameter_grid.rds"))

# ------------------------------------------------------------------
# 2. Aggregate scenario-level means
# ------------------------------------------------------------------
cat("Aggregating scenario-level means...\n")
scenario_means <- results[, .(
  kappa_mean = mean(kappa, na.rm = TRUE),
  PABAK_mean = mean(PABAK),
  AC1_mean   = mean(AC1, na.rm = TRUE),
  P0_mean    = mean(P0),
  pos_agree_mean = mean(pos_agreement, na.rm = TRUE),
  neg_agree_mean = mean(neg_agreement, na.rm = TRUE),
  prev_idx_mean  = mean(prevalence_index),
  bias_idx_mean  = mean(bias_index),
  kappa_mcse = sd(kappa, na.rm = TRUE) / sqrt(sum(!is.na(kappa)))
), by = scenario_id]

scenario_means <- merge(scenario_means, grid, by = "scenario_id")

# ------------------------------------------------------------------
# 3. Apply GRASS thresholds (from Table 1 in paper)
# ------------------------------------------------------------------
# For bias-dominant scenarios at prevalences not in the original table,
# use nearest tabulated value
thresholds <- data.table(
  prevalence = c(0.01, 0.05, 0.10, 0.20, 0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60, 0.80, 0.90, 0.95, 0.99),
  kappa_thresh  = c(0.060, 0.240, 0.385, 0.492, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612, 0.482, 0.376, 0.241, 0.060),
  PABAK_thresh  = c(0.639, 0.639, 0.636, 0.638, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.638, 0.638, 0.639, 0.640),
  AC1_thresh    = c(0.759, 0.759, 0.740, 0.704, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.703, 0.730, 0.759, 0.759)
)

scenario_means <- merge(scenario_means, thresholds, by = "prevalence")

# ------------------------------------------------------------------
# 4. Classify and evaluate
# ------------------------------------------------------------------
scenario_means[, ground_truth := quality_band == "high"]
scenario_means[, kappa_classifies_high := kappa_mean >= kappa_thresh]
scenario_means[, PABAK_classifies_high := PABAK_mean >= PABAK_thresh]
scenario_means[, AC1_classifies_high := AC1_mean >= AC1_thresh]

scenario_means[, kappa_correct := kappa_classifies_high == ground_truth]
scenario_means[, PABAK_correct := PABAK_classifies_high == ground_truth]
scenario_means[, AC1_correct := AC1_classifies_high == ground_truth]

# ------------------------------------------------------------------
# 5. Summary by regime
# ------------------------------------------------------------------
cat("\n=== GRASS Threshold Performance by Regime ===\n\n")
regime_summary <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct),
  n_scenarios    = .N
), by = regime]
print(regime_summary)

cat("\n=== By Regime and Quality Band ===\n\n")
regime_quality <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct),
  n = .N
), by = .(regime, quality_band)]
print(regime_quality[order(regime, quality_band)])

cat("\n=== By Regime and Prevalence ===\n\n")
regime_prev <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct),
  n = .N
), by = .(regime, prevalence)]
print(regime_prev[order(regime, prevalence)])

# ------------------------------------------------------------------
# 6. Discordance analysis by regime
# ------------------------------------------------------------------
scenario_means[, all_agree := (kappa_classifies_high == PABAK_classifies_high) &
                               (PABAK_classifies_high == AC1_classifies_high)]

cat("\n=== Metric Agreement by Regime ===\n\n")
discord_summary <- scenario_means[, .(
  n_agree = sum(all_agree),
  n_discord = sum(!all_agree),
  pct_agree = 100 * mean(all_agree)
), by = regime]
print(discord_summary)

# ------------------------------------------------------------------
# 7. Clean profile names for plotting
# ------------------------------------------------------------------
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("(^|\\s)(\\w)", "\\1\\U\\2", x, perl = TRUE)
  x <- gsub("^Pd ", "PD-", x)
  x <- gsub("^Bd ", "BD-", x)
  x
}

# Profile ordering
profile_order <- unique(scenario_means[, .(profile_name, quality_band, min_accuracy, regime)])
profile_order <- profile_order[order(
  factor(regime, levels = c("prevalence_dominant", "bias_dominant")),
  -factor(quality_band, levels = c("high", "medium", "low")),
  -min_accuracy, profile_name
)]
profile_order[, display_name := paste0(clean_name(profile_name), "  [",
                                        toupper(substr(quality_band, 1, 1)), "]")]

# ------------------------------------------------------------------
# 8. Tri-color accuracy heatmap by regime
# ------------------------------------------------------------------
long <- rbindlist(list(
  scenario_means[, .(scenario_id, prevalence, profile_name, quality_band, N, regime,
                     metric = "kappa", correct = kappa_correct)],
  scenario_means[, .(scenario_id, prevalence, profile_name, quality_band, N, regime,
                     metric = "PABAK", correct = PABAK_correct)],
  scenario_means[, .(scenario_id, prevalence, profile_name, quality_band, N, regime,
                     metric = "AC1", correct = AC1_correct)]
))

# Aggregate across N
agg <- long[, .(prop_correct = mean(correct)),
            by = .(prevalence, profile_name, quality_band, metric, regime)]

agg[, display_name := paste0(clean_name(profile_name), "  [",
                              toupper(substr(quality_band, 1, 1)), "]")]
agg[, display_name := factor(display_name, levels = profile_order$display_name)]
agg[, metric := factor(metric, levels = c("kappa", "PABAK", "AC1"))]
agg[, regime_label := ifelse(regime == "prevalence_dominant",
                              "Prevalence-Dominant", "Bias-Dominant")]

fig1 <- ggplot(agg, aes(
    x = as.numeric(factor(prevalence)) + (as.numeric(metric) - 2) * 0.28,
    y = display_name,
    fill = prop_correct)) +
  geom_tile(width = 0.26, height = 0.85, color = "grey85", linewidth = 0.15) +
  facet_wrap(~ regime_label, scales = "free", ncol = 2) +
  scale_fill_gradient2(
    low = "#C1272D", mid = "#FFF3B0", high = "#006837",
    midpoint = 0.5, limits = c(0, 1),
    name = "Classification\naccuracy"
  ) +
  scale_x_continuous(
    breaks = seq_along(unique(agg$prevalence)),
    labels = function(x) {
      prevs <- sort(unique(agg$prevalence))
      ifelse(x <= length(prevs), as.character(prevs[x]), "")
    },
    name = "Prevalence"
  ) +
  scale_y_discrete(name = NULL) +
  labs(
    title = "GRASS threshold accuracy: prevalence-dominant vs. bias-dominant regimes",
    subtitle = expression(paste(
      "Each cell: ", kappa, " | PABAK | AC1 (left to right). Aggregated across sample sizes."))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6.5, family = "mono"),
    axis.text.x = element_text(size = 8),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    plot.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

# ------------------------------------------------------------------
# 9. Best metric by regime (simple summary)
# ------------------------------------------------------------------
best_by_regime <- scenario_means[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct)
), by = .(prevalence, profile_name, quality_band, regime)]

best_by_regime[, best_metric := ifelse(
  kappa_accuracy >= PABAK_accuracy & kappa_accuracy >= AC1_accuracy, "kappa",
  ifelse(PABAK_accuracy >= AC1_accuracy, "PABAK", "AC1")
)]
best_by_regime[, all_equal := (kappa_accuracy == PABAK_accuracy) &
                               (PABAK_accuracy == AC1_accuracy)]
best_by_regime[all_equal == TRUE, best_metric := "All equal"]

best_by_regime[, display_name := paste0(clean_name(profile_name), "  [",
                                         toupper(substr(quality_band, 1, 1)), "]")]
best_by_regime[, display_name := factor(display_name, levels = profile_order$display_name)]
best_by_regime[, regime_label := ifelse(regime == "prevalence_dominant",
                                         "Prevalence-Dominant", "Bias-Dominant")]

fig2 <- ggplot(best_by_regime, aes(x = factor(prevalence), y = display_name, fill = best_metric)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(~ regime_label, scales = "free", ncol = 2) +
  scale_fill_manual(
    values = c("kappa" = "#E41A1C", "PABAK" = "#377EB8",
               "AC1" = "#4DAF4A", "All equal" = "#D9D9D9"),
    name = "Best metric"
  ) +
  labs(
    x = "Prevalence", y = NULL,
    title = "Which metric best recovers ground truth, by regime?",
    subtitle = "Grey = all three metrics perform equally."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6.5, family = "mono"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

# ------------------------------------------------------------------
# 10. Save
# ------------------------------------------------------------------
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("output/figures/fig_balanced_accuracy.png", fig1,
       width = 18, height = 12, dpi = 300)
ggsave("output/figures/fig_balanced_accuracy.pdf", fig1,
       width = 18, height = 12)

ggsave("output/figures/fig_balanced_by_regime.png", fig2,
       width = 16, height = 10, dpi = 300)
ggsave("output/figures/fig_balanced_by_regime.pdf", fig2,
       width = 16, height = 10)

saveRDS(list(
  scenario_means = scenario_means,
  regime_summary = regime_summary,
  regime_quality = regime_quality,
  regime_prev    = regime_prev,
  discord_summary = discord_summary,
  best_by_regime = best_by_regime
), "output/balanced_sim/balanced_results.rds")

cat("\nFigures saved to output/figures/fig_balanced_*.{png,pdf}\n")
cat("Results saved to output/balanced_sim/balanced_results.rds\n")
