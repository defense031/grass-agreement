#!/usr/bin/env Rscript
# run_visualizations.R — Publication-quality figures for PABAK study

source("R/00_packages.R")

summaries   <- readRDS("output/scenario_summaries.rds")
utility_th  <- readRDS("output/thresholds/utility_thresholds.rds")
roc_th      <- readRDS("output/thresholds/roc_thresholds.rds")
cluster_res <- readRDS("output/thresholds/clustering_results.rds")

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Consistent theme
theme_pabak <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------
# Figure 1: Kappa vs PABAK divergence heatmap
# ------------------------------------------------------------------
message("Figure 1: Divergence heatmap...")

# Average divergence across sample sizes for cleaner view
div_summary <- summaries[, .(
  divergence = mean(divergence_mean),
  kappa = mean(kappa_mean, na.rm = TRUE),
  PABAK = mean(PABAK_mean)
), by = .(profile_name, prevalence)]

# Order profiles by mean divergence
prof_order <- div_summary[, .(mean_div = mean(divergence)), by = profile_name]
setorder(prof_order, mean_div)
div_summary[, profile_name := factor(profile_name, levels = prof_order$profile_name)]

p1 <- ggplot(div_summary, aes(x = factor(prevalence), y = profile_name, fill = divergence)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", divergence)), size = 3, color = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0.3, name = "Mean |PABAK - kappa|") +
  labs(title = "Kappa vs PABAK Divergence by Profile and Prevalence",
       x = "Prevalence", y = "Rater Profile") +
  theme_pabak +
  theme(legend.position = "right")

ggsave("output/figures/fig1_divergence_heatmap.pdf", p1, width = 10, height = 7)
ggsave("output/figures/fig1_divergence_heatmap.png", p1, width = 10, height = 7, dpi = 300)

# ------------------------------------------------------------------
# Figure 2: Paradox frequency by prevalence
# ------------------------------------------------------------------
message("Figure 2: Paradox frequency...")

results <- readRDS("output/sim_results/all_results.rds")
grid    <- readRDS("output/parameter_grid.rds")

# Merge prevalence into results via scenario_id
results_with_prev <- merge(
  results[, .(scenario_id, P0, kappa)],
  grid[, c("scenario_id", "prevalence", "N")],
  by = "scenario_id"
)

paradox_rates <- results_with_prev[!is.na(kappa), .(
  paradox_rate = 100 * mean(P0 > 0.80 & kappa < 0.40)
), by = .(prevalence, N)]

p2 <- ggplot(paradox_rates, aes(x = prevalence, y = paradox_rate,
                                  color = factor(N), group = factor(N))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Set1", name = "Sample Size (N)") +
  scale_x_continuous(breaks = unique(paradox_rates$prevalence)) +
  labs(title = "Paradox Frequency: High Agreement (P0 > 0.80) but Low Kappa (< 0.40)",
       x = "Prevalence", y = "% of Replications") +
  theme_pabak

ggsave("output/figures/fig2_paradox_frequency.pdf", p2, width = 9, height = 5.5)
ggsave("output/figures/fig2_paradox_frequency.png", p2, width = 9, height = 5.5, dpi = 300)

rm(results, results_with_prev)  # free memory

# ------------------------------------------------------------------
# Figure 3: Context-aware kappa thresholds vs prevalence
# ------------------------------------------------------------------
message("Figure 3: Kappa thresholds vs prevalence...")

roc_kappa_high <- roc_th[target == "high" & metric == "kappa", ]

p3 <- ggplot(roc_kappa_high, aes(x = prevalence, y = threshold)) +
  geom_line(linewidth = 1.2, color = "#d73027") +
  geom_point(size = 3, color = "#d73027") +
  geom_hline(yintercept = 0.61, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_hline(yintercept = 0.41, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  annotate("text", x = 0.85, y = 0.63, label = 'L-K "substantial"',
           size = 3.5, color = "grey40") +
  annotate("text", x = 0.85, y = 0.43, label = 'L-K "moderate"',
           size = 3.5, color = "grey60") +
  scale_x_continuous(breaks = unique(roc_kappa_high$prevalence)) +
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.7, 0.1)) +
  labs(title = "Context-Aware Kappa Threshold for 'High' Quality\nvs Fixed Landis-Koch Cut Points",
       x = "Prevalence", y = "Kappa Threshold") +
  theme_pabak

ggsave("output/figures/fig3_kappa_thresholds.pdf", p3, width = 8, height = 5.5)
ggsave("output/figures/fig3_kappa_thresholds.png", p3, width = 8, height = 5.5, dpi = 300)

# ------------------------------------------------------------------
# Figure 4: L-K accuracy by prevalence (kappa vs PABAK)
# ------------------------------------------------------------------
message("Figure 4: Landis-Koch accuracy...")

lk_acc <- data.frame(
  prevalence = rep(sort(unique(summaries$prevalence)), 2),
  metric = rep(c("kappa", "PABAK"), each = 9),
  accuracy = NA_real_
)

for (i in seq_len(nrow(lk_acc))) {
  p <- lk_acc$prevalence[i]
  sub <- summaries[prevalence == p, ]
  if (lk_acc$metric[i] == "kappa") {
    lk_acc$accuracy[i] <- 100 * mean(sub$lk_kappa_band == sub$quality_band, na.rm = TRUE)
  } else {
    lk_acc$accuracy[i] <- 100 * mean(sub$lk_PABAK_band == sub$quality_band, na.rm = TRUE)
  }
}

p4 <- ggplot(lk_acc, aes(x = prevalence, y = accuracy, color = metric, group = metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 33.3, linetype = "dotted", color = "grey50") +
  annotate("text", x = 0.75, y = 36, label = "Random chance (3 classes)",
           size = 3, color = "grey50") +
  scale_color_manual(values = c("kappa" = "#d73027", "PABAK" = "#4575b4"),
                     name = "Metric") +
  scale_x_continuous(breaks = unique(lk_acc$prevalence)) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
  labs(title = "Landis-Koch Label Accuracy vs Ground Truth by Prevalence",
       subtitle = "L-K labels mapped to quality bands (low/medium/high) and compared to Se/Sp-defined truth",
       x = "Prevalence", y = "Classification Accuracy (%)") +
  theme_pabak

ggsave("output/figures/fig4_lk_accuracy.pdf", p4, width = 9, height = 5.5)
ggsave("output/figures/fig4_lk_accuracy.png", p4, width = 9, height = 5.5, dpi = 300)

# ------------------------------------------------------------------
# Figure 5: Utility thresholds by cost ratio
# ------------------------------------------------------------------
message("Figure 5: Utility threshold surfaces...")

util_pabak <- utility_th[metric == "PABAK" & rule == "require" &
                            loss_threshold == 0.10, ]

# Replace NA thresholds with 1.0 (infeasible: no PABAK can achieve acceptable loss)
util_pabak[, threshold_plot := fifelse(is.na(threshold), 1.0, threshold)]
util_pabak[, feasible := !is.na(threshold)]

p5 <- ggplot(util_pabak, aes(x = prevalence, y = threshold_plot,
                               color = factor(cost_ratio), group = factor(cost_ratio))) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(shape = feasible), size = 2.5) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     labels = c("TRUE" = "Feasible", "FALSE" = "Infeasible"),
                     name = NULL) +
  scale_color_brewer(palette = "Dark2", name = expression(C[FN] / C[FP])) +
  scale_x_continuous(breaks = sort(unique(util_pabak$prevalence))) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "grey60") +
  annotate("text", x = 0.50, y = 1.03, label = "Infeasible region (no acceptable PABAK)",
           size = 3, color = "grey50") +
  labs(title = "Context-Aware PABAK Thresholds by Cost Ratio",
       subtitle = "Minimum PABAK for acceptable expected loss (<0.10), require-agreement rule",
       x = "Prevalence", y = "PABAK Threshold") +
  theme_pabak

ggsave("output/figures/fig5_utility_thresholds.pdf", p5, width = 9, height = 5.5)
ggsave("output/figures/fig5_utility_thresholds.png", p5, width = 9, height = 5.5, dpi = 300)

# ------------------------------------------------------------------
# Figure 6: Composite panel (key findings)
# ------------------------------------------------------------------
message("Figure 6: Composite panel...")

composite <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Context-Aware Interpretation Framework for PABAK and Kappa",
    subtitle = "Monte Carlo simulation: 3.1M replications across 495 scenarios",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, color = "grey30")
    )
  )

ggsave("output/figures/fig6_composite.pdf", composite, width = 16, height = 12)
ggsave("output/figures/fig6_composite.png", composite, width = 16, height = 12, dpi = 300)

message("\n=== All figures saved to output/figures/ ===")
message(sprintf("   %d PDF + %d PNG files",
                length(list.files("output/figures", pattern = "\\.pdf$")),
                length(list.files("output/figures", pattern = "\\.png$"))))
