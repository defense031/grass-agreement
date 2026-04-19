#!/usr/bin/env Rscript
# run_best_comps.R — Publication-quality metric comparison figures
#
# Figure 1: Tri-color accuracy heatmap (aggregated across N)
# Figure 2: Best metric heatmap faceted by N

source("R/00_packages.R")
library(data.table)
library(ggplot2)

res <- readRDS("output/metric_recommendation.rds")
d <- res$scenario_data

# ------------------------------------------------------------------
# Common setup
# ------------------------------------------------------------------

# Clean profile names
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  # Title case
  x <- gsub("(^|\\s)(\\w)", "\\1\\U\\2", x, perl = TRUE)
  # Fix specific labels
  x <- gsub("Vlow", "Very Low", x)
  x <- gsub("Vhigh", "Very High", x)
  x <- gsub("Asym ", "Asym. ", x)
  x <- gsub("Med$", "Moderate", x)
  x
}

# Profile ordering: high at top, low at bottom, alphabetical within
profile_order <- unique(d[, .(profile_name, quality_band, min_accuracy)])
profile_order <- profile_order[order(
  -factor(quality_band, levels = c("high", "medium", "low")),
  -min_accuracy, profile_name
)]
profile_order[, display_name := paste0(clean_name(profile_name), "  [",
                                        toupper(substr(quality_band, 1, 1)), "]")]

# ------------------------------------------------------------------
# Figure 1: Tri-color accuracy heatmap
# ------------------------------------------------------------------

# Long format: one row per scenario x metric
long <- rbindlist(list(
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "kappa", correct = kappa_correct)],
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "PABAK", correct = PABAK_correct)],
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "AC1", correct = AC1_correct)]
))

# Aggregate across N
agg <- long[, .(prop_correct = mean(correct)),
            by = .(prevalence, profile_name, quality_band, metric)]

agg[, display_name := paste0(clean_name(profile_name), "  [",
                              toupper(substr(quality_band, 1, 1)), "]")]
agg[, display_name := factor(display_name, levels = profile_order$display_name)]
agg[, metric := factor(metric, levels = c("kappa", "PABAK", "AC1"))]

# Horizontal separator positions between quality bands
band_breaks <- c(
  sum(profile_order$quality_band == "high") + 0.5,
  sum(profile_order$quality_band %in% c("high", "medium")) + 0.5
)

fig1 <- ggplot(agg, aes(
    x = as.numeric(factor(prevalence)) + (as.numeric(metric) - 2) * 0.28,
    y = display_name,
    fill = prop_correct)) +
  geom_tile(width = 0.26, height = 0.85, color = "grey85", linewidth = 0.15) +
  # Band separators
  geom_hline(yintercept = band_breaks, color = "grey30", linewidth = 0.6, linetype = "dashed") +
  scale_fill_gradient2(
    low = "#C1272D", mid = "#FFF3B0", high = "#006837",
    midpoint = 0.5, limits = c(0, 1),
    name = "Classification\naccuracy",
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  scale_x_continuous(
    breaks = seq_along(unique(agg$prevalence)),
    labels = sort(unique(agg$prevalence)),
    name = "Prevalence"
  ) +
  scale_y_discrete(name = NULL) +
  # Metric labels at top
  annotate("text", x = 0.5, y = nrow(profile_order) + 1.2,
           label = expression(kappa), size = 3, hjust = 0.5) +
  annotate("text", x = 0.78, y = nrow(profile_order) + 1.2,
           label = "P", size = 3, hjust = 0.5) +
  annotate("text", x = 1.06, y = nrow(profile_order) + 1.2,
           label = "A", size = 3, hjust = 0.5) +
  labs(
    title = "Metric classification accuracy across prevalence and rater profiles",
    subtitle = expression(paste(
      "Each cell: ", kappa, " | PABAK | AC1 (left to right). Aggregated across sample sizes. ",
      "[H] = High, [M] = Medium, [L] = Low quality."))
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7.5, family = "mono"),
    axis.text.x = element_text(size = 9),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

ggsave("output/figures/best_comps/fig_accuracy_heatmap.png", fig1,
       width = 13, height = 10, dpi = 300)
ggsave("output/figures/best_comps/fig_accuracy_heatmap.pdf", fig1,
       width = 13, height = 10)

# ------------------------------------------------------------------
# Figure 2: Best metric heatmap, faceted by N
# ------------------------------------------------------------------

# Per prevalence x profile x N
best_by_n <- d[, .(
  kappa_accuracy = mean(kappa_correct),
  PABAK_accuracy = mean(PABAK_correct),
  AC1_accuracy   = mean(AC1_correct)
), by = .(prevalence, profile_name, quality_band, N)]

best_by_n[, best_metric := ifelse(
  kappa_accuracy >= PABAK_accuracy & kappa_accuracy >= AC1_accuracy, "kappa",
  ifelse(PABAK_accuracy >= AC1_accuracy, "PABAK", "AC1")
)]
best_by_n[, all_equal := (kappa_accuracy == PABAK_accuracy) &
                          (PABAK_accuracy == AC1_accuracy)]
best_by_n[all_equal == TRUE, best_metric := "All equal"]

best_by_n[, display_name := paste0(clean_name(profile_name), "  [",
                                    toupper(substr(quality_band, 1, 1)), "]")]
best_by_n[, display_name := factor(display_name, levels = profile_order$display_name)]

# Numeric facet ordering
best_by_n[, N_label := factor(paste0("N = ", N),
                               levels = paste0("N = ", c(10, 50, 100, 500, 1000)))]

fig2 <- ggplot(best_by_n, aes(x = factor(prevalence), y = display_name, fill = best_metric)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_hline(yintercept = band_breaks, color = "grey30", linewidth = 0.5, linetype = "dashed") +
  facet_wrap(~ N_label, nrow = 1) +
  scale_fill_manual(
    values = c("kappa" = "#E41A1C", "PABAK" = "#377EB8",
               "AC1" = "#4DAF4A", "All equal" = "#D9D9D9"),
    name = "Best metric"
  ) +
  labs(
    x = "Prevalence",
    y = NULL,
    title = "Which metric best recovers ground-truth rater quality?",
    subtitle = "Best = highest classification accuracy. Grey = all three metrics perform equally. [H] = High, [M] = Medium, [L] = Low quality."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6.5, family = "mono"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm")
  )

ggsave("output/figures/best_comps/fig_best_metric_by_n.png", fig2,
       width = 18, height = 10, dpi = 300)
ggsave("output/figures/best_comps/fig_best_metric_by_n.pdf", fig2,
       width = 18, height = 10)

cat("Saved to output/figures/best_comps/\n")
