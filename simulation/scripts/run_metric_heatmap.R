#!/usr/bin/env Rscript
# run_metric_heatmap.R — Tri-color heatmap showing which metrics correctly
# classify rater quality at each prevalence x profile combination
#
# Each cell is split into three vertical bars (kappa, PABAK, AC1)
# colored by whether that metric correctly classified the scenario.

source("R/00_packages.R")
library(data.table)
library(ggplot2)

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
res <- readRDS("output/metric_recommendation.rds")
d <- res$scenario_data

# ------------------------------------------------------------------
# 2. Build long-format data for tri-color cells
# ------------------------------------------------------------------
# Each scenario gets three rows (one per metric)
long <- rbindlist(list(
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "kappa", correct = kappa_correct)],
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "PABAK", correct = PABAK_correct)],
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "AC1", correct = AC1_correct)]
))

# Aggregate across sample sizes: proportion correct per prevalence x profile x metric
agg <- long[, .(
  prop_correct = mean(correct),
  n = .N
), by = .(prevalence, profile_name, quality_band, metric)]

# ------------------------------------------------------------------
# 3. Order profiles by quality band
# ------------------------------------------------------------------
profile_order <- unique(d[, .(profile_name, quality_band, min_accuracy)])
profile_order <- profile_order[order(
  factor(quality_band, levels = c("low", "medium", "high")),
  min_accuracy, profile_name
)]

# Clean up profile names for display
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)
  x
}

profile_order[, display_name := paste0(clean_name(profile_name), "  [",
                                        toupper(substr(quality_band, 1, 1)), "]")]

agg[, display_name := paste0(clean_name(profile_name), "  [",
                              toupper(substr(quality_band, 1, 1)), "]")]
agg[, display_name := factor(display_name, levels = profile_order$display_name)]

# Metric ordering: kappa, PABAK, AC1 (left to right in each cell)
agg[, metric := factor(metric, levels = c("kappa", "PABAK", "AC1"))]

# ------------------------------------------------------------------
# 4. Plot: tri-color heatmap
# ------------------------------------------------------------------
# Use position_dodge to create side-by-side bars within each cell

p <- ggplot(agg, aes(x = factor(prevalence), y = display_name, fill = prop_correct)) +
  geom_tile(aes(x = as.numeric(factor(prevalence)) +
                  (as.numeric(metric) - 2) * 0.28),
            width = 0.27, height = 0.9, color = "grey90", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#D7191C", mid = "#FFFFBF", high = "#1A9641",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "Classification\naccuracy"
  ) +
  scale_x_discrete(name = "Prevalence") +
  scale_y_discrete(name = NULL) +
  labs(
    title = "Metric classification accuracy by prevalence and rater profile",
    subtitle = expression(
      paste("Each cell: ", kappa, " | PABAK | AC1 (left to right). ",
            "Green = correct, Red = incorrect. Aggregated across sample sizes.")
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7, family = "mono"),
    axis.text.x = element_text(size = 9),
    plot.subtitle = element_text(size = 9),
    legend.position = "right"
  )

# ------------------------------------------------------------------
# 5. Alternative: faceted by N to show sample size effect
# ------------------------------------------------------------------
agg_by_n <- long[, .(
  prop_correct = mean(correct)
), by = .(prevalence, profile_name, quality_band, metric, N)]

agg_by_n[, display_name := paste0(clean_name(profile_name), "  [",
                                    toupper(substr(quality_band, 1, 1)), "]")]
agg_by_n[, display_name := factor(display_name, levels = profile_order$display_name)]
agg_by_n[, metric := factor(metric, levels = c("kappa", "PABAK", "AC1"))]

# Simpler version: show only discordant scenarios to reduce visual noise
# A scenario is discordant if not all metrics have the same correctness
discord_check <- agg_by_n[, .(all_same = length(unique(prop_correct)) == 1),
                           by = .(prevalence, profile_name, N)]
discord_scenarios <- discord_check[all_same == FALSE]

agg_discord <- merge(agg_by_n, discord_scenarios[, .(prevalence, profile_name, N)],
                     by = c("prevalence", "profile_name", "N"))

p_discord <- ggplot(agg_discord, aes(x = factor(prevalence), y = display_name)) +
  geom_tile(aes(x = as.numeric(factor(prevalence)) +
                  (as.numeric(metric) - 2) * 0.28,
                fill = factor(prop_correct)),
            width = 0.27, height = 0.9, color = "grey80", linewidth = 0.3) +
  scale_fill_manual(
    values = c("0" = "#D7191C", "1" = "#1A9641"),
    labels = c("Incorrect", "Correct"),
    name = "Classification"
  ) +
  facet_wrap(~ paste0("N = ", N), nrow = 1) +
  scale_x_discrete(name = "Prevalence") +
  scale_y_discrete(name = NULL) +
  labs(
    title = "Discordant scenarios: which metric gets it right?",
    subtitle = expression(
      paste("Each cell: ", kappa, " | PABAK | AC1 (left to right). ",
            "Only scenarios where metrics disagree are shown.")
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7, family = "mono"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ------------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------------
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("output/figures/fig_metric_heatmap_all.png", p,
       width = 14, height = 10, dpi = 300)
ggsave("output/figures/fig_metric_heatmap_all.pdf", p,
       width = 14, height = 10)

ggsave("output/figures/fig_metric_heatmap_discord.png", p_discord,
       width = 16, height = 8, dpi = 300)
ggsave("output/figures/fig_metric_heatmap_discord.pdf", p_discord,
       width = 16, height = 8)

cat("Saved: output/figures/fig_metric_heatmap_all.{png,pdf}\n")
cat("Saved: output/figures/fig_metric_heatmap_discord.{png,pdf}\n")
