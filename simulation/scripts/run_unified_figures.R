#!/usr/bin/env Rscript
# run_unified_figures.R — Regenerate ALL figures from unified simulation
# Also derives Medium-quality thresholds and generates both table versions.

source("R/00_packages.R")
library(data.table)
library(ggplot2)

res <- readRDS("output/unified_sim/unified_results.rds")
d <- res$scenario_means
grid <- as.data.table(readRDS("output/unified_sim/parameter_grid.rds"))
grid[, regime := ifelse(PI_theoretical > BI_theoretical, "prevalence_dominant",
                 ifelse(BI_theoretical > PI_theoretical, "bias_dominant", "balanced"))]

dir.create("figures", showWarnings = FALSE)

# ==================================================================
# Figure 1: Three metrics across prevalence (Se=Sp=0.90, N=100)
# ==================================================================
cat("Generating fig_three_metrics...\n")

# Find the symmetric_high profile at N=100
three_met <- d[profile_name == "symmetric_high" & N == 100,
               .(prevalence, kappa_mean, PABAK_mean, AC1_mean)]
setorder(three_met, prevalence)

three_long <- melt(three_met, id.vars = "prevalence",
                   variable.name = "metric", value.name = "value")
three_long[, metric := gsub("_mean", "", metric)]
three_long[, metric := factor(metric, levels = c("kappa", "PABAK", "AC1"),
                               labels = c(expression(kappa), "PABAK", "AC1"))]

# Need character labels for ggplot
three_long2 <- melt(three_met, id.vars = "prevalence",
                    variable.name = "metric", value.name = "value")
three_long2[, metric := gsub("_mean", "", metric)]
three_long2[, metric := factor(metric, levels = c("kappa", "PABAK", "AC1"))]

fig_three <- ggplot(three_long2, aes(x = prevalence, y = value, color = metric, shape = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.61, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.15, y = 0.63, label = 'L-K "substantial" (0.61)',
           size = 3, color = "grey40") +
  scale_color_manual(values = c("kappa" = "#E41A1C", "PABAK" = "#377EB8", "AC1" = "#4DAF4A"),
                     labels = c(expression(kappa), "PABAK", "AC1")) +
  scale_shape_manual(values = c("kappa" = 16, "PABAK" = 17, "AC1" = 15),
                     labels = c(expression(kappa), "PABAK", "AC1")) +
  scale_x_continuous(breaks = sort(unique(three_met$prevalence))) +
  labs(x = "Prevalence", y = "Metric value",
       color = "Metric", shape = "Metric") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position.inside = c(0.5, 0.15),
    legend.background = element_rect(fill = "white", color = "grey80"),
    panel.grid.minor = element_blank()
  )

ggsave("figures/fig_three_metrics.png", fig_three, width = 8, height = 5.5, dpi = 300)
ggsave("figures/fig_three_metrics.pdf", fig_three, width = 8, height = 5.5)

# ==================================================================
# Figure 2: Tri-color accuracy heatmap (from unified analysis)
# ==================================================================
cat("Generating fig_accuracy_heatmap...\n")

clean_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("(^|\\s)(\\w)", "\\1\\U\\2", x, perl = TRUE)
  x
}

profile_order <- unique(d[, .(profile_name, quality_band, min_accuracy)])
profile_order <- profile_order[order(
  -factor(quality_band, levels = c("high", "medium", "low")),
  -min_accuracy, profile_name
)]
profile_order[, display_name := paste0(clean_name(profile_name), "  [",
                                        toupper(substr(quality_band, 1, 1)), "]")]

long <- rbindlist(list(
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "kappa", correct = kappa_correct)],
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
        metric = "PABAK", correct = PABAK_correct)],
  d[, .(scenario_id, prevalence, profile_name, quality_band, N,
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

ggsave("figures/fig_accuracy_heatmap.png", fig_accuracy, width = 16, height = 14, dpi = 300)
ggsave("figures/fig_accuracy_heatmap.pdf", fig_accuracy, width = 16, height = 14)

# ==================================================================
# Figure 3: Best metric heatmap
# ==================================================================
cat("Generating fig_best_metric...\n")

best <- d[, .(
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

ggsave("figures/fig_best_metric.png", fig_best, width = 14, height = 14, dpi = 300)
ggsave("figures/fig_best_metric.pdf", fig_best, width = 14, height = 14)

# ==================================================================
# Derive Medium-quality thresholds
# ==================================================================
cat("Deriving Medium thresholds...\n")

derive_threshold <- function(metric_vals, labels) {
  if (length(unique(labels)) < 2) return(list(threshold = NA, J = NA))
  candidates <- sort(unique(metric_vals))
  best_j <- -Inf; best_t <- NA
  for (t in candidates) {
    pred <- metric_vals >= t
    se <- sum(pred & labels) / max(sum(labels), 1)
    sp <- sum(!pred & !labels) / max(sum(!labels), 1)
    j <- se + sp - 1
    if (j > best_j) { best_j <- j; best_t <- t }
  }
  list(threshold = best_t, J = best_j)
}

# Medium = min(Se,Sp) >= 0.70 vs < 0.70
d[, is_medium_or_high := quality_band %in% c("high", "medium")]

prevalences <- sort(unique(d$prevalence))
med_thresh_list <- list()
for (p in prevalences) {
  sub <- d[prevalence == p]
  k <- derive_threshold(sub$kappa_mean, sub$is_medium_or_high)
  pb <- derive_threshold(sub$PABAK_mean, sub$is_medium_or_high)
  ac <- derive_threshold(sub$AC1_mean, sub$is_medium_or_high)
  med_thresh_list[[as.character(p)]] <- data.table(
    prevalence = p,
    kappa_med = k$threshold, J_kappa_med = k$J,
    PABAK_med = pb$threshold, J_PABAK_med = pb$J,
    AC1_med = ac$threshold, J_AC1_med = ac$J
  )
}
med_thresholds <- rbindlist(med_thresh_list)

cat("\nMedium thresholds:\n")
print(med_thresholds, digits = 3)

# Combine with High thresholds
high_thresholds <- res$thresholds
combined <- merge(high_thresholds, med_thresholds, by = "prevalence")

cat("\nCombined High + Medium thresholds:\n")
print(combined, digits = 3)

saveRDS(combined, "output/unified_sim/combined_thresholds.rds")

cat("\nAll figures and thresholds generated.\n")
