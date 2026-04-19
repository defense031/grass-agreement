#!/usr/bin/env Rscript
# run_thresholds.R — Derive context-aware thresholds from simulation results

source("R/00_packages.R")
source("R/06_threshold_derivation.R")

summaries <- readRDS("output/scenario_summaries.rds")
config    <- yaml::read_yaml("config.yaml")

dir.create("output/thresholds", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Strategy 1: Decision-utility thresholds
# ------------------------------------------------------------------
message("Strategy 1: Decision-utility mapping...")

# Test multiple loss thresholds
utility_results <- list()
for (lt in c(0.05, 0.10, 0.15, 0.20)) {
  res <- derive_utility_thresholds(summaries, config$cost_ratios, loss_threshold = lt)
  res$loss_threshold <- lt
  utility_results[[length(utility_results) + 1]] <- res
}
utility_all <- data.table::rbindlist(utility_results)
saveRDS(utility_all, "output/thresholds/utility_thresholds.rds")

message(sprintf("   %d threshold entries across %d loss thresholds",
                nrow(utility_all), length(unique(utility_all$loss_threshold))))

# Print summary for default settings (cr=1, loss<0.10, require rule)
default_util <- utility_all[cost_ratio == 1 & rule == "require" &
                              loss_threshold == 0.10 & metric == "PABAK", ]
message("\n   Default utility thresholds (cr=1, require-agreement, loss<0.10, PABAK):")
for (i in seq_len(nrow(default_util))) {
  r <- default_util[i, ]
  message(sprintf("     prev=%.2f: PABAK >= %.3f  (J=%.3f, Se=%.2f, Sp=%.2f, n_accept=%d/%d)",
                  r$prevalence, r$threshold, r$youden_j,
                  r$sensitivity, r$specificity,
                  r$n_acceptable, r$n_acceptable + r$n_unacceptable))
}

# ------------------------------------------------------------------
# Strategy 2: ROC thresholds for quality band classification
# ------------------------------------------------------------------
message("\nStrategy 2: ROC-style thresholds...")

roc_results <- derive_roc_thresholds(summaries)
saveRDS(roc_results, "output/thresholds/roc_thresholds.rds")

message("   Thresholds for 'high' quality classification by PABAK:")
roc_high_pabak <- roc_results[target == "high" & metric == "PABAK", ]
for (i in seq_len(nrow(roc_high_pabak))) {
  r <- roc_high_pabak[i, ]
  message(sprintf("     prev=%.2f: PABAK >= %.3f  (J=%.3f, Se=%.2f, Sp=%.2f)",
                  r$prevalence, r$threshold, r$youden_j,
                  r$sensitivity, r$specificity))
}

message("\n   Thresholds for 'high' quality classification by kappa:")
roc_high_kappa <- roc_results[target == "high" & metric == "kappa", ]
for (i in seq_len(nrow(roc_high_kappa))) {
  r <- roc_high_kappa[i, ]
  message(sprintf("     prev=%.2f: kappa >= %.3f  (J=%.3f, Se=%.2f, Sp=%.2f)",
                  r$prevalence, r$threshold, r$youden_j,
                  r$sensitivity, r$specificity))
}

# ------------------------------------------------------------------
# Strategy 3: Empirical clustering
# ------------------------------------------------------------------
message("\nStrategy 3: Empirical clustering...")

cluster_km <- derive_clustering_thresholds(summaries, k = 3, method = "kmeans")
cluster_hc <- derive_clustering_thresholds(summaries, k = 3, method = "hclust")
cluster_all <- rbind(cluster_km, cluster_hc)
saveRDS(cluster_all, "output/thresholds/clustering_results.rds")

message("   Clustering performance (ARI / Purity) by prevalence:")
message("   kmeans:")
for (i in seq_len(nrow(cluster_km))) {
  r <- cluster_km[i, ]
  message(sprintf("     prev=%.2f: ARI=%.3f  Purity=%.3f  (n=%d)",
                  r$prevalence, r$ari, r$purity, r$n_scenarios))
}
message("   hclust (Ward.D2):")
for (i in seq_len(nrow(cluster_hc))) {
  r <- cluster_hc[i, ]
  message(sprintf("     prev=%.2f: ARI=%.3f  Purity=%.3f  (n=%d)",
                  r$prevalence, r$ari, r$purity, r$n_scenarios))
}

# ------------------------------------------------------------------
# Compare context-aware vs Landis-Koch
# ------------------------------------------------------------------
message("\n=== COMPARISON: Context-Aware vs Landis-Koch ===")

lk_label <- function(x) {
  cut(x, breaks = c(-Inf, 0, 0.20, 0.40, 0.60, 0.80, Inf),
      labels = c("poor", "slight", "fair", "moderate", "substantial", "almost_perfect"),
      right = TRUE)
}

summaries[, lk_kappa := lk_label(kappa_mean)]
summaries[, lk_PABAK := lk_label(PABAK_mean)]

# How often does L-K kappa label match ground-truth quality band?
# Map L-K to 3 bands: poor/slight/fair -> "low", moderate -> "medium", substantial/almost_perfect -> "high"
lk_to_band <- function(lk) {
  ifelse(lk %in% c("poor", "slight", "fair"), "low",
         ifelse(lk == "moderate", "medium", "high"))
}

summaries[, lk_kappa_band := lk_to_band(as.character(lk_kappa))]
summaries[, lk_PABAK_band := lk_to_band(as.character(lk_PABAK))]

kappa_accuracy <- mean(summaries$lk_kappa_band == summaries$quality_band, na.rm = TRUE)
pabak_accuracy <- mean(summaries$lk_PABAK_band == summaries$quality_band, na.rm = TRUE)

message(sprintf("   L-K(kappa) matches ground truth: %.1f%%", 100 * kappa_accuracy))
message(sprintf("   L-K(PABAK) matches ground truth: %.1f%%", 100 * pabak_accuracy))

# By prevalence
message("\n   L-K accuracy by prevalence:")
for (p in sort(unique(summaries$prevalence))) {
  sub <- summaries[prevalence == p, ]
  k_acc <- 100 * mean(sub$lk_kappa_band == sub$quality_band, na.rm = TRUE)
  p_acc <- 100 * mean(sub$lk_PABAK_band == sub$quality_band, na.rm = TRUE)
  message(sprintf("     prev=%.2f:  kappa=%.1f%%  PABAK=%.1f%%", p, k_acc, p_acc))
}

saveRDS(summaries, "output/scenario_summaries.rds")

message("\n=== Threshold derivation complete ===")
message("Saved: output/thresholds/utility_thresholds.rds")
message("Saved: output/thresholds/roc_thresholds.rds")
message("Saved: output/thresholds/clustering_results.rds")
message("Updated: output/scenario_summaries.rds")
