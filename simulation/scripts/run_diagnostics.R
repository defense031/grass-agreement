#!/usr/bin/env Rscript
# run_diagnostics.R — MC diagnostics, convergence checks, result documentation
#
# Produces: output/diagnostics_report.txt
#           output/scenario_summaries.rds

source("R/00_packages.R")

results <- readRDS("output/sim_results/all_results.rds")
grid    <- readRDS("output/parameter_grid.rds")

sink("output/diagnostics_report.txt")

cat("================================================================\n")
cat("  PABAK SIMULATION — MC DIAGNOSTICS REPORT\n")
cat(sprintf("  Generated: %s\n", Sys.time()))
cat("================================================================\n\n")

# ------------------------------------------------------------------
# 1. Overall summary
# ------------------------------------------------------------------
cat("1. OVERALL SUMMARY\n")
cat(sprintf("   Total replications:  %s\n", format(nrow(results), big.mark = ",")))
cat(sprintf("   Total scenarios:     %d\n", nrow(grid)))
cat(sprintf("   Standard (5000):     %d\n", sum(grid$n_reps == 5000)))
cat(sprintf("   Stress   (20000):    %d\n", sum(grid$n_reps == 20000)))
cat(sprintf("   File size:           %.1f MB\n\n", file.size("output/sim_results/all_results.rds") / 1e6))

# ------------------------------------------------------------------
# 2. Scenario-level summaries
# ------------------------------------------------------------------
cat("2. COMPUTING SCENARIO-LEVEL SUMMARIES...\n")

summaries <- results[, .(
  n_reps        = .N,
  # P0
  P0_mean       = mean(P0),
  P0_sd         = sd(P0),
  P0_mcse       = sd(P0) / sqrt(.N),
  # Kappa
  kappa_mean    = mean(kappa, na.rm = TRUE),
  kappa_sd      = sd(kappa, na.rm = TRUE),
  kappa_mcse    = sd(kappa, na.rm = TRUE) / sqrt(sum(!is.na(kappa))),
  kappa_na_pct  = 100 * sum(is.na(kappa)) / .N,
  # PABAK
  PABAK_mean    = mean(PABAK),
  PABAK_sd      = sd(PABAK),
  PABAK_mcse    = sd(PABAK) / sqrt(.N),
  # Positive / negative agreement
  pos_agree_mean = mean(pos_agreement, na.rm = TRUE),
  pos_agree_na   = 100 * sum(is.na(pos_agreement)) / .N,
  neg_agree_mean = mean(neg_agreement, na.rm = TRUE),
  neg_agree_na   = 100 * sum(is.na(neg_agreement)) / .N,
  # Prevalence / bias indices
  prev_idx_mean  = mean(prevalence_index),
  bias_idx_mean  = mean(bias_index),
  # Divergence: |PABAK - kappa|
  divergence_mean = mean(abs(PABAK - kappa), na.rm = TRUE),
  divergence_max  = max(abs(PABAK - kappa), na.rm = TRUE)
), by = scenario_id]

# Merge grid info
summaries <- merge(summaries, grid, by = "scenario_id", suffixes = c("", ".grid"))

saveRDS(summaries, "output/scenario_summaries.rds")
cat(sprintf("   Saved: output/scenario_summaries.rds (%d rows)\n\n", nrow(summaries)))

# ------------------------------------------------------------------
# 3. MC standard error diagnostics
# ------------------------------------------------------------------
cat("3. MONTE CARLO STANDARD ERROR DIAGNOSTICS\n\n")

cat("   MCSE(kappa) distribution:\n")
q_mcse <- quantile(summaries$kappa_mcse, c(0, 0.25, 0.50, 0.75, 0.90, 0.95, 1.0), na.rm = TRUE)
for (nm in names(q_mcse)) {
  cat(sprintf("     %5s: %.6f\n", nm, q_mcse[nm]))
}

flagged <- summaries[kappa_mcse > 0.01, ]
cat(sprintf("\n   Scenarios with MCSE(kappa) > 0.01: %d / %d\n", nrow(flagged), nrow(summaries)))
if (nrow(flagged) > 0) {
  cat("   Flagged scenarios (top 10):\n")
  flagged_top <- head(flagged[order(-kappa_mcse)], 10)
  for (i in seq_len(nrow(flagged_top))) {
    r <- flagged_top[i, ]
    cat(sprintf("     scenario %3d: %-20s prev=%.2f N=%4d  MCSE=%.4f  kappa_mean=%.4f\n",
                r$scenario_id, r$profile_name, r$prevalence, r$N, r$kappa_mcse, r$kappa_mean))
  }
}

# ------------------------------------------------------------------
# 4. Convergence: MC means vs theoretical
# ------------------------------------------------------------------
cat("\n4. CONVERGENCE CHECK: MC MEANS vs THEORETICAL\n\n")

summaries[, P0_error := abs(P0_mean - P0_theoretical)]
summaries[, kappa_error := abs(kappa_mean - kappa_theoretical)]

cat("   |MC_mean(P0) - theoretical P0| distribution:\n")
q_p0 <- quantile(summaries$P0_error, c(0.50, 0.90, 0.95, 0.99, 1.0))
for (nm in names(q_p0)) cat(sprintf("     %5s: %.6f\n", nm, q_p0[nm]))

cat("\n   |MC_mean(kappa) - theoretical kappa| distribution:\n")
q_k <- quantile(summaries$kappa_error, c(0.50, 0.90, 0.95, 0.99, 1.0), na.rm = TRUE)
for (nm in names(q_k)) cat(sprintf("     %5s: %.6f\n", nm, q_k[nm]))

big_miss <- summaries[kappa_error > 0.02, ]
cat(sprintf("\n   Scenarios with |kappa_error| > 0.02: %d\n", nrow(big_miss)))
if (nrow(big_miss) > 0) {
  cat("   These are expected at extreme prevalence + small N:\n")
  for (i in seq_len(min(nrow(big_miss), 10))) {
    r <- big_miss[i, ]
    cat(sprintf("     scenario %3d: %-20s prev=%.2f N=%4d  MC=%.4f  theo=%.4f  err=%.4f\n",
                r$scenario_id, r$profile_name, r$prevalence, r$N,
                r$kappa_mean, r$kappa_theoretical, r$kappa_error))
  }
}

# ------------------------------------------------------------------
# 5. Degenerate cases
# ------------------------------------------------------------------
cat("\n5. DEGENERATE CASE FREQUENCY\n\n")

cat("   Kappa NA rate by scenario (top 10):\n")
high_na <- head(summaries[order(-kappa_na_pct)], 10)
for (i in seq_len(nrow(high_na))) {
  r <- high_na[i, ]
  cat(sprintf("     scenario %3d: %-20s prev=%.2f N=%4d  kappa_NA=%.1f%%\n",
              r$scenario_id, r$profile_name, r$prevalence, r$N, r$kappa_na_pct))
}

cat("\n   Positive agreement NA rate (top 10):\n")
high_pos_na <- head(summaries[order(-pos_agree_na)], 10)
for (i in seq_len(nrow(high_pos_na))) {
  r <- high_pos_na[i, ]
  cat(sprintf("     scenario %3d: %-20s prev=%.2f N=%4d  pos_agree_NA=%.1f%%\n",
              r$scenario_id, r$profile_name, r$prevalence, r$N, r$pos_agree_na))
}

# ------------------------------------------------------------------
# 6. Kappa vs PABAK divergence overview
# ------------------------------------------------------------------
cat("\n6. KAPPA vs PABAK DIVERGENCE OVERVIEW\n\n")

cat("   Mean |PABAK - kappa| distribution across scenarios:\n")
q_div <- quantile(summaries$divergence_mean, c(0, 0.25, 0.50, 0.75, 0.90, 0.95, 1.0), na.rm = TRUE)
for (nm in names(q_div)) cat(sprintf("     %5s: %.4f\n", nm, q_div[nm]))

cat("\n   Top 10 scenarios by mean divergence:\n")
top_div <- head(summaries[order(-divergence_mean)], 10)
for (i in seq_len(nrow(top_div))) {
  r <- top_div[i, ]
  cat(sprintf("     scenario %3d: %-20s prev=%.2f N=%4d  kappa=%.3f  PABAK=%.3f  div=%.3f\n",
              r$scenario_id, r$profile_name, r$prevalence, r$N,
              r$kappa_mean, r$PABAK_mean, r$divergence_mean))
}

# ------------------------------------------------------------------
# 7. Paradox frequency: high P0 but low kappa
# ------------------------------------------------------------------
cat("\n7. PARADOX FREQUENCY (P0 > 0.80 AND kappa < 0.40)\n\n")

paradox_reps <- results[P0 > 0.80 & kappa < 0.40 & !is.na(kappa), ]
cat(sprintf("   Replications meeting paradox criteria: %s / %s (%.2f%%)\n",
            format(nrow(paradox_reps), big.mark = ","),
            format(nrow(results), big.mark = ","),
            100 * nrow(paradox_reps) / nrow(results)))

# By prevalence band
paradox_by_prev <- results[!is.na(kappa), .(
  total = .N,
  paradox = sum(P0 > 0.80 & kappa < 0.40),
  pct = 100 * sum(P0 > 0.80 & kappa < 0.40) / .N
), by = .(prevalence = results[!is.na(kappa), ]$scenario_id)]

# Merge prevalence from grid
paradox_by_prev <- merge(paradox_by_prev, grid[, c("scenario_id", "prevalence")],
                          by.x = "prevalence", by.y = "scenario_id")
paradox_by_prev2 <- paradox_by_prev[, .(
  total = sum(total),
  paradox = sum(paradox),
  pct = 100 * sum(paradox) / sum(total)
), by = prevalence.y]
setorder(paradox_by_prev2, prevalence.y)

cat("\n   Paradox rate by prevalence level:\n")
for (i in seq_len(nrow(paradox_by_prev2))) {
  r <- paradox_by_prev2[i, ]
  cat(sprintf("     prev=%.2f: %6.2f%% (%s / %s)\n",
              r$prevalence.y, r$pct,
              format(r$paradox, big.mark = ","),
              format(r$total, big.mark = ",")))
}

# ------------------------------------------------------------------
# 8. Landis-Koch mislabeling analysis
# ------------------------------------------------------------------
cat("\n8. LANDIS-KOCH MISLABELING ANALYSIS\n\n")

# Apply Landis-Koch labels to kappa and PABAK
lk_label <- function(x) {
  cut(x, breaks = c(-Inf, 0, 0.20, 0.40, 0.60, 0.80, Inf),
      labels = c("poor", "slight", "fair", "moderate", "substantial", "almost_perfect"),
      right = TRUE)
}

# At scenario level using means
summaries[, lk_kappa := lk_label(kappa_mean)]
summaries[, lk_PABAK := lk_label(PABAK_mean)]
summaries[, lk_disagree := lk_kappa != lk_PABAK]

cat(sprintf("   Scenarios where L-K label differs for kappa vs PABAK: %d / %d (%.1f%%)\n",
            sum(summaries$lk_disagree, na.rm = TRUE), nrow(summaries),
            100 * mean(summaries$lk_disagree, na.rm = TRUE)))

cat("\n   Cross-tabulation of L-K labels (kappa vs PABAK):\n")
ct <- table(kappa = summaries$lk_kappa, PABAK = summaries$lk_PABAK)
print(ct)

# Most dramatic mislabels: kappa says "poor/slight" but PABAK says "substantial/almost_perfect"
dramatic <- summaries[lk_kappa %in% c("poor", "slight") & lk_PABAK %in% c("substantial", "almost_perfect"), ]
cat(sprintf("\n   Dramatic mislabels (kappa=poor/slight, PABAK=substantial+): %d scenarios\n", nrow(dramatic)))
if (nrow(dramatic) > 0) {
  for (i in seq_len(min(nrow(dramatic), 10))) {
    r <- dramatic[i, ]
    cat(sprintf("     scenario %3d: %-20s prev=%.2f N=%4d  kappa=%.3f(%s)  PABAK=%.3f(%s)\n",
                r$scenario_id, r$profile_name, r$prevalence, r$N,
                r$kappa_mean, r$lk_kappa, r$PABAK_mean, r$lk_PABAK))
  }
}

cat("\n================================================================\n")
cat("  END OF DIAGNOSTICS REPORT\n")
cat("================================================================\n")

sink()

cat("Diagnostics report written to: output/diagnostics_report.txt\n")
cat("Scenario summaries saved to:   output/scenario_summaries.rds\n")
