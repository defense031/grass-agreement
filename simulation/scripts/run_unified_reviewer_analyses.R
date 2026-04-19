#!/usr/bin/env Rscript
# run_unified_reviewer_analyses.R — Jackknife stability + CI coverage + sensitivity
# on the unified 2,850-scenario grid
#
# Produces: output/unified_sim/jackknife_stability.rds
#           output/unified_sim/jackknife_detail.rds
#           output/unified_sim/ci_coverage.rds
#           output/unified_sim/sensitivity_analysis.rds

library(data.table)
library(yaml)
source("R/00_packages.R")
source("R/03_ground_truth.R")
source("R/04_grid.R")

# Load unified data
grid <- as.data.table(readRDS("output/unified_sim/parameter_grid.rds"))
d <- readRDS("output/unified_sim/sim_results/all_results.rds")
setDT(d)

# Merge grid info
grid_cols <- grid[, .(scenario_id, prevalence, profile_name, quality_band,
                       kappa_theoretical, min_accuracy, grid_N = N)]
d2 <- merge(d, grid_cols, by = "scenario_id")

# =====================================================================
# PART 1: Jackknife cross-validation of kappa thresholds
# =====================================================================
cat("=== JACKKNIFE THRESHOLD STABILITY ===\n")

derive_kappa_threshold <- function(scenario_summaries) {
  prevs <- sort(unique(scenario_summaries$prevalence))
  results <- list()
  for (p in prevs) {
    sub <- scenario_summaries[prevalence == p]
    sub[, is_high := quality_band == "high"]
    if (length(unique(sub$is_high)) < 2) {
      results[[length(results) + 1]] <- data.table(prevalence = p, threshold = NA_real_, J = NA_real_)
      next
    }
    vals <- sort(unique(sub$kappa_mean))
    best_j <- -Inf; best_t <- NA_real_
    for (t in vals) {
      pred <- sub$kappa_mean >= t
      se <- sum(pred & sub$is_high) / max(sum(sub$is_high), 1)
      sp <- sum(!pred & !sub$is_high) / max(sum(!sub$is_high), 1)
      j <- se + sp - 1
      if (j > best_j) { best_j <- j; best_t <- t }
    }
    results[[length(results) + 1]] <- data.table(prevalence = p, threshold = best_t, J = round(best_j, 3))
  }
  rbindlist(results)
}

# Scenario-level summaries
summaries <- d2[, .(kappa_mean = mean(kappa, na.rm = TRUE),
                     PABAK_mean = mean(PABAK),
                     AC1_mean = mean(AC1, na.rm = TRUE)),
                 by = .(scenario_id, prevalence, profile_name, quality_band, grid_N)]

# Full-grid thresholds
full_thresholds <- derive_kappa_threshold(summaries)
setnames(full_thresholds, c("threshold", "J"), c("full_threshold", "full_J"))

# Leave-one-profile-out
profiles <- unique(summaries$profile_name)
cat(sprintf("Running jackknife across %d profiles...\n", length(profiles)))

jackknife_results <- list()
for (prof in profiles) {
  sub <- summaries[profile_name != prof]
  thresh <- derive_kappa_threshold(sub)
  thresh[, dropped_profile := prof]
  jackknife_results[[length(jackknife_results) + 1]] <- thresh
}
jk <- rbindlist(jackknife_results)

jk_summary <- jk[, .(min_thresh = min(threshold, na.rm = TRUE),
                       max_thresh = max(threshold, na.rm = TRUE),
                       mean_thresh = mean(threshold, na.rm = TRUE),
                       sd_thresh = sd(threshold, na.rm = TRUE)),
                   by = prevalence]
jk_summary <- merge(jk_summary, full_thresholds, by = "prevalence")
jk_summary[, max_deviation := pmax(abs(full_threshold - min_thresh), abs(full_threshold - max_thresh))]

cat("\nJackknife results:\n")
print(jk_summary[, .(prevalence, full_threshold = round(full_threshold, 3),
                       jk_min = round(min_thresh, 3), jk_max = round(max_thresh, 3),
                       max_dev = round(max_deviation, 3))])

dir.create("output/unified_sim", showWarnings = FALSE, recursive = TRUE)
saveRDS(jk_summary, "output/unified_sim/jackknife_stability.rds")
saveRDS(jk, "output/unified_sim/jackknife_detail.rds")

# =====================================================================
# PART 2: CI coverage proportions
# =====================================================================
cat("\n=== CI COVERAGE PROBABILITIES ===\n")

d2[, wald_covers := (!is.na(kappa_wald_lower)) &
     (kappa_wald_lower <= kappa_theoretical) &
     (kappa_wald_upper >= kappa_theoretical)]

d2[, logit_covers := (!is.na(kappa_wilson_lower)) &
     (kappa_wilson_lower <= kappa_theoretical) &
     (kappa_wilson_upper >= kappa_theoretical)]

# Wald boundary violations
d2[, wald_violates := (!is.na(kappa_wald_lower)) &
     (kappa_wald_lower < -1 | kappa_wald_upper > 1)]
d2[, logit_violates := (!is.na(kappa_wilson_lower)) &
     (kappa_wilson_lower < -1 | kappa_wilson_upper > 1)]

cat(sprintf("Wald boundary violations: %.1f%%\n", 100 * mean(d2$wald_violates, na.rm = TRUE)))
cat(sprintf("Logit boundary violations: %.1f%%\n", 100 * mean(d2$logit_violates, na.rm = TRUE)))

# Coverage by prevalence x N (only for original 9 prevalence levels for table comparability)
orig_prevs <- c(0.01, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.99)
coverage <- d2[prevalence %in% orig_prevs, .(
  wald_coverage = round(100 * mean(wald_covers), 1),
  logit_coverage = round(100 * mean(logit_covers), 1),
  n_reps = .N
), by = .(prevalence, grid_N)]
setorder(coverage, prevalence, grid_N)

# Wide format for paper table
wald_wide <- dcast(coverage, prevalence ~ grid_N, value.var = "wald_coverage")
logit_wide <- dcast(coverage, prevalence ~ grid_N, value.var = "logit_coverage")

cat("\nWald coverage (%):\n")
print(wald_wide)
cat("\nLogit coverage (%):\n")
print(logit_wide)

saveRDS(coverage, "output/unified_sim/ci_coverage.rds")

# =====================================================================
# PART 3: Sensitivity analysis — strict and lenient ground-truth defs
# =====================================================================
cat("\n=== SENSITIVITY ANALYSIS ===\n")

# L-K label function
lk_label <- function(x) {
  cut(x, breaks = c(-Inf, 0, 0.20, 0.40, 0.60, 0.80, Inf),
      labels = c("poor", "slight", "fair", "moderate", "substantial", "almost_perfect"),
      right = TRUE)
}
lk_to_band <- function(lk) {
  ifelse(lk %in% c("substantial", "almost_perfect"), "high",
  ifelse(lk %in% c("moderate", "fair"), "medium", "low"))
}

# Three ground-truth definitions
definitions <- list(
  standard = list(high = 0.90, medium = 0.70, label = "Standard (0.90/0.70)"),
  strict   = list(high = 0.95, medium = 0.75, label = "Strict (0.95/0.75)"),
  lenient  = list(high = 0.85, medium = 0.65, label = "Lenient (0.85/0.65)")
)

sensitivity_results <- list()
for (def_name in names(definitions)) {
  def <- definitions[[def_name]]

  # Recompute quality band for each definition
  gt <- grid[, .(scenario_id, min_accuracy)]
  gt[, gt_band := ifelse(min_accuracy >= def$high, "high",
                  ifelse(min_accuracy >= def$medium, "medium", "low"))]

  summ_def <- merge(summaries[, .(scenario_id, kappa_mean)], gt[, .(scenario_id, gt_band)], by = "scenario_id")
  summ_def[, lk_kappa := lk_label(kappa_mean)]
  summ_def[, lk_band := lk_to_band(lk_kappa)]
  summ_def[, lk_correct := lk_band == gt_band]

  acc <- mean(summ_def$lk_correct)
  n_high <- sum(summ_def$gt_band == "high")
  n_med  <- sum(summ_def$gt_band == "medium")
  n_low  <- sum(summ_def$gt_band == "low")

  cat(sprintf("\n%s: L-K kappa accuracy = %.1f%% (H=%d, M=%d, L=%d)\n",
              def$label, 100 * acc, n_high, n_med, n_low))

  sensitivity_results[[def_name]] <- data.table(
    definition = def$label,
    accuracy = round(100 * acc, 1),
    n_high = n_high, n_medium = n_med, n_low = n_low
  )
}

sens_dt <- rbindlist(sensitivity_results)
print(sens_dt)

saveRDS(sens_dt, "output/unified_sim/sensitivity_analysis.rds")

cat("\nAll analyses complete.\n")
