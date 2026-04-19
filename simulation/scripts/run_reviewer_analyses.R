# run_reviewer_analyses.R — Jackknife threshold stability + CI coverage
# Addresses reviewer concerns about threshold generalizability and CI quality

library(data.table)
library(yaml)

# Source project files
source("R/00_packages.R")
source("R/03_ground_truth.R")
source("R/04_grid.R")

cfg <- read_yaml("config.yaml")
grid <- build_parameter_grid(cfg)
setDT(grid)

# Load replicate-level data
d <- readRDS("output/sim_results/all_results.rds")
setDT(d)

# Merge with grid for profile_name, prevalence, quality_band, kappa_theoretical
grid_cols <- grid[, .(scenario_id, prevalence, profile_name, quality_band,
                       kappa_theoretical, grid_N = N)]
d2 <- merge(d, grid_cols, by = "scenario_id")

# =====================================================================
# PART 1: Jackknife cross-validation of ROC-optimal kappa thresholds
# =====================================================================

cat("=== JACKKNIFE THRESHOLD STABILITY ===\n")

# Helper: derive kappa High-vs-not-High threshold via Youden's J
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
    results[[length(results) + 1]] <- data.table(prevalence = p, threshold = best_t, J = round(best_j, 2))
  }
  rbindlist(results)
}

# Compute scenario-level summaries
summaries_full <- d2[, .(kappa_mean = mean(kappa), PABAK_mean = mean(PABAK),
                          AC1_mean = mean(AC1)),
                      by = .(scenario_id, prevalence, profile_name, quality_band, grid_N)]

# Full-grid thresholds (baseline)
full_thresholds <- derive_kappa_threshold(summaries_full)
setnames(full_thresholds, c("threshold", "J"), c("full_threshold", "full_J"))

# Leave-one-profile-out
profiles <- unique(summaries_full$profile_name)
jackknife_results <- list()

for (prof in profiles) {
  sub <- summaries_full[profile_name != prof]
  thresh <- derive_kappa_threshold(sub)
  thresh[, dropped_profile := prof]
  jackknife_results[[length(jackknife_results) + 1]] <- thresh
}

jk <- rbindlist(jackknife_results)

# Compute deviations
jk_summary <- jk[, .(min_thresh = min(threshold, na.rm = TRUE),
                       max_thresh = max(threshold, na.rm = TRUE),
                       mean_thresh = mean(threshold, na.rm = TRUE),
                       sd_thresh = sd(threshold, na.rm = TRUE)),
                   by = prevalence]
jk_summary <- merge(jk_summary, full_thresholds, by = "prevalence")
jk_summary[, max_deviation := pmax(abs(full_threshold - min_thresh), abs(full_threshold - max_thresh))]

cat("\nJackknife results (kappa threshold for High quality):\n")
print(jk_summary[, .(prevalence, full_threshold = round(full_threshold, 3),
                       jk_min = round(min_thresh, 3), jk_max = round(max_thresh, 3),
                       max_dev = round(max_deviation, 3))])

saveRDS(jk_summary, "output/thresholds/jackknife_stability.rds")
saveRDS(jk, "output/thresholds/jackknife_detail.rds")

# =====================================================================
# PART 2: CI coverage proportions
# =====================================================================

cat("\n=== CI COVERAGE PROBABILITIES ===\n")

# Wald coverage
d2[, wald_covers := (!is.na(kappa_wald_lower)) &
     (kappa_wald_lower <= kappa_theoretical) &
     (kappa_wald_upper >= kappa_theoretical)]

# Logit (stored as "wilson" in the data) coverage
d2[, logit_covers := (!is.na(kappa_wilson_lower)) &
     (kappa_wilson_lower <= kappa_theoretical) &
     (kappa_wilson_upper >= kappa_theoretical)]

# Coverage by prevalence x N
coverage <- d2[, .(wald_coverage = round(mean(wald_covers), 3),
                    logit_coverage = round(mean(logit_covers), 3),
                    n_reps = .N),
                by = .(prevalence, grid_N)]
setorder(coverage, prevalence, grid_N)

cat("\nCoverage by prevalence and N:\n")
print(coverage)

# Wide format for the paper table
wald_wide <- dcast(coverage, prevalence ~ grid_N, value.var = "wald_coverage")
logit_wide <- dcast(coverage, prevalence ~ grid_N, value.var = "logit_coverage")

cat("\nWald coverage (wide):\n")
print(wald_wide)
cat("\nLogit coverage (wide):\n")
print(logit_wide)

# Overall
cat(sprintf("\nOverall Wald coverage: %.3f\n", d2[, mean(wald_covers)]))
cat(sprintf("Overall Logit coverage: %.3f\n", d2[, mean(logit_covers)]))

saveRDS(coverage, "output/thresholds/ci_coverage.rds")

cat("\nDone.\n")
