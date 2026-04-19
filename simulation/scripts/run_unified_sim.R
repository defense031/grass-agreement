#!/usr/bin/env Rscript
# run_unified_sim.R — Unified full-factorial simulation
# 39 profiles x 15 prevalences x 5 sample sizes = 2,925 scenarios
#
# Produces: output/unified_sim/parameter_grid.rds
#           output/unified_sim/sim_results/all_results.rds

source("R/00_packages.R")
source("R/01_data_generating.R")
source("R/02_metrics.R")
source("R/03_ground_truth.R")
source("R/04_grid.R")
source("R/05_runner.R")

config <- yaml::read_yaml("config_unified.yaml")
cat("Building unified grid...\n")
grid <- build_parameter_grid(config)

# Add regime classification based on theoretical PI/BI
grid$regime <- ifelse(grid$PI_theoretical > grid$BI_theoretical, "prevalence_dominant",
               ifelse(grid$BI_theoretical > grid$PI_theoretical, "bias_dominant", "balanced"))

cat(sprintf("Total scenarios: %d\n", nrow(grid)))
cat(sprintf("  Unique profiles: %d\n", length(unique(grid$profile_name))))
cat(sprintf("  Prevalence levels: %d\n", length(unique(grid$prevalence))))
cat(sprintf("  Sample sizes: %d\n", length(unique(grid$N))))
cat(sprintf("\nRegime distribution:\n"))
cat(sprintf("  Prevalence-dominant: %d (%.1f%%)\n",
    sum(grid$regime == "prevalence_dominant"),
    100 * mean(grid$regime == "prevalence_dominant")))
cat(sprintf("  Bias-dominant:       %d (%.1f%%)\n",
    sum(grid$regime == "bias_dominant"),
    100 * mean(grid$regime == "bias_dominant")))
cat(sprintf("  Balanced (PI = BI):  %d (%.1f%%)\n",
    sum(grid$regime == "balanced"),
    100 * mean(grid$regime == "balanced")))
cat(sprintf("\nQuality band distribution:\n"))
print(table(grid$quality_band))
cat(sprintf("\nTotal replications: %s\n", format(sum(grid$n_reps), big.mark = ",")))

# Save grid
dir.create("output/unified_sim/sim_results", recursive = TRUE, showWarnings = FALSE)
saveRDS(grid, "output/unified_sim/parameter_grid.rds")

# Run simulation
cat("\nStarting simulation...\n")
t0 <- Sys.time()
results <- run_all_scenarios(grid, config)
t1 <- Sys.time()

cat(sprintf("\nSimulation complete: %s replications in %.1f minutes\n",
    format(nrow(results), big.mark = ","), as.numeric(t1 - t0, units = "mins")))

saveRDS(results, "output/unified_sim/sim_results/all_results.rds")
cat("Results saved to output/unified_sim/\n")
