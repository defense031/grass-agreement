#!/usr/bin/env Rscript
# run_balanced_sim.R — Balanced bias stress test (supplementary simulation)
#
# Does NOT modify or replace the original simulation.
# Produces: output/balanced_sim/parameter_grid.rds
#           output/balanced_sim/sim_results/all_results.rds

source("R/00_packages.R")
source("R/01_data_generating.R")
source("R/02_metrics.R")
source("R/03_ground_truth.R")
source("R/04_grid.R")
source("R/05_runner.R")

# Load balanced config
config <- yaml::read_yaml("config_balanced.yaml")

# Build grid
cat("Building balanced grid...\n")
grid <- build_balanced_grid(config)

cat(sprintf("Total scenarios: %d\n", nrow(grid)))
cat(sprintf("  Prevalence-dominant: %d\n", sum(grid$regime == "prevalence_dominant")))
cat(sprintf("  Bias-dominant:       %d\n", sum(grid$regime == "bias_dominant")))
cat(sprintf("  Unique profiles:     %d\n", length(unique(grid$profile_name))))

# Verify regime assignments
cat("\nRegime verification:\n")
cat(sprintf("  PD scenarios with PI > BI: %d / %d\n",
    sum(grid$regime == "prevalence_dominant" & grid$PI_theoretical > grid$BI_theoretical),
    sum(grid$regime == "prevalence_dominant")))
cat(sprintf("  BD scenarios with BI > PI: %d / %d\n",
    sum(grid$regime == "bias_dominant" & grid$BI_theoretical > grid$PI_theoretical),
    sum(grid$regime == "bias_dominant")))

# Save grid
dir.create("output/balanced_sim/sim_results", recursive = TRUE, showWarnings = FALSE)
saveRDS(grid, "output/balanced_sim/parameter_grid.rds")
cat("\nGrid saved to output/balanced_sim/parameter_grid.rds\n")

# Run simulation
cat("\nStarting simulation...\n")
t0 <- Sys.time()
results <- run_all_scenarios(grid, config)
t1 <- Sys.time()

cat(sprintf("\nSimulation complete: %s replications in %.1f minutes\n",
    format(nrow(results), big.mark = ","), as.numeric(t1 - t0, units = "mins")))

# Save results
saveRDS(results, "output/balanced_sim/sim_results/all_results.rds")
cat("Results saved to output/balanced_sim/sim_results/all_results.rds\n")
