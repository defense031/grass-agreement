#!/usr/bin/env Rscript
# run_simulation.R — Main entry point for PABAK Monte Carlo simulation
#
# Usage: Rscript scripts/run_simulation.R
# Or source interactively from the project root.

# Source all function files
source("R/00_packages.R")
source("R/01_data_generating.R")
source("R/02_metrics.R")
source("R/03_ground_truth.R")
source("R/04_grid.R")
source("R/05_runner.R")

# Load configuration
config <- yaml::read_yaml("config.yaml")

# Build parameter grid
message("Building parameter grid...")
grid <- build_parameter_grid(config)
message(sprintf("Grid: %d scenarios (%d standard + %d stress)",
                nrow(grid),
                sum(grid$n_reps == config$n_reps),
                sum(grid$n_reps == config$n_reps_stress)))

# Save grid for reference
dir.create("output", showWarnings = FALSE)
saveRDS(grid, "output/parameter_grid.rds")

# Run simulation
t0 <- Sys.time()
results <- run_all_scenarios(grid, config)
t1 <- Sys.time()

# Save consolidated results
saveRDS(results, "output/sim_results/all_results.rds")

# Summary
message(sprintf("\n=== Simulation Complete ==="))
message(sprintf("Total replications: %s", format(nrow(results), big.mark = ",")))
message(sprintf("Total scenarios:    %d", nrow(grid)))
message(sprintf("Elapsed time:       %s", format(round(difftime(t1, t0, units = "mins"), 1))))
message(sprintf("Results saved to:   output/sim_results/all_results.rds"))
message(sprintf("Grid saved to:      output/parameter_grid.rds"))
message(sprintf("Session: R %s.%s on %s",
                R.version$major, R.version$minor, Sys.info()["nodename"]))
