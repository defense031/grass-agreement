# 05_runner.R — Monte Carlo simulation runner
#
# run_scenario(): inner loop for a single scenario (n_reps replications)
# run_all_scenarios(): outer orchestration with future/furrr parallelism

run_scenario <- function(scenario_row, base_seed) {
  # Reproducible seeding per scenario
  set.seed(base_seed + scenario_row$scenario_id * 100000L)

  n_reps <- scenario_row$n_reps
  N      <- scenario_row$N
  prev   <- scenario_row$prevalence
  Se1    <- scenario_row$Se1
  Sp1    <- scenario_row$Sp1
  Se2    <- scenario_row$Se2
  Sp2    <- scenario_row$Sp2

  # Pre-allocate list for speed
  results <- vector("list", n_reps)

  for (i in seq_len(n_reps)) {
    tab <- generate_binary_ratings(N, prev, Se1, Sp1, Se2, Sp2)
    metrics <- compute_agreement_metrics(tab)
    results[[i]] <- as.list(metrics)
  }

  # Fast row-binding via data.table
  dt <- data.table::rbindlist(results)
  dt[, rep_id := seq_len(.N)]
  dt[, scenario_id := scenario_row$scenario_id]

  dt
}

run_all_scenarios <- function(grid, config) {
  base_seed <- config$seed

  # Set up parallel backend
  n_cores <- parallel::detectCores() - 1L
  n_cores <- max(n_cores, 1L)
  future::plan(future::multisession, workers = n_cores)
  on.exit(future::plan(future::sequential), add = TRUE)

  message(sprintf("Running %d scenarios across %d cores...", nrow(grid), n_cores))

  # Split grid into batches for fault tolerance
  batch_size <- 50L
  n_batches  <- ceiling(nrow(grid) / batch_size)
  batch_ids  <- rep(seq_len(n_batches), each = batch_size, length.out = nrow(grid))

  # Ensure output directory exists
  out_dir <- file.path("output", "sim_results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  all_results <- vector("list", n_batches)

  for (b in seq_len(n_batches)) {
    batch_rows <- which(batch_ids == b)
    batch_grid <- grid[batch_rows, ]

    # Convert each row to a list for furrr
    row_list <- split(batch_grid, seq_len(nrow(batch_grid)))

    batch_results <- furrr::future_map(
      row_list,
      function(row) run_scenario(row, base_seed),
      .options = furrr::furrr_options(seed = NULL)
    )

    batch_dt <- data.table::rbindlist(batch_results)

    # Save per-batch for fault tolerance
    batch_file <- file.path(out_dir, sprintf("batch_%03d.rds", b))
    saveRDS(batch_dt, batch_file)

    all_results[[b]] <- batch_dt

    message(sprintf("Batch %d/%d complete (%d scenarios)",
                    b, n_batches, nrow(batch_grid)))
  }

  # Consolidate
  final <- data.table::rbindlist(all_results)
  message(sprintf("Simulation complete: %d total replications across %d scenarios.",
                  nrow(final), nrow(grid)))

  final
}
