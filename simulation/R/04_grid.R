# 04_grid.R — Parameter grid construction
#
# Builds the full factorial grid (495 scenarios) from config.yaml
# and pre-computes ground-truth columns for each scenario.

build_parameter_grid <- function(config) {
  # Extract profile Se/Sp values into a data.frame
  profiles <- do.call(rbind, lapply(names(config$rater_profiles), function(pname) {
    p <- config$rater_profiles[[pname]]
    data.frame(
      profile_name = pname,
      Se1 = p$Se1, Sp1 = p$Sp1,
      Se2 = p$Se2, Sp2 = p$Sp2,
      stringsAsFactors = FALSE
    )
  }))

  # Full factorial: prevalence x sample_size x profile
  grid <- expand.grid(
    prevalence   = config$prevalences,
    N            = config$sample_sizes,
    profile_name = profiles$profile_name,
    stringsAsFactors = FALSE
  )

  # Join Se/Sp values
  grid <- merge(grid, profiles, by = "profile_name", sort = FALSE)

  # Add scenario ID
  grid$scenario_id <- seq_len(nrow(grid))

  # Determine n_reps: stress conditions get more reps
  stress_prev <- config$stress_conditions$prevalences
  stress_n    <- config$stress_conditions$sample_sizes
  grid$n_reps <- ifelse(
    grid$prevalence %in% stress_prev & grid$N %in% stress_n,
    config$n_reps_stress,
    config$n_reps
  )

  # Pre-compute ground truth for each scenario
  gt_list <- mapply(
    compute_ground_truth,
    Se1 = grid$Se1, Sp1 = grid$Sp1,
    Se2 = grid$Se2, Sp2 = grid$Sp2,
    prevalence = grid$prevalence,
    MoreArgs = list(cost_ratios = config$cost_ratios),
    SIMPLIFY = FALSE
  )
  gt_df <- do.call(rbind, gt_list)

  # Bind ground truth columns to grid
  grid <- cbind(grid, gt_df)

  # Reorder columns for readability
  front_cols <- c("scenario_id", "profile_name", "prevalence", "N",
                  "Se1", "Sp1", "Se2", "Sp2", "n_reps")
  other_cols <- setdiff(names(grid), front_cols)
  grid <- grid[, c(front_cols, other_cols)]

  grid
}

# Build balanced grid for the bias stress test (supplementary simulation).
# Non-factorial: tiered prevalence sets by regime and quality band.
build_balanced_grid <- function(config) {

  # Helper: extract profiles from a named list into a data.frame
  extract_profiles <- function(profile_list) {
    do.call(rbind, lapply(names(profile_list), function(pname) {
      p <- profile_list[[pname]]
      data.frame(profile_name = pname,
                 Se1 = p$Se1, Sp1 = p$Sp1,
                 Se2 = p$Se2, Sp2 = p$Sp2,
                 stringsAsFactors = FALSE)
    }))
  }

  # --- Prevalence-dominant sub-grid ---
  pd_profiles <- extract_profiles(config$prevalence_dominant$profiles)
  pd_grid <- expand.grid(
    prevalence   = config$prevalence_dominant$prevalences,
    N            = config$sample_sizes,
    profile_name = pd_profiles$profile_name,
    stringsAsFactors = FALSE
  )
  pd_grid <- merge(pd_grid, pd_profiles, by = "profile_name", sort = FALSE)
  pd_grid$regime <- "prevalence_dominant"
  pd_grid$tier   <- "prevalence_dominant"


  # --- Bias-dominant sub-grids (one per quality tier) ---
  bd_grids <- list()
  for (tier_name in c("high", "medium", "low")) {
    tier_cfg  <- config$bias_dominant[[tier_name]]
    tier_prof <- extract_profiles(tier_cfg$profiles)
    tier_grid <- expand.grid(
      prevalence   = tier_cfg$prevalences,
      N            = config$sample_sizes,
      profile_name = tier_prof$profile_name,
      stringsAsFactors = FALSE
    )
    tier_grid <- merge(tier_grid, tier_prof, by = "profile_name", sort = FALSE)
    tier_grid$regime <- "bias_dominant"
    tier_grid$tier   <- paste0("bias_dominant_", tier_name)
    bd_grids[[tier_name]] <- tier_grid
  }

  # --- Combine ---
  grid <- rbind(pd_grid, do.call(rbind, bd_grids))
  grid$scenario_id <- seq_len(nrow(grid))

  # Determine n_reps
  stress_prev <- config$stress_conditions$prevalences
  stress_n    <- config$stress_conditions$sample_sizes
  grid$n_reps <- ifelse(
    grid$prevalence %in% stress_prev & grid$N %in% stress_n,
    config$n_reps_stress,
    config$n_reps
  )

  # Pre-compute ground truth
  gt_list <- mapply(
    compute_ground_truth,
    Se1 = grid$Se1, Sp1 = grid$Sp1,
    Se2 = grid$Se2, Sp2 = grid$Sp2,
    prevalence = grid$prevalence,
    MoreArgs = list(cost_ratios = config$cost_ratios),
    SIMPLIFY = FALSE
  )
  gt_df <- do.call(rbind, gt_list)
  grid <- cbind(grid, gt_df)

  # --- Verify regime assignments ---
  pd_check <- grid$regime == "prevalence_dominant" & grid$PI_theoretical <= grid$BI_theoretical
  bd_check <- grid$regime == "bias_dominant" & grid$BI_theoretical <= grid$PI_theoretical
  if (any(pd_check)) {
    warning(sprintf("%d prevalence-dominant scenarios have PI <= BI!", sum(pd_check)))
  }
  if (any(bd_check)) {
    warning(sprintf("%d bias-dominant scenarios have BI <= PI!", sum(bd_check)))
  }

  # Reorder columns
  front_cols <- c("scenario_id", "profile_name", "regime", "tier",
                  "prevalence", "N", "Se1", "Sp1", "Se2", "Sp2", "n_reps")
  other_cols <- setdiff(names(grid), front_cols)
  grid <- grid[, c(front_cols, other_cols)]

  grid
}
