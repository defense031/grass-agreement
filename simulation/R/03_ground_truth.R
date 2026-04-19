# 03_ground_truth.R — Ground-truth labels and theoretical metrics
#
# Ground truth is defined at the SCENARIO level (Se/Sp/prevalence),
# not at the replicate level. A high-accuracy profile is "high quality"
# regardless of any single noisy replicate.

compute_ground_truth <- function(Se1, Sp1, Se2, Sp2, prevalence, cost_ratios) {

  # --- Quality band from minimum operating characteristic ---
  min_acc <- min(Se1, Sp1, Se2, Sp2)
  quality_band <- if (min_acc >= 0.90) {
    "high"
  } else if (min_acc >= 0.70) {
    "medium"
  } else {
    "low"
  }

  # --- Theoretical population P0 (conditional independence model) ---
  # P(agree) = P(both correct) + P(both wrong)
  # P(agree | Y=1) = Se1*Se2 + (1-Se1)*(1-Se2)
  # P(agree | Y=0) = Sp1*Sp2 + (1-Sp1)*(1-Sp2)
  P0_agree_pos <- Se1 * Se2 + (1 - Se1) * (1 - Se2)
  P0_agree_neg <- Sp1 * Sp2 + (1 - Sp1) * (1 - Sp2)
  P0_theoretical <- prevalence * P0_agree_pos + (1 - prevalence) * P0_agree_neg

  # --- Theoretical marginals for Pe ---
  # P(R1=1) = prevalence*Se1 + (1-prevalence)*(1-Sp1)
  # P(R2=1) = prevalence*Se2 + (1-prevalence)*(1-Sp2)
  p_R1_pos <- prevalence * Se1 + (1 - prevalence) * (1 - Sp1)
  p_R2_pos <- prevalence * Se2 + (1 - prevalence) * (1 - Sp2)
  p_R1_neg <- 1 - p_R1_pos
  p_R2_neg <- 1 - p_R2_pos

  Pe_theoretical <- p_R1_pos * p_R2_pos + p_R1_neg * p_R2_neg

  # --- Theoretical kappa ---
  if (Pe_theoretical == 1) {
    kappa_theoretical <- ifelse(P0_theoretical == 1, 1, 0)
  } else {
    kappa_theoretical <- (P0_theoretical - Pe_theoretical) / (1 - Pe_theoretical)
  }

  PABAK_theoretical <- 2 * P0_theoretical - 1

  # --- Expected decision loss ---
  # Two decision rules:
  #   "require_agreement": classify positive only if BOTH raters say 1
  #     Se_eff = Se1 * Se2, Sp_eff = 1 - (1-Sp1)*(1-Sp2)
  #   "any_positive": classify positive if EITHER rater says 1
  #     Se_eff = 1 - (1-Se1)*(1-Se2), Sp_eff = Sp1 * Sp2

  Se_require <- Se1 * Se2
  Sp_require <- 1 - (1 - Sp1) * (1 - Sp2)
  Se_any     <- 1 - (1 - Se1) * (1 - Se2)
  Sp_any     <- Sp1 * Sp2

  # --- Theoretical PI and BI ---
  # Cell probabilities under conditional independence
  prob_a <- prevalence * (1 - Se1) * (1 - Se2) + (1 - prevalence) * Sp1 * Sp2
  prob_b <- prevalence * (1 - Se1) * Se2 + (1 - prevalence) * Sp1 * (1 - Sp2)
  prob_c <- prevalence * Se1 * (1 - Se2) + (1 - prevalence) * (1 - Sp1) * Sp2
  prob_d <- prevalence * Se1 * Se2 + (1 - prevalence) * (1 - Sp1) * (1 - Sp2)
  PI_theoretical <- abs(prob_d - prob_a)
  BI_theoretical <- abs(prob_b - prob_c)
  delta_theoretical <- PI_theoretical^2 - BI_theoretical^2

  # Build result row
  result <- data.frame(
    quality_band       = quality_band,
    min_accuracy       = min_acc,
    P0_theoretical     = P0_theoretical,
    Pe_theoretical     = Pe_theoretical,
    kappa_theoretical  = kappa_theoretical,
    PABAK_theoretical  = PABAK_theoretical,
    PI_theoretical     = PI_theoretical,
    BI_theoretical     = BI_theoretical,
    delta_theoretical  = delta_theoretical,
    Se_require         = Se_require,
    Sp_require         = Sp_require,
    Se_any             = Se_any,
    Sp_any             = Sp_any,
    stringsAsFactors   = FALSE
  )

  # Expected loss for each cost ratio:
  # loss = prevalence * (1 - Se_eff) * cost_ratio + (1 - prevalence) * (1 - Sp_eff)
  for (cr in cost_ratios) {
    cr_label <- gsub("\\.", "p", as.character(cr))
    result[[paste0("loss_require_cr", cr_label)]] <-
      prevalence * (1 - Se_require) * cr + (1 - prevalence) * (1 - Sp_require)
    result[[paste0("loss_any_cr", cr_label)]] <-
      prevalence * (1 - Se_any) * cr + (1 - prevalence) * (1 - Sp_any)
  }

  result
}
