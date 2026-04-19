# 02_metrics.R — Agreement metric computation
#
# Hand-rolled on 2x2 integer matrices for speed (~2.5M calls).
# Validated against irr::kappa2() before full simulation run.
#
# Cell layout for tab (from table(factor(X1,0:1), factor(X2,0:1))):
#         X2=0  X2=1
#   X1=0 [ a     b  ]
#   X1=1 [ c     d  ]
#
# a = both say 0, d = both say 1, b = R1=0/R2=1, c = R1=1/R2=0

compute_agreement_metrics <- function(tab) {
  a <- tab[1, 1]
  b <- tab[1, 2]
  c <- tab[2, 1]
  d <- tab[2, 2]
  N <- a + b + c + d

  # Observed agreement
  P0 <- (a + d) / N

  # Expected agreement (Cohen's model, marginal-based)
  row1 <- a + b  # R1 says 0
  row2 <- c + d  # R1 says 1
  col1 <- a + c  # R2 says 0
  col2 <- b + d  # R2 says 1
  Pe_cohen <- (row1 * col1 + row2 * col2) / (N * N)

  # Cohen's kappa
  if (Pe_cohen == 1) {
    kappa <- ifelse(P0 == 1, 1, 0)
  } else {
    kappa <- (P0 - Pe_cohen) / (1 - Pe_cohen)
  }

  # PABAK (= Brennan-Prediger for binary)
  PABAK <- 2 * P0 - 1

  # Positive agreement: 2d / (2d + b + c)
  pos_denom <- 2 * d + b + c
  pos_agreement <- ifelse(pos_denom == 0, NA_real_, 2 * d / pos_denom)

  # Negative agreement: 2a / (2a + b + c)
  neg_denom <- 2 * a + b + c
  neg_agreement <- ifelse(neg_denom == 0, NA_real_, 2 * a / neg_denom)

  # Byrt et al. prevalence and bias indices
  prevalence_index <- abs(d - a) / N
  bias_index <- abs(b - c) / N

  # Gwet's AC1 for binary ratings
  # Pe_AC1 = 2 * pi_hat * (1 - pi_hat), where pi_hat = average positive rate
  pi_hat <- (row2 / N + col2 / N) / 2  # average P(rater says 1)
  Pe_ac1 <- 2 * pi_hat * (1 - pi_hat)
  if (Pe_ac1 == 1) {
    AC1 <- ifelse(P0 == 1, 1, 0)
  } else {
    AC1 <- (P0 - Pe_ac1) / (1 - Pe_ac1)
  }

  # Kappa SE — Wald (Fleiss large-sample formula)
  p1p <- row1 / N
  p2p <- row2 / N
  pp1 <- col1 / N
  pp2 <- col2 / N

  sum_term <- p1p * pp1 * (p1p + pp1) + p2p * pp2 * (p2p + pp2)
  se_num <- Pe_cohen + Pe_cohen^2 - sum_term
  se_denom <- N * (1 - Pe_cohen)^2

  if (se_denom <= 0 || se_num < 0) {
    kappa_se_wald <- NA_real_
  } else {
    kappa_se_wald <- sqrt(se_num / se_denom)
  }

  kappa_wald_lower <- kappa - 1.96 * ifelse(is.na(kappa_se_wald), 0, kappa_se_wald)
  kappa_wald_upper <- kappa + 1.96 * ifelse(is.na(kappa_se_wald), 0, kappa_se_wald)

  # Kappa CI — Wilson-type (variance-stabilized via log transform)
  # Uses the delta-method on the logit of (kappa+1)/2 to respect [-1,1] bounds
  if (!is.na(kappa_se_wald) && kappa_se_wald > 0 &&
      kappa > -1 && kappa < 1) {
    kappa_01 <- (kappa + 1) / 2  # map to (0,1)
    logit_k <- log(kappa_01 / (1 - kappa_01))
    # Delta method: SE(logit) = SE(kappa) / (kappa_01 * (1-kappa_01) * 2)
    se_logit <- kappa_se_wald / (kappa_01 * (1 - kappa_01) * 2)
    logit_lower <- logit_k - 1.96 * se_logit
    logit_upper <- logit_k + 1.96 * se_logit
    # Back-transform
    kappa_wilson_lower <- 2 * plogis(logit_lower) - 1
    kappa_wilson_upper <- 2 * plogis(logit_upper) - 1
  } else {
    kappa_wilson_lower <- kappa_wald_lower
    kappa_wilson_upper <- kappa_wald_upper
  }

  c(N       = N,
    P0      = P0,
    Pe      = Pe_cohen,
    kappa   = kappa,
    PABAK   = PABAK,
    AC1     = AC1,
    pos_agreement    = pos_agreement,
    neg_agreement    = neg_agreement,
    prevalence_index = prevalence_index,
    bias_index       = bias_index,
    kappa_se_wald    = kappa_se_wald,
    kappa_wald_lower = kappa_wald_lower,
    kappa_wald_upper = kappa_wald_upper,
    kappa_wilson_lower = kappa_wilson_lower,
    kappa_wilson_upper = kappa_wilson_upper)
}
