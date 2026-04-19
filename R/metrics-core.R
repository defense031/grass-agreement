# Agreement metric computation.
#
# Hand-rolled on 2x2 integer matrices for speed and validated against
# irr::kappa2() during simulation development.
#
# Cell layout for tab (from table(factor(R1, 0:1), factor(R2, 0:1))):
#         R2=0  R2=1
#   R1=0 [  a    b  ]
#   R1=1 [  c    d  ]
#
# a = both say 0, d = both say 1, b = R1=0 / R2=1, c = R1=1 / R2=0.

compute_agreement_metrics <- function(tab) {
  a <- tab[1, 1]
  b <- tab[1, 2]
  c <- tab[2, 1]
  d <- tab[2, 2]
  N <- a + b + c + d

  P0 <- (a + d) / N

  row1 <- a + b
  row2 <- c + d
  col1 <- a + c
  col2 <- b + d
  Pe_cohen <- (row1 * col1 + row2 * col2) / (N * N)

  if (Pe_cohen == 1) {
    kappa <- ifelse(P0 == 1, 1, 0)
  } else {
    kappa <- (P0 - Pe_cohen) / (1 - Pe_cohen)
  }

  PABAK <- 2 * P0 - 1

  pos_denom <- 2 * d + b + c
  pos_agreement <- ifelse(pos_denom == 0, NA_real_, 2 * d / pos_denom)

  neg_denom <- 2 * a + b + c
  neg_agreement <- ifelse(neg_denom == 0, NA_real_, 2 * a / neg_denom)

  prevalence_index <- abs(d - a) / N
  bias_index <- abs(b - c) / N

  pi_hat <- (row2 / N + col2 / N) / 2
  Pe_ac1 <- 2 * pi_hat * (1 - pi_hat)
  if (Pe_ac1 == 1) {
    AC1 <- ifelse(P0 == 1, 1, 0)
  } else {
    AC1 <- (P0 - Pe_ac1) / (1 - Pe_ac1)
  }

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

  if (!is.na(kappa_se_wald) && kappa_se_wald > 0 &&
      kappa > -1 && kappa < 1) {
    kappa_01 <- (kappa + 1) / 2
    logit_k <- log(kappa_01 / (1 - kappa_01))
    se_logit <- kappa_se_wald / (kappa_01 * (1 - kappa_01) * 2)
    logit_lower <- logit_k - 1.96 * se_logit
    logit_upper <- logit_k + 1.96 * se_logit
    kappa_wilson_lower <- 2 * plogis(logit_lower) - 1
    kappa_wilson_upper <- 2 * plogis(logit_upper) - 1
  } else {
    kappa_wilson_lower <- kappa_wald_lower
    kappa_wilson_upper <- kappa_wald_upper
  }

  c(N = N,
    P0 = P0,
    Pe = Pe_cohen,
    kappa = kappa,
    PABAK = PABAK,
    AC1 = AC1,
    pos_agreement = pos_agreement,
    neg_agreement = neg_agreement,
    prevalence_index = prevalence_index,
    bias_index = bias_index,
    kappa_se_wald = kappa_se_wald,
    kappa_wald_lower = kappa_wald_lower,
    kappa_wald_upper = kappa_wald_upper,
    kappa_wilson_lower = kappa_wilson_lower,
    kappa_wilson_upper = kappa_wilson_upper)
}
