# Generate the internal `reference_binary` sysdata object.
#
# For each rater-quality band q in {0.70, 0.80, 0.85, 0.90} and each
# prevalence in the standard 23-point grid, compute the expected value
# of Cohen's kappa, PABAK, and Gwet's AC1 on the Se = Sp = q diagonal
# under conditional independence of the two raters given the true class.
#
# This is the analytical calibration reference: the metric value a pair
# of raters at exactly Se = Sp = q should achieve in expectation at the
# given prevalence. Youden's J at the diagonal is constant: J = 2q - 1.
#
# NOTE on estimand change from pre-0.1.2: prior releases shipped a
# simulation-derived ROC threshold (a classifier cutpoint separating
# scenarios of different ground-truth quality bands). The 0.1.2
# reference is a different estimand — the expected metric value at an
# exact calibration point — and is numerically different. The new
# estimand has a closed-form interpretation, extends uniformly to
# additional bands (0.80, 0.90), and is what most users expect when
# they see "reference at Se = Sp = 0.85".

library(usethis)

prevalence_grid <- c(0.01, 0.05, seq(0.10, 0.90, by = 0.05), 0.95, 0.99)
quality_levels  <- c(0.70, 0.80, 0.85, 0.90)

# Closed-form agreement metrics under conditional independence of R1, R2
# given the true class D, with Se = Sp = q. Derived once, evaluated at
# each (p, q) grid point.
metric_at_diagonal <- function(p, q) {
  m <- q * p + (1 - q) * (1 - p)
  P_agree_11 <- p * q^2       + (1 - p) * (1 - q)^2
  P_agree_00 <- p * (1 - q)^2 + (1 - p) * q^2
  P_agree <- P_agree_11 + P_agree_00
  Pe_cohen <- m^2 + (1 - m)^2
  Pe_ac1   <- 2 * m * (1 - m)
  kappa <- if (abs(1 - Pe_cohen) < .Machine$double.eps)
    if (abs(P_agree - 1) < .Machine$double.eps) 1 else 0
  else (P_agree - Pe_cohen) / (1 - Pe_cohen)
  PABAK <- 2 * P_agree - 1
  AC1 <- if (abs(1 - Pe_ac1) < .Machine$double.eps)
    if (abs(P_agree - 1) < .Machine$double.eps) 1 else 0
  else (P_agree - Pe_ac1) / (1 - Pe_ac1)
  list(kappa = kappa, PABAK = PABAK, AC1 = AC1)
}

rows <- list()
for (q in quality_levels) {
  J <- 2 * q - 1
  for (p in prevalence_grid) {
    vals <- metric_at_diagonal(p, q)
    for (met in c("kappa", "PABAK", "AC1")) {
      rows[[length(rows) + 1]] <- data.frame(
        prevalence      = p,
        reference_level = q,
        metric          = met,
        reference       = vals[[met]],
        J               = J,
        stringsAsFactors = FALSE
      )
    }
  }
}
reference_binary <- do.call(rbind, rows)
reference_binary <- reference_binary[order(reference_binary$reference_level,
                                            reference_binary$metric,
                                            reference_binary$prevalence), ]
row.names(reference_binary) <- NULL

cat("reference_binary:\n")
cat("  rows =", nrow(reference_binary), "\n")
cat("  bands =", paste(unique(reference_binary$reference_level), collapse = ", "), "\n")
cat("  metrics =", paste(unique(reference_binary$metric), collapse = ", "), "\n")
cat("\nhead(at q = 0.85):\n")
print(head(reference_binary[reference_binary$reference_level == 0.85, ], 10))

usethis::use_data(reference_binary, internal = TRUE, overwrite = TRUE)
