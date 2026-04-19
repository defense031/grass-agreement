# 01_data_generating.R — Binary rating data-generating process
#
# Generates a 2x2 agreement table from a known rater-accuracy model.
# Latent true class Y ~ Bernoulli(prevalence).
# Each rater's observed rating is drawn conditional on Y via Se/Sp.
# Raters are conditionally independent given Y.

generate_binary_ratings <- function(N, prevalence, Se1, Sp1, Se2, Sp2) {
  # Latent true class
  Y <- rbinom(N, 1L, prevalence)
  n_pos <- sum(Y)
  n_neg <- N - n_pos

  # Rater 1
  X1 <- integer(N)
  if (n_pos > 0) X1[Y == 1L] <- rbinom(n_pos, 1L, Se1)
  if (n_neg > 0) X1[Y == 0L] <- rbinom(n_neg, 1L, 1 - Sp1)

  # Rater 2
  X2 <- integer(N)
  if (n_pos > 0) X2[Y == 1L] <- rbinom(n_pos, 1L, Se2)
  if (n_neg > 0) X2[Y == 0L] <- rbinom(n_neg, 1L, 1 - Sp2)

  # Always return a full 2x2 table (factor levels guarantee no missing cells)
  table(factor(X1, levels = 0:1), factor(X2, levels = 0:1))
}
