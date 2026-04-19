# Published 2x2 agreement tables used as test fixtures.
#
# Cell layout for every matrix below:
#         R2=0  R2=1
#   R1=0 [ a    b  ]
#   R1=1 [ c    d  ]
#
# Constructed via matrix(c(a, c, b, d), nrow = 2) so column-major fill matches
# the same (R1, R2) cell labeling used in grass internals.

make_tab <- function(a, b, c, d) {
  matrix(as.integer(c(a, c, b, d)), nrow = 2,
         dimnames = list(R1 = c("0", "1"), R2 = c("0", "1")))
}

# Cohen (1960), Table 1 worked example. N = 200.
# Expected: kappa ~= 0.76, PABAK = 0.76, AC1 = 0.76 (all agree at balanced P).
fixture_cohen_1960 <- make_tab(a = 88, b = 14, c = 10, d = 88)

# Feinstein & Cicchetti (1990) high-prevalence example. High P0, negative
# kappa — a case where high raw agreement coexists with kappa near zero.
# a=40, b=5, c=5, d=0. N=50. P0=0.80, kappa ~= -0.111, PABAK = 0.60.
fixture_fc_highP <- make_tab(a = 40, b = 5, c = 5, d = 0)

# Byrt, Bishop & Carlin (1993) example. PABAK markedly higher than kappa
# when prevalence is extreme.  N = 125.  a=118, b=5, c=2, d=0 gives P0 = 0.944
# but kappa near 0 because one cell is empty.
fixture_bbc_1993 <- make_tab(a = 118, b = 5, c = 2, d = 0)

# Landis & Koch (1977) Table 1 example. N = 100. a=45, b=10, c=5, d=40.
fixture_landis_koch <- make_tab(a = 45, b = 10, c = 5, d = 40)

# Perfect agreement.
fixture_perfect <- make_tab(a = 50, b = 0, c = 0, d = 50)

# Complete disagreement (balanced).
fixture_disagree <- make_tab(a = 0, b = 50, c = 50, d = 0)

# One-cell dominant (extreme prevalence, tiny minority class).
fixture_extreme <- make_tab(a = 180, b = 2, c = 3, d = 15)

all_fixtures <- list(
  cohen_1960   = fixture_cohen_1960,
  fc_highP     = fixture_fc_highP,
  bbc_1993     = fixture_bbc_1993,
  landis_koch  = fixture_landis_koch,
  perfect      = fixture_perfect,
  disagree     = fixture_disagree,
  extreme      = fixture_extreme
)

# Tolerance used in expect_equal for metric comparisons.
tol_metric <- 1e-4
