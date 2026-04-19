source(test_path("fixtures", "published-tables.R"))

# Hand-computed expected values from the cell counts.
# P0 = (a+d)/N;  Pe = (row1*col1 + row2*col2)/N^2
# kappa = (P0 - Pe)/(1 - Pe);  PABAK = 2*P0 - 1
# AC1 uses Pe_AC1 = 2 * pi_hat * (1 - pi_hat) with pi_hat = ((row2 + col2)/2)/N

expected_values <- function(tab) {
  a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
  N <- a + b + c + d
  P0 <- (a + d) / N
  row1 <- a + b; row2 <- c + d
  col1 <- a + c; col2 <- b + d
  Pe <- (row1 * col1 + row2 * col2) / N^2
  kappa <- if (Pe == 1) 0 else (P0 - Pe) / (1 - Pe)
  PABAK <- 2 * P0 - 1
  pi_hat <- (row2 / N + col2 / N) / 2
  Pe_ac1 <- 2 * pi_hat * (1 - pi_hat)
  AC1 <- if (Pe_ac1 == 1) 0 else (P0 - Pe_ac1) / (1 - Pe_ac1)
  list(N = N, P0 = P0, Pe = Pe, kappa = kappa, PABAK = PABAK, AC1 = AC1)
}

test_that("Cohen 1960 Table 1 reproduces kappa ~= 0.76", {
  v <- compute_agreement_metrics(fixture_cohen_1960)
  expect_equal(unname(v["N"]), 200)
  expect_equal(unname(v["P0"]), 0.88, tolerance = tol_metric)
  expect_equal(unname(v["kappa"]), 0.7601, tolerance = 5e-4)
  expect_equal(unname(v["PABAK"]), 0.76, tolerance = tol_metric)
  expect_equal(unname(v["AC1"]), 0.76, tolerance = tol_metric)
})

test_that("Feinstein & Cicchetti 1990 example: high P0 with near-zero/negative kappa", {
  v <- compute_agreement_metrics(fixture_fc_highP)
  expect_equal(unname(v["P0"]),    0.80, tolerance = tol_metric)
  expect_equal(unname(v["kappa"]), -0.1111, tolerance = 1e-3)
  expect_equal(unname(v["PABAK"]), 0.60, tolerance = tol_metric)
  expect_true(v["PABAK"] - v["kappa"] > 0.5)
})

test_that("Byrt-Bishop-Carlin 1993: high PABAK, low kappa at extreme prevalence", {
  v <- compute_agreement_metrics(fixture_bbc_1993)
  expect_true(v["PABAK"] > 0.8)
  expect_true(v["kappa"] < 0.2)
})

test_that("Perfect agreement yields all metrics = 1", {
  v <- compute_agreement_metrics(fixture_perfect)
  expect_equal(unname(v["P0"]),    1)
  expect_equal(unname(v["kappa"]), 1)
  expect_equal(unname(v["PABAK"]), 1)
  expect_equal(unname(v["AC1"]),   1)
})

test_that("Balanced complete disagreement yields P0 = 0, negative kappa/PABAK", {
  v <- compute_agreement_metrics(fixture_disagree)
  expect_equal(unname(v["P0"]),    0)
  expect_equal(unname(v["PABAK"]), -1)
  expect_true(v["kappa"] <= -0.9)
})

test_that("Metric formulas agree with hand computation across all fixtures", {
  for (nm in names(all_fixtures)) {
    tab <- all_fixtures[[nm]]
    v <- compute_agreement_metrics(tab)
    exp <- expected_values(tab)
    expect_equal(unname(v["N"]),     exp$N,     info = nm)
    expect_equal(unname(v["P0"]),    exp$P0,    tolerance = tol_metric, info = nm)
    expect_equal(unname(v["kappa"]), exp$kappa, tolerance = tol_metric, info = nm)
    expect_equal(unname(v["PABAK"]), exp$PABAK, tolerance = tol_metric, info = nm)
    expect_equal(unname(v["AC1"]),   exp$AC1,   tolerance = tol_metric, info = nm)
  }
})

test_that("grass_metrics has consistent N between $n slot and values['N']", {
  # Divergence guard: these two are independently maintained (lifted from the
  # simulation panel vs. ergonomic top-level slot). If they ever drift, tests
  # catch it before a user does.
  for (nm in names(all_fixtures)) {
    m <- grass_compute(all_fixtures[[nm]], format = "matrix")
    expect_equal(m$n, unname(m$values["N"]), info = nm)
  }
})
