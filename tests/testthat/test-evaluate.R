source(test_path("fixtures", "published-tables.R"))

test_that("grass_report returns a grass_result with expected slots", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  expect_s3_class(r, "grass_result")
  expect_s3_class(r$metrics, "grass_metrics")
  expect_true(r$regime %in% c("balanced", "prevalence-dominated",
                              "bias-dominated", "mixed"))
  expect_s3_class(r$reference, "grass_reference")
  expect_s3_class(r$distance, "data.frame")
  expect_setequal(r$distance$metric, c("kappa", "PABAK", "AC1"))
})

test_that("grass_report has no verdict column anywhere", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  expect_null(r$verdict)
  expect_false("verdict" %in% names(r$distance))
  expect_false("verdict" %in% names(as.data.frame(r)))
})

test_that("Regime classification matches the PI / BI structure", {
  r_bal  <- grass_report(fixture_cohen_1960, format = "matrix")
  r_prev <- grass_report(fixture_fc_highP,   format = "matrix")
  r_ext  <- grass_report(fixture_bbc_1993,   format = "matrix")

  expect_equal(r_bal$regime,  "balanced")
  expect_equal(r_prev$regime, "prevalence-dominated")
  expect_equal(r_ext$regime,  "prevalence-dominated")
})

test_that("reference = 'none' drops the reference curve and distance", {
  r <- grass_report(fixture_cohen_1960, format = "matrix", reference = "none")
  expect_null(r$reference)
  expect_null(r$distance)
})

test_that("User-supplied prevalence overrides the marginal estimate", {
  r <- grass_report(fixture_cohen_1960, format = "matrix", prevalence = 0.1)
  expect_equal(r$prevalence, 0.1)
  expect_equal(r$prevalence_source, "user_supplied")
})

test_that("grass_prevalence agrees with the internal estimator", {
  p <- grass_prevalence(fixture_cohen_1960, format = "matrix")
  expect_equal(p, 0.5, tolerance = 1e-10)
})

test_that("Distance column is a signed numeric, not a category", {
  r <- grass_report(fixture_fc_highP, format = "matrix")
  expect_true(is.numeric(r$distance$distance))
  expect_equal(length(r$distance$distance), 3)
})
