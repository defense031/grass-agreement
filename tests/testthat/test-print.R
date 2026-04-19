source(test_path("fixtures", "published-tables.R"))

test_that("print.grass_result shows metrics, skew diagnostics, and regime", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  out <- capture.output(print(r))
  full <- paste(out, collapse = "\n")
  expect_match(full, "grass report")
  expect_match(full, "Skew diagnostics")
  expect_match(full, "Regime")
  expect_false(grepl("\\[ABOVE\\]|\\[BELOW\\]", full))
})

test_that("print.grass_result runs on all fixtures", {
  for (nm in names(all_fixtures)) {
    tab <- all_fixtures[[nm]]
    r <- suppressWarnings(grass_report(tab, format = "matrix"))
    expect_output(print(r), "grass report", info = nm)
  }
})

test_that("print.grass_metrics shows the 2x2 table", {
  m <- grass_compute(fixture_cohen_1960, format = "matrix")
  expect_output(print(m), "2x2 table")
})

test_that("print.grass_reference labels curves as references", {
  t <- grass_reference(0.3)
  expect_output(print(t), "reference")
})

test_that("as.data.frame.grass_result returns context columns, no verdict", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  df <- as.data.frame(r)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("N", "prevalence", "regime",
                    "kappa", "PABAK", "AC1",
                    "prevalence_index", "bias_index") %in% names(df)))
  expect_false("verdict" %in% names(df))
})
