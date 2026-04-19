source(test_path("fixtures", "published-tables.R"))

test_that("grass_format_report returns a single character string", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  s <- grass_format_report(r)
  expect_type(s, "character")
  expect_length(s, 1)
})

test_that("grass_format_report includes all expected components", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  s <- grass_format_report(r)
  for (tok in c("PABAK", "AC1", "PI", "BI", "N", "prevalence", "regime")) {
    expect_true(grepl(tok, s, fixed = TRUE), info = tok)
  }
})

test_that("ascii = TRUE emits 'kappa' instead of Unicode", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  s <- grass_format_report(r, ascii = TRUE)
  expect_true(grepl("kappa", s, fixed = TRUE))
  expect_false(grepl("\u03ba", s, fixed = TRUE))
})

test_that("ascii = FALSE emits Unicode kappa", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  s <- grass_format_report(r, ascii = FALSE)
  expect_true(grepl("\u03ba", s, fixed = TRUE))
})

test_that("grass_format_report errors on non-grass_result input", {
  expect_error(grass_format_report(list(a = 1)), "grass_result")
})
