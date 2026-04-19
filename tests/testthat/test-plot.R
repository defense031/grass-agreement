source(test_path("fixtures", "published-tables.R"))

test_that("plot methods return ggplot objects when ggplot2 is available", {
  skip_if_not_installed("ggplot2")
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  expect_s3_class(plot(r),                              "ggplot")
  expect_s3_class(plot(r, type = "regime"),             "ggplot")
  expect_s3_class(plot(r, labels = "inline"),           "ggplot")
  expect_s3_class(plot(r, labels = "legend"),           "ggplot")
  expect_s3_class(plot(r$reference),                    "ggplot")
})

test_that("plot.grass_metrics errors with a redirect message", {
  skip_if_not_installed("ggplot2")
  m <- grass_compute(fixture_cohen_1960, format = "matrix")
  expect_error(plot(m), "grass_report")
})

test_that("theme_grass returns a ggplot2 theme", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(theme_grass(), "theme")
})
