source(test_path("fixtures", "published-tables.R"))

# ---- response = alias ------------------------------------------------

test_that("`response =` is an alias for `rater_cols =` in wide format", {
  set.seed(1)
  df <- data.frame(
    xray = 1:20,
    r1   = rbinom(20, 1, 0.3),
    r2   = rbinom(20, 1, 0.3)
  )
  a <- suppressWarnings(
    grass_compute(df, format = "wide", id_col = "xray",
                  rater_cols = c("r1", "r2"))
  )
  b <- suppressWarnings(
    grass_compute(df, format = "wide", id_col = "xray",
                  response = c("r1", "r2"))
  )
  expect_equal(a$values, b$values)
})

test_that("`response =` and `rater_cols =` conflict errors", {
  set.seed(2)
  df <- data.frame(
    r1 = rbinom(20, 1, 0.3),
    r2 = rbinom(20, 1, 0.3),
    other = rbinom(20, 1, 0.3)
  )
  expect_error(
    grass_compute(df, format = "wide",
                  rater_cols = c("r1", "r2"),
                  response = c("r1", "other")),
    "disagree"
  )
})

test_that("`response =` works in long format as alias for `rating =`", {
  long_df <- data.frame(
    subject = rep(1:10, 2),
    rater   = rep(c("r1", "r2"), each = 10),
    score   = rbinom(20, 1, 0.3)
  )
  a <- suppressWarnings(grass_compute(long_df, format = "long", rating = "score"))
  b <- suppressWarnings(grass_compute(long_df, format = "long", response = "score"))
  expect_equal(a$values, b$values)
})

# ---- ci_width on grass_format_report ---------------------------------

test_that("ci_width = TRUE appends width and descriptor", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  s <- grass_format_report(r, ci_width = TRUE)
  expect_true(grepl("CI width =", s))
  expect_true(grepl("tight|moderate|wide", s))
})

test_that("ci_width = FALSE (default) does not append width", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  s <- grass_format_report(r)
  expect_false(grepl("CI width", s))
})

# ---- tidy.grass_result long form -------------------------------------

test_that("tidy.grass_result returns 9 rows with expected schema", {
  r <- grass_report(fixture_cohen_1960, format = "matrix")
  tdf <- grass:::tidy.grass_result(r)
  expect_s3_class(tdf, "data.frame")
  expect_equal(nrow(tdf), 9)
  expect_setequal(tdf$metric, c("kappa", "PABAK", "AC1"))
  expect_setequal(tdf$quantity, c("estimate", "reference", "distance"))
  expect_true(all(c("metric", "quantity", "value", "conf.low",
                    "conf.high", "n", "prevalence") %in% names(tdf)))
  # conf.low / conf.high are populated only for kappa estimate
  kappa_est <- tdf[tdf$metric == "kappa" & tdf$quantity == "estimate", ]
  expect_false(is.na(kappa_est$conf.low))
  pabak_ref <- tdf[tdf$metric == "PABAK" & tdf$quantity == "reference", ]
  expect_true(is.na(pabak_ref$conf.low))
})

# ---- grass_report_by -------------------------------------------------

test_that("grass_report_by returns one row per group with .cohort column", {
  set.seed(3)
  df <- data.frame(
    xray_id = sprintf("CXR-%03d", 1:60),
    rater_A = sample(c("abnormal", "normal"), 60, replace = TRUE),
    rater_B = sample(c("abnormal", "normal"), 60, replace = TRUE),
    cohort  = rep(c("site1", "site2", "site3"), each = 20)
  )
  out <- suppressWarnings(
    grass_report_by(df, cohort, id_col = "xray_id", positive = "abnormal")
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3)
  expect_true(".cohort" %in% names(out))
  expect_setequal(out$.cohort, c("site1", "site2", "site3"))
  # regime_note should be absent (compact = TRUE inside)
  expect_false("regime_note" %in% names(out))
})

test_that("grass_report_by accepts a string for `group`", {
  set.seed(4)
  df <- data.frame(
    r1 = rbinom(30, 1, 0.3),
    r2 = rbinom(30, 1, 0.3),
    cohort = rep(c("A", "B", "C"), each = 10)
  )
  out <- suppressWarnings(grass_report_by(df, "cohort"))
  expect_equal(nrow(out), 3)
  expect_setequal(out$.cohort, c("A", "B", "C"))
})

# ---- id_col hint in wide-format error --------------------------------

test_that("wide-format error softly suggests id_col when one column is non-binary", {
  df <- data.frame(
    subject_id = sprintf("S-%03d", 1:10),
    a = c(1,0,1,0,1,1,0,1,0,1),
    b = c(0,1,1,0,1,0,1,1,0,0)
  )
  msg <- tryCatch(grass_compute(df, format = "wide"),
                  error = function(e) conditionMessage(e))
  expect_true(grepl("subject_id", msg))
  expect_true(grepl("looks like an identifier", msg))
})

test_that("wide-format error hint is silent when ambiguous (no clear id column)", {
  df <- data.frame(
    a = c(1,0,1,0,1),
    b = c(0,1,1,0,0),
    c = c(1,1,0,0,1)
  )
  msg <- tryCatch(grass_compute(df, format = "wide"),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("looks like an identifier", msg))
})

# ---- reset_grass_warnings --------------------------------------------

test_that("reset_grass_warnings clears the once-per-session cache", {
  env <- grass:::.grass_env
  env$msg_seen[["test_key"]] <- TRUE
  expect_true(isTRUE(env$msg_seen[["test_key"]]))
  reset_grass_warnings()
  expect_null(env$msg_seen[["test_key"]])
})
