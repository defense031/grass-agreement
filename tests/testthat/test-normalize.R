source(test_path("fixtures", "published-tables.R"))

# Helper: from a 2x2 count table, reconstruct the implied rater vectors in
# each of the four supported shapes.
reconstruct_from_table <- function(tab) {
  a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
  r1 <- c(rep(0L, a + b), rep(1L, c + d))
  r2 <- c(rep(0L, a), rep(1L, b), rep(0L, c), rep(1L, d))
  list(r1 = r1, r2 = r2)
}

test_that("coerce_binary handles logical, 0/1 numeric, factor, character", {
  out <- coerce_binary(c(TRUE, FALSE, TRUE, NA))
  expect_equal(out$values, c(1L, 0L, 1L, NA_integer_))

  out <- coerce_binary(c(0, 1, 1, 0))
  expect_equal(out$values, c(0L, 1L, 1L, 0L))

  out <- coerce_binary(factor(c("yes", "no", "yes", "no")))
  expect_equal(out$values, c(1L, 0L, 1L, 0L))
  expect_equal(out$positive_level, "yes")

  out <- coerce_binary(c("case", "control", "case"))
  expect_equal(out$positive_level, "case")

  out <- coerce_binary(c("A", "B", "A", "B"), positive = "B")
  expect_equal(out$positive_level, "B")
  expect_equal(out$values, c(0L, 1L, 0L, 1L))
})

test_that("coerce_binary rejects >2 levels and numeric non-binary", {
  expect_error(coerce_binary(factor(c("a", "b", "c"))), "binary")
  expect_error(coerce_binary(c(0, 1, 2)), "0/1")
})

test_that("normalize_matrix rejects non-2x2", {
  expect_error(normalize_input(matrix(1:6, nrow = 2), format = "matrix"), "2x2")
})

test_that("Cross-format invariance: all four formats yield identical metrics", {
  skip_on_cran()
  for (nm in names(all_fixtures)) {
    tab <- all_fixtures[[nm]]
    rv <- reconstruct_from_table(tab)

    m_matrix <- grass_compute(tab, format = "matrix")

    m_paired <- grass_compute(list(rv$r1, rv$r2), format = "paired")

    wide_df <- data.frame(rater1 = rv$r1, rater2 = rv$r2)
    m_wide <- grass_compute(wide_df, format = "wide")

    N <- length(rv$r1)
    long_df <- data.frame(
      subject = rep(seq_len(N), 2),
      rater   = rep(c("rater1", "rater2"), each = N),
      rating  = c(rv$r1, rv$r2)
    )
    m_long <- grass_compute(long_df, format = "long")

    core <- c("N", "P0", "kappa", "PABAK", "AC1",
              "prevalence_index", "bias_index")
    expect_equal(m_matrix$values[core], m_paired$values[core], info = paste(nm, "paired"))
    expect_equal(m_matrix$values[core], m_wide$values[core],   info = paste(nm, "wide"))
    expect_equal(m_matrix$values[core], m_long$values[core],   info = paste(nm, "long"))
  }
})

test_that("Long format errors on wrong rater count", {
  df <- data.frame(subject = 1:3,
                   rater = c("A", "A", "A"),
                   rating = c(1, 0, 1))
  expect_error(grass_compute(df, format = "long"), "2 distinct raters")
})

test_that("Long format errors on duplicate subject x rater", {
  df <- data.frame(subject = c(1, 1, 2, 2, 1),
                   rater = c("A", "B", "A", "B", "A"),
                   rating = c(1, 0, 1, 1, 0))
  expect_error(grass_compute(df, format = "long"), "duplicated")
})

test_that("NA pairs are dropped with a warning", {
  r1 <- c(1, 0, NA, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
  r2 <- c(1, 0, 1, NA, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1)
  expect_warning(res <- grass_compute(list(r1, r2), format = "paired"), "Dropping")
  expect_equal(res$n, 13)
  expect_equal(res$n_dropped, 2)
})

test_that(">50% NA pairs errors", {
  r1 <- c(NA, NA, NA, 1)
  r2 <- c(0, 1, 1, 1)
  expect_error(grass_compute(list(r1, r2), format = "paired"), "50%")
})

test_that("Wide format with extra columns requires rater_cols", {
  set.seed(1)
  df <- data.frame(
    id = 1:20,
    rater1 = rbinom(20, 1, 0.5),
    rater2 = rbinom(20, 1, 0.5),
    other  = letters[1:20]
  )
  expect_error(grass_compute(df, format = "wide"), "exactly two rater columns")
  res <- grass_compute(df, format = "wide", rater_cols = c("rater1", "rater2"))
  expect_s3_class(res, "grass_metrics")
  # id_col dropping: remove `id`, detect remaining two raters + `other`.
  res2 <- grass_compute(df, format = "wide", id_col = "id",
                        rater_cols = c("rater1", "rater2"))
  expect_s3_class(res2, "grass_metrics")
})

test_that("Small N triggers warning", {
  r1 <- c(1, 0, 1, 0, 1)
  r2 <- c(1, 0, 0, 0, 1)
  expect_warning(grass_compute(list(r1, r2), format = "paired"), "Small sample")
})
