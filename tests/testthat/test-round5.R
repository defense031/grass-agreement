source(test_path("fixtures", "published-tables.R"))

# ---- spec-dispatch architecture ---------------------------------------

test_that("grass_spec_binary constructs with valid reference_level", {
  s <- grass_spec_binary()
  expect_s3_class(s, "grass_spec_binary")
  expect_s3_class(s, "grass_spec")
  expect_equal(s$family, "binary")
  expect_equal(s$reference_level, 0.85)
  for (lvl in c(0.70, 0.80, 0.85, 0.90)) {
    expect_equal(grass_spec_binary(reference_level = lvl)$reference_level, lvl)
  }
  expect_null(grass_spec_binary(reference_level = NULL)$reference_level)
})

test_that("grass_spec_binary rejects invalid reference_level", {
  expect_error(grass_spec_binary(reference_level = 0.75), "must be one of")
  expect_error(grass_spec_binary(reference_level = "high"), "must be one of")
  expect_error(grass_spec_binary(reference_level = c(0.85, 0.90)), "must be one of")
})

test_that("stub spec constructors return placeholder specs", {
  for (ctor in list(grass_spec_ordinal, grass_spec_multirater, grass_spec_continuous)) {
    s <- ctor()
    expect_s3_class(s, "grass_spec")
  }
})

test_that("stub specs error at dispatch with a ?grass_roadmap pointer", {
  tab <- matrix(c(88, 10, 14, 88), nrow = 2,
                dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
  for (spec in list(grass_spec_ordinal(), grass_spec_multirater(),
                    grass_spec_continuous())) {
    expect_error(grass_report(tab, format = "matrix", spec = spec),
                 "grass_roadmap")
  }
})

test_that("grass_report returns spec slot on result", {
  tab <- matrix(c(88, 10, 14, 88), nrow = 2,
                dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
  r <- grass_report(tab, format = "matrix",
                    spec = grass_spec_binary(reference_level = 0.80))
  expect_s3_class(r$spec, "grass_spec_binary")
  expect_equal(r$spec$reference_level, 0.80)
})

test_that("reference_level = NULL skips the reference attachment", {
  tab <- matrix(c(88, 10, 14, 88), nrow = 2,
                dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
  r <- grass_report(tab, format = "matrix",
                    spec = grass_spec_binary(reference_level = NULL))
  expect_null(r$reference)
  expect_null(r$distance)
})

# ---- reference_level bands --------------------------------------------

test_that("Each reference_level matches the analytical diagonal closed-form", {
  # At p=0.5, Se=Sp=q gives P_agree = q^2 + (1-q)^2 and Pe_cohen = 0.5.
  # kappa = 2*(q^2 + (1-q)^2 - 0.5) = 2q^2 - 2q + ... simplifies to
  # (q-0.5)^2 * 4 basically; cross-check numerically.
  expected <- function(p, q) {
    m <- q*p + (1-q)*(1-p)
    P_agree <- p*(q^2 + (1-q)^2) + (1-p)*(q^2 + (1-q)^2)
    Pe <- m^2 + (1-m)^2
    (P_agree - Pe) / (1 - Pe)
  }
  for (q in c(0.70, 0.80, 0.85, 0.90)) {
    r <- grass_reference(0.5, reference_level = q)
    k <- r$reference$reference[r$reference$metric == "kappa"]
    expect_equal(k, expected(0.5, q), tolerance = 1e-6,
                 info = sprintf("kappa at q=%.2f", q))
  }
})

# ---- .parallel on grass_report_by -------------------------------------

test_that(".parallel matches sequential output", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("progressr")
  set.seed(11)
  df <- data.frame(
    xray = sprintf("CXR-%03d", 1:60),
    rA   = sample(c("abnormal","normal"), 60, replace = TRUE),
    rB   = sample(c("abnormal","normal"), 60, replace = TRUE),
    cohort = rep(c("s1","s2","s3"), each = 20)
  )
  r_seq <- grass_report_by(df, cohort, id_col = "xray", positive = "abnormal")
  future::plan(future::sequential)
  progressr::handlers("void")
  r_par <- progressr::with_progress(
    grass_report_by(df, cohort, id_col = "xray", positive = "abnormal",
                    .parallel = TRUE)
  )
  for (col in c("N","kappa","PABAK","AC1","prevalence")) {
    expect_equal(r_seq[[col]], r_par[[col]], info = col)
  }
})

test_that(".parallel errors cleanly when future.apply is unavailable", {
  # Simulate missing future.apply via namespace unloading check.
  has_future <- requireNamespace("future.apply", quietly = TRUE)
  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  skip_if(has_future && has_progressr,
          "test requires at least one of future.apply / progressr to be absent")
  # When at least one of them is absent, .parallel = TRUE must error.
  set.seed(12)
  df <- data.frame(
    rA = rbinom(20, 1, 0.3),
    rB = rbinom(20, 1, 0.3),
    cohort = rep(c("s1","s2"), each = 10)
  )
  expect_error(grass_report_by(df, cohort, .parallel = TRUE), "install.packages")
})

# ---- grass_methods() --------------------------------------------------

test_that("grass_methods returns a non-trivial paragraph", {
  tab <- matrix(c(88, 10, 14, 88), nrow = 2,
                dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
  r <- grass_report(tab, format = "matrix")
  out <- grass_methods(r)
  expect_type(out, "character")
  expect_length(out, 1)
  expect_match(out, "GRASS framework", fixed = TRUE)
  expect_match(out, "N = 200", fixed = TRUE)
  expect_match(out, "Se = Sp = 0.85", fixed = TRUE)
})

test_that("grass_methods honours all three formats", {
  tab <- matrix(c(88, 10, 14, 88), nrow = 2,
                dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
  r <- grass_report(tab, format = "matrix")
  md <- grass_methods(r, format = "markdown")
  tex <- grass_methods(r, format = "latex")
  plain <- grass_methods(r, format = "plain")
  expect_match(md, "\u03ba", fixed = TRUE)          # Unicode kappa
  expect_match(tex, "\\kappa", fixed = TRUE)        # LaTeX macro
  expect_match(plain, "kappa", fixed = TRUE)         # word
  expect_match(md, "*balanced*", fixed = TRUE)
  expect_match(tex, "\\emph{balanced}", fixed = TRUE)
})

test_that("grass_methods describes reference_level = NULL gracefully", {
  tab <- matrix(c(88, 10, 14, 88), nrow = 2,
                dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
  r <- grass_report(tab, format = "matrix",
                    spec = grass_spec_binary(reference_level = NULL))
  out <- grass_methods(r)
  expect_match(out, "No prevalence-conditioned reference", fixed = TRUE)
})

# ---- legacy reference = ... still works with deprecation --------------

test_that("legacy reference = 'high' maps to reference_level = 0.85", {
  tab <- matrix(c(88, 10, 14, 88), nrow = 2,
                dimnames = list(R1 = c("0","1"), R2 = c("0","1")))
  reset_grass_warnings()
  expect_message(
    grass_report(tab, format = "matrix", reference = "high"),
    "deprecated"
  )
  reset_grass_warnings()
  r <- suppressMessages(grass_report(tab, format = "matrix", reference = "medium"))
  expect_equal(r$spec$reference_level, 0.70)
})
