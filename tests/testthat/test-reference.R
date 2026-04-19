test_that("grass_reference returns a 3-row data.frame", {
  r <- grass_reference(0.5)
  expect_s3_class(r, "grass_reference")
  expect_equal(nrow(r$reference), 3)
  expect_setequal(r$reference$metric, c("kappa", "PABAK", "AC1"))
  expect_equal(r$prevalence, 0.5)
  expect_equal(r$quality, "high")
})

test_that("grass_reference interpolates between stored points", {
  t1 <- grass_reference(0.10)
  t2 <- grass_reference(0.12)
  t3 <- grass_reference(0.15)
  k1 <- t1$reference$reference[t1$reference$metric == "kappa"]
  k2 <- t2$reference$reference[t2$reference$metric == "kappa"]
  k3 <- t3$reference$reference[t3$reference$metric == "kappa"]
  expect_true(k2 >= min(k1, k3) - 1e-8 && k2 <= max(k1, k3) + 1e-8)
})

test_that("medium quality returns lower reference values", {
  t_hi  <- grass_reference(0.5, quality = "high")
  t_med <- grass_reference(0.5, quality = "medium")
  for (m in c("kappa", "PABAK", "AC1")) {
    th_hi  <- t_hi$reference$reference[t_hi$reference$metric == m]
    th_med <- t_med$reference$reference[t_med$reference$metric == m]
    expect_true(th_med < th_hi, info = m)
  }
})

test_that("Out-of-range prevalence clamps with warning; invalid errors", {
  expect_warning(grass_reference(0.001), "clamping")
  expect_warning(grass_reference(0.995), "clamping")
  expect_error(grass_reference(-0.1), "lie in")
  expect_error(grass_reference(1.5), "lie in")
  expect_error(grass_reference(NA), "finite")
  expect_error(grass_reference(c(0.2, 0.3)), "single")
})

test_that("grass_reference_table returns the full long-form table (all four bands)", {
  tbl <- grass_reference_table()
  expect_s3_class(tbl, "data.frame")
  # 21 prevalence points x 4 reference_levels x 3 metrics = 252. The
  # stored grid spans [0.01, 0.99] in 0.05 steps with end-point ticks at
  # 0.01/0.05 and 0.95/0.99.
  expect_equal(nrow(tbl), 21 * 4 * 3)
  expect_true(all(c("prevalence", "reference_level", "metric", "reference", "J")
                  %in% names(tbl)))
  expect_setequal(unique(tbl$reference_level), c(0.70, 0.80, 0.85, 0.90))
  expect_setequal(unique(tbl$metric), c("kappa", "PABAK", "AC1"))
})

test_that("grass_reference_table(reference_level=) filters to one band", {
  one <- grass_reference_table(reference_level = 0.80)
  expect_equal(nrow(one), 21 * 3)
  expect_equal(unique(one$reference_level), 0.80)
})

test_that("new reference_level bands are monotonically ordered", {
  # At a mid prevalence, higher Se=Sp should give a higher reference
  # value for each metric.
  mid <- 0.5
  vals <- lapply(c(0.70, 0.80, 0.85, 0.90), function(q) {
    r <- grass_reference(mid, reference_level = q)$reference
    setNames(r$reference, r$metric)
  })
  for (m in c("kappa", "PABAK", "AC1")) {
    seq_vals <- vapply(vals, function(v) v[[m]], numeric(1))
    expect_true(all(diff(seq_vals) > 0), info = m)
  }
})
