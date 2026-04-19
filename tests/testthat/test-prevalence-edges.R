test_that("validate_prevalence errors on exact 0 and 1 (degenerate)", {
  expect_error(grass_reference(0), "degenerate")
  expect_error(grass_reference(1), "degenerate")
})

test_that("validate_prevalence warn-clamps small and large (non-degenerate)", {
  expect_warning(grass_reference(0.005), "clamping to 0.01")
  expect_warning(grass_reference(0.995), "clamping to 0.99")
})

test_that("validate_prevalence accepts values strictly inside (0, 1)", {
  expect_silent(grass_reference(0.5))
  expect_silent(grass_reference(0.01))
  expect_silent(grass_reference(0.99))
})
