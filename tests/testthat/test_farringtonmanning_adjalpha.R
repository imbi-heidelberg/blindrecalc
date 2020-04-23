context("test adjusted alpha for Farrington-Manning test")

test_that("errors are defined correctly", {
  d1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1)
  expect_error(adjusted_alpha(d1, n1 = 30, nuisance = 1.01, TRUE))

  d2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1, n_max = 301)
  expect_error(adjusted_alpha(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(adjusted_alpha(d2, 20, 0.5, TRUE, "exact"))

})


test_that("adjusted alpha works for recalculation", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1,
                              delta = 0, delta_NI = 0.1, n_max = 300)

  expect_lte(adjusted_alpha(d, n1 = 20, nuisance = 0.5, recalculation = TRUE), d@alpha)

})


test_that("adjusted alpha works without recalculation", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1,
                              delta = 0, delta_NI = 0.1, n_max = 300)
  n <- n_fix(d, 0.5)
  expect_lte(
    adjusted_alpha(d, n1 = n/2, nuisance = 0.5, nuis_ass = 0.5, recalculation = FALSE, allocation = "exact"),
    d@alpha
  )


  expect_lte(
    adjusted_alpha(d, n1 = n/2, nuisance = 0.5, nuis_ass = 0.5, recalculation = FALSE, allocation = "approximate"),
    d@alpha
  )

})

