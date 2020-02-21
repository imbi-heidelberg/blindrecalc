context("test toer for the Farrington-Manning test")

test_that("toer gives same values as pow for delta = 0", {
  design1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1,
    delta = 0, delta_NI = 0, n_max = 250)
  toer1 <- toer(design1, n1 = 80, nuisance = 0.4, recalculation = FALSE)
  power1 <- pow(design1, n1 = 80, nuisance = 0.4, recalculation = FALSE)

  expect_equal(toer1, power1)
})
