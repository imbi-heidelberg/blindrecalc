context("test pow for the Farrington-Manning test")

test_that("pow works for the fixed design (Friede et al. 2007)", {
  design1 <- setupFarringtonManning(alpha = 0.05, beta = 0.2,
    r = 1, delta = 0, delta_NI = 0.1)
  pow1 <- pow(design1, n1 = 658, nuisance = 0.3, recalculation = FALSE)
  pow2 <- pow(design1, n1 = 780, nuisance = 0.5, recalculation = FALSE)

  design2 <- setupFarringtonManning(alpha = 0.05, beta = 0.2,
    r = 2, delta = 0, delta_NI = 0.1)
  n1d2 <- n_fix(design2, nuisance = 0.4, rounded = FALSE)
  pow3 <- pow(design2, n1 = n1d2, nuisance = 0.4,
    allocation = "approximate", recalculation = FALSE)
  n2d2 <- n_fix(design2, nuisance = 0.6, rounded = FALSE)
  pow4 <- pow(design2, n1 = n2d2, nuisance = 0.6,
    allocation = "approximate", recalculation = FALSE)

  design3 <- setupFarringtonManning(alpha = 0.05, beta = 0.2,
    r = 3, delta = 0, delta_NI = 0.1)
  n1d3 <- n_fix(design3, nuisance = 0.3, rounded = FALSE)
  pow5 <- pow(design3, n1 = n1d3, nuisance = 0.3,
    allocation = "approximate", recalculation = FALSE)
  n2d3 <- n_fix(design3, nuisance = 0.4, rounded = FALSE)
  pow6 <- pow(design3, n1 = n2d3, nuisance = 0.4,
    allocation = "approximate", recalculation = FALSE)


  expect_equal(pow1, 0.8005762, tolerance = 1e-7)
  expect_equal(pow2, 0.7948942, tolerance = 1e-7)
  expect_equal(pow3, 0.8005885, tolerance = 1e-7)
  expect_equal(pow4, 0.8000787, tolerance = 1e-7)
  expect_equal(pow5, 0.8009629, tolerance = 1e-7)
  expect_equal(pow6, 0.8026954, tolerance = 1e-7)
})
