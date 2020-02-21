context("test n_fix for the Farrington-Manning test")

test_that("n_fix works (Friede et al. 2007)", {
  design1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
    r = 1, delta = 0, delta_NI = 0.1)
  ss1 <- n_fix(design1, nuisance = 0.3)
  ss2 <- n_fix(design1, nuisance = 0.5)

  design2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
    r = 2, delta = 0, delta_NI = 0.1)
  ss3_temp <- n_fix(design2, nuisance = 0.4, rounded = FALSE)
  ss3 <- ceiling(ss3_temp / (design2@r + 1)) +
    ceiling(design2@r * ss3_temp / (design2@r + 1))
  ss4_temp <- n_fix(design2, nuisance = 0.7, rounded = FALSE)
  ss4 <- ceiling(ss4_temp / (design2@r + 1)) +
    ceiling(design2@r * ss4_temp / (design2@r + 1))

  design3 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
    r = 3, delta = 0, delta_NI = 0.1)
  ss5_temp <- n_fix(design3, nuisance = 0.5, rounded = FALSE)
  ss5 <- ceiling(ss5_temp / (design3@r + 1)) +
    ceiling(design3@r * ss5_temp / (design3@r + 1))
  ss6_temp <- n_fix(design3, nuisance = 0.6, rounded = FALSE)
  ss6 <- ceiling(ss6_temp / (design3@r + 1)) +
    ceiling(design3@r * ss6_temp / (design3@r + 1))

  expect_equal(ss1, 658)
  expect_equal(ss2, 780)
  expect_equal(ss3, 857)
  expect_equal(ss4, 707)
  expect_equal(ss5, 1035)
  expect_equal(ss6, 966)
})
