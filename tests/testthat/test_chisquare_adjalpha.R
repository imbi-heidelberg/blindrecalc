context("test adjusted_alpha for the Chi-Squared test")

test_that("adjusted_alpha gives actual level at most nominal level", {
  nuis_vec <- seq(0.2, 0.8, by = 0.1)
  gamma <- 0.001

  design1   <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
  n1        <- n_fix(design1, nuisance = 0.3)
  adjalpha1 <- adjusted_alpha(design1, n1 = n1, nuisance = nuis_vec,
                              gamma = gamma, nuis_ass = 0.3, precision = 0.001,
                              recalculation = FALSE)
  design1@alpha <- adjalpha1
  n1_new        <- n_fix(design1, nuisance = 0.3)
  alpha_max1    <- max(toer(design1, n1 = n1_new, nuisance = nuis_vec,
                         recalculation = FALSE, allocation = "exact"))

  design2   <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
  n2        <- n_fix(design2, nuisance = 0.3)
  adjalpha2 <- adjusted_alpha(design2, n1 = n2, nuisance = nuis_vec,
                              gamma = gamma, nuis_ass = 0.3, precision = 0.001,
                              recalculation = FALSE, allocation = "approximate")
  design2@alpha <- adjalpha2
  n2_new        <- n_fix(design2, nuisance = 0.3)
  alpha_max2 <- max(toer(design2, n1 = n2_new, nuisance = nuis_vec,
                         recalculation = FALSE, allocation = "approximate"))


  design3 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 3, delta = 0.3)
  adjalpha3 <- adjusted_alpha(design3, n1 = 50, nuisance = nuis_vec,
                              gamma = gamma, precision = 0.005,
                              recalculation = TRUE, allocation = "approximate")
  design3@alpha <- adjalpha3
  alpha_max3 <- max(toer(design3, n1 = 50, nuisance = nuis_vec,
                         recalculation = TRUE, allocation = "approximate"))

  expect_lte(alpha_max1, 0.025 - gamma)
  expect_lte(alpha_max2, 0.025 - gamma)
  expect_lte(alpha_max3, 0.025 - gamma)
})



test_that("errors are thrown correctly", {
  d1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 2, delta = 0.1, n_max = 300)
  expect_error(adjusted_alpha(d1, n1 = 21, nuisance = 1.1, TRUE))
  expect_error(adjusted_alpha(d1, n1 = 20, nuisance = 0.8, TRUE))
  d1@n_max <- 301
  expect_error(adjusted_alpha(d1, n1 = 21, nuisance = 0.8, TRUE))
})


