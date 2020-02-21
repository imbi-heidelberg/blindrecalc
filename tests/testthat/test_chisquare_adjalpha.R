context("test adjusted_alpha for the Chi-Squared test")

test_that("adjusted_alpha gives actual level at most nominal level", {
  nuis_vec <- seq(0.2, 0.8, by = 0.1)
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1,
    delta = 0.1)
  adjalpha1 <- adjusted_alpha(design1, n1 = 100, nuisance = nuis_vec,
    precision = 0.001, recalculation = FALSE)
  design1@alpha <- adjalpha1
  alpha_max1 <- max(toer(design1, n1 = 100, nuisance = nuis_vec,
    recalculation = FALSE))

  design2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 3,
    delta = 0.3)
  adjalpha2 <- adjusted_alpha(design2, n1 = 50, nuisance = nuis_vec,
    precision = 0.005, recalculation = TRUE, allocation = "approximate")
  design2@alpha <- adjalpha2
  alpha_max2 <- max(toer(design2, n1 = 50, nuisance = nuis_vec,
    recalculation = TRUE, allocation = "approximate"))

  expect_lte(alpha_max1, 0.025)
  expect_lte(alpha_max2, 0.025)
})
