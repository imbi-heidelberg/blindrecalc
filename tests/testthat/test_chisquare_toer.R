context("test toer for the ChiSquare-test")

test_that("toer gives same values as pow for delta = 0", {
  design1 <- setupChiSquare(alpha = 0.05, beta = 0.2, r = 1,
    delta = 0)
  size1 <- toer(design1, n1 = 100, nuisance = 0.2, recalculation = FALSE)
  pow1 <- pow(design1, n1 = 100, nuisance = 0.2, recalculation = FALSE)

  design2 <- setupChiSquare(alpha = 0.05, beta = 0.2, r = 1,
    delta = 0, n_max = 150)
  size2 <- toer(design2, n1 = 75, nuisance = 0.4, recalculation = TRUE,
    allocation = "approximate")
  pow2 <- toer(design2, n1 = 75, nuisance = 0.4, recalculation = TRUE,
    allocation = "approximate")

  expect_equal(size1, pow1)
  skip_on_cran()
  expect_equal(size2, pow2)
})
