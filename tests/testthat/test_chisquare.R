context("test ChiSquare-test")

test_that("sample size is calculated correctly", {
  design <- setupChiSquare(alpha = 0.05, beta = 0.2, r = 1,
    delta = 0.2)

  expect_equal(

  ss1 <- ceiling(getn_chisq(p = 0.2, delta = 0.2, r = 1,
    alpha = 0.05, power = 0.8, variance = "heterogeneous"))

  ss2 <- ceiling(getn_chisq(p = (7 / 30), delta = 0.2, r = 2,
    alpha = 0.05, power = 0.8, variance = "heterogeneous"))

  ss3 <- getn_chisq(p = 0.1, delta = 0.3, r = 1, alpha = 0.05,
    power = 0.8, variance = "heterogeneous")

  expect_equal(ss1, 124)
  expect_equal(ss2, 145)
  expect_equal(ss3, NA)
})
