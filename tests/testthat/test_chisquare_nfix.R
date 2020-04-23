context("test n_fix for the ChiSquare-test")

test_that("n_fix works with variance = heterogeneous (Kieser 2020, Table 5.2)", {
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
  ss1 <- n_fix(design1, nuisance = c(0.2, 0.3))

  design2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 2, delta = 0.2)
  ss2 <- n_fix(design2, nuisance = c((7 / 30), (1 / 3)))

  expect_equal(ss1, c(124, 164))
  expect_equal(ss2, c(147, 189))
})

test_that("n_fix works with variance = homogeneous (Friede & Kieser 2004)", {
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
  ss1 <- n_fix(design1, nuisance = c(0.15, 0.3), variance = "homogeneous")

  design2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 3, delta = 0.2)
  ss2 <- n_fix(design2, nuisance = c(0.2, 0.55), variance = "homogeneous")

  expect_equal(ss1, c(102, 166))
  expect_equal(ss2, c(168, 260))
})

test_that("n_fix returns NA if delta is unachievable", {
  design_na <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.3)
  ss_na1 <- n_fix(design_na, nuisance = 0.1)
  ss_na2 <- n_fix(design_na, nuisance = 0.9)

  expect_equal(ss_na1, NA)
  expect_equal(ss_na2, NA)

  expect_error(n_fix(design_na, nuisance = 1.2))
})

