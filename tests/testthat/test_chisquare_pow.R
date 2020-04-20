context("test pow for the ChiSquare-test")

test_that("pow works for fixed design (Kieser 2020, Table 5.2)", {
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1,
    delta = 0.2)
  pow1 <- pow(design1, n1 = 164, nuisance = 0.3, recalculation = FALSE)
  pow2 <- pow(design1, n1 = 186, nuisance = 0.4, recalculation = FALSE)

  design2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 2,
    delta = 0.2)
  pow3 <- pow(design2, n1 = 207, nuisance = (19 / 30),
    recalculation = FALSE)
  pow4 <- pow(design2, n1 = 189, nuisance = (1 / 3),
    recalculation = FALSE)

  expect_equal(pow1, 0.8074344, tolerance = 1e-7)
  expect_equal(pow2, 0.7991127, tolerance = 1e-7)
  expect_equal(pow3, 0.8028214, tolerance = 1e-7)
  expect_equal(pow4, 0.8114369, tolerance = 1e-7)
})

test_that("pow works for recalculation design (Friede & Kieser 2004)", {
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1,
    delta = 0.2)
  pow1 <- pow(design1, n1 = 80, nuisance = 0.3, recalculation = TRUE,
    allocation = "kf_approx", variance = "homogeneous")
  pow2 <- pow(design1, n1 = 120, nuisance = 0.3, recalculation = TRUE,
    allocation = "kf_approx", variance = "homogeneous")

  design2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 3,
    delta = 0.2)
  pow3 <- pow(design2, n1 = 80, nuisance = 0.75, recalculation = TRUE,
    allocation = "kf_approx", variance = "homogeneous")
  pow4 <- pow(design2, n1 = 120, nuisance = 0.75, recalculation = TRUE,
    allocation = "kf_approx", variance = "homogeneous")

  skip_on_cran()
  expect_equal(pow1, 0.8020152, tolerance = 1e-7)
  expect_equal(pow2, 0.8055927, tolerance = 1e-7)
  expect_equal(pow3, 0.7749653, tolerance = 1e-7)
  expect_equal(pow4, 0.7789805, tolerance = 1e-7)
})
