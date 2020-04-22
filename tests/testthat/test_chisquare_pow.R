context("test pow for the ChiSquare-test")

test_that("pow works for fixed design (Kieser 2020, Table 5.2)", {
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
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
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
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

  expect_equal(pow1, 0.8020152, tolerance = 1e-7)
  expect_equal(pow2, 0.8055927, tolerance = 1e-7)
  expect_equal(pow3, 0.7749653, tolerance = 1e-7)
  expect_equal(pow4, 0.7789805, tolerance = 1e-7)
})




test_that("errors are thrown correctly", {
  d1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.1)
  expect_error(pow(d1, n1 = 20, nuisance = 1.1, TRUE))

  d2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 2, delta = 0.1, n_max = 301)
  expect_error(pow(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(pow(d2, 20, 0.5, TRUE, "exact"))

  expect_error(pow(d2, d2@n_max + 1, 0.7, TRUE, "approximate"))

  expect_error(pow(d2, c(20, 30), c(0.6, 0.7), TRUE, "approximate"))

})


test_that("vectorization in n1 works", {
  skip_on_cran()
  d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.1)
  expect_equal(
    pow(d, n1 = c(10, 20), nuisance = 0.25, recalculation = TRUE, allocation = "approximate"),
    sapply(c(10, 20), function(n) {
      pow(d, n1 = n, nuisance = 0.25, recalculation = TRUE, allocation = "approximate")
    })
  )

  expect_equal(
    pow(d, n1 = c(10, 20), nuisance = 0.25, recalculation = FALSE, allocation = "approximate"),
    sapply(c(10, 20), function(n) {
      pow(d, n1 = n, nuisance = 0.25, recalculation = FALSE, allocation = "approximate")
    })
  )

})


test_that("vectorization in nuisance works", {
  skip_on_cran()
  d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.1)
  expect_equal(
    pow(d, n1 = 15, nuisance = c(0.4, 0.6), recalculation = TRUE, allocation = "approximate"),
    sapply(c(0.4, 0.6), function(p) {
      pow(d, n1 = 15, nuisance = p, recalculation = TRUE, allocation = "approximate")
    })
  )

  expect_equal(
    pow(d, n1 = 15, nuisance = c(0.4, 0.6), recalculation = FALSE, allocation = "exact"),
    sapply(c(0.4, 0.6), function(p) {
      pow(d, n1 = 15, nuisance = p, recalculation = FALSE, allocation = "exact")
    })
  )

})

