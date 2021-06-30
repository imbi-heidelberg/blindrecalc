context("test pow for the ChiSquare-test")

test_that("pow works for fixed design", {
  # Compare with Kieser 2020, Table 5.2
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
  pow1 <- pow(design1, n1 = 164, nuisance = 0.3, recalculation = FALSE)
  pow2 <- pow(design1, n1 = 186, nuisance = 0.4, recalculation = FALSE)

  design2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 2,
    delta = 0.2)
  pow3 <- pow(design2, n1 = 207, nuisance = (19 / 30),
    recalculation = FALSE)
  pow4 <- pow(design2, n1 = 189, nuisance = (1 / 3),
    recalculation = FALSE)

  # Compare with PASS 16.0.3 (Z-Test Pooled)
  design4 <- setupChiSquare(alpha = 0.025, beta = 0.5, r = 1, delta = 0.2)
  pow5 <- pow(design4, n1 = 62, nuisance = 0.2, allocation = "exact",
    recalculation = FALSE)
  pow6 <- pow(design4, n1 = 78, nuisance = 0.3, allocation = "exact",
    recalculation = FALSE)
  pow7 <- pow(design4, n1 = 94, nuisance = 0.4, allocation = "exact",
    recalculation = FALSE)
  pow8 <- pow(design4, n1 = 96, nuisance = 0.5, allocation = "exact",
    recalculation = FALSE)

  design5 <- setupChiSquare(alpha = 0.025, beta = 0.8, r = 1, delta = 0.2)
  pow9 <- pow(design5, n1 = 26, nuisance = 0.2, allocation = "exact",
    recalculation = FALSE)
  pow10 <- pow(design5, n1 = 26, nuisance = 0.3, allocation = "exact",
    recalculation = FALSE)
  pow11 <- pow(design5, n1 = 26, nuisance = 0.4, allocation = "exact",
    recalculation = FALSE)
  pow12 <- pow(design5, n1 = 26, nuisance = 0.5, allocation = "exact",
    recalculation = FALSE)

  expect_equal(pow1, 0.8074344, tolerance = 1e-7)
  expect_equal(pow2, 0.7991127, tolerance = 1e-7)
  expect_equal(pow3, 0.8028214, tolerance = 1e-7)
  expect_equal(pow4, 0.8114369, tolerance = 1e-7)
  expect_equal(pow5, 0.5159118, tolerance = 1e-7)
  expect_equal(pow6, 0.5032302, tolerance = 1e-7)
  expect_equal(pow7, 0.5041563, tolerance = 1e-7)
  expect_equal(pow8, 0.5111488, tolerance = 1e-7)
  expect_equal(pow9, 0.2252016, tolerance = 1e-7)
  expect_equal(pow10, 0.2125732, tolerance = 1e-7)
  expect_equal(pow11, 0.2210516, tolerance = 1e-7)
  expect_equal(pow12, 0.225603, tolerance = 1e-7)
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

test_that("power of the recalculation design is close to target power", {
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2,
    r = 1, delta = 0.3)
  pow1 <- pow(design1, n1 = 40, nuisance = c(0.2, 0.3, 0.4, 0.5),
    recalculation = TRUE)

  design2 <- setupChiSquare(alpha = 0.025, beta = 0.4,
    r = 1, delta = 0.2)
  pow2 <- pow(design2, n1 = 40, nuisance = c(0.2, 0.3, 0.4, 0.5),
    recalculation = TRUE)

  design3 <- setupChiSquare(alpha = 0.025, beta = 0.6,
    r = 1, delta = 0.2)
  pow3 <- pow(design3, n1 = 30, nuisance = c(0.2, 0.3, 0.4, 0.5),
    recalculation = TRUE)

  expect_equal(pow1, rep(0.8, 4), tolerance = 0.025)
  expect_equal(pow2, rep(0.6, 4), tolerance = 0.025)
  expect_equal(pow3, rep(0.4, 4), tolerance = 0.025)
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
  d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.25)
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
  d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.25)
  expect_equal(
    pow(d, n1 = 20, nuisance = c(0.4, 0.6), recalculation = TRUE, allocation = "approximate"),
    sapply(c(0.4, 0.6), function(p) {
      pow(d, n1 = 20, nuisance = p, recalculation = TRUE, allocation = "approximate")
    })
  )

  expect_equal(
    pow(d, n1 = 20, nuisance = c(0.4, 0.6), recalculation = FALSE, allocation = "exact"),
    sapply(c(0.4, 0.6), function(p) {
      pow(d, n1 = 20, nuisance = p, recalculation = FALSE, allocation = "exact")
    })
  )

})

