context("test pow for the Farrington-Manning test")

test_that("pow works for the fixed design", {
  # Compare with Friede et al., 2007
  design1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
                                    r = 1, delta = 0, delta_NI = 0.1)
  pow1 <- pow(design1, n1 = 658, nuisance = 0.3, recalculation = FALSE)
  pow2 <- pow(design1, n1 = 780, nuisance = 0.5, recalculation = FALSE)

  design2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
                                    r = 2, delta = 0, delta_NI = 0.1)
  n1d2 <- n_fix(design2, nuisance = 0.4, rounded = FALSE)
  pow3 <- pow(design2, n1 = n1d2, nuisance = 0.4,
              allocation = "approximate", recalculation = FALSE)
  n2d2 <- n_fix(design2, nuisance = 0.6, rounded = FALSE)
  pow4 <- pow(design2, n1 = n2d2, nuisance = 0.6,
              allocation = "approximate", recalculation = FALSE)

  design3 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
                                    r = 3, delta = 0, delta_NI = 0.1)
  n1d3 <- n_fix(design3, nuisance = 0.3, rounded = FALSE)
  pow5 <- pow(design3, n1 = n1d3, nuisance = 0.3,
              allocation = "approximate", recalculation = FALSE)
  n2d3 <- n_fix(design3, nuisance = 0.4, rounded = FALSE)
  pow6 <- pow(design3, n1 = n2d3, nuisance = 0.4,
              allocation = "approximate", recalculation = FALSE)

  # Compare with PASS 16.0.3
  design4 <- setupFarringtonManning(alpha = 0.025, beta = 0.5,
    r = 1, delta = 0, delta_NI = 0.2)
  pow7 <- pow(design4, n1 = 64, nuisance = 0.2, allocation = "exact",
    recalculation = FALSE)
  pow8 <- pow(design4, n1 = 80, nuisance = 0.3, allocation = "exact",
    recalculation = FALSE)
  pow9 <- pow(design4, n1 = 90, nuisance = 0.4, allocation = "exact",
    recalculation = FALSE)
  pow10 <- pow(design4, n1 = 92, nuisance = 0.5, allocation = "exact",
    recalculation = FALSE)

  design5 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
    r = 1, delta = 0, delta_NI = 0.2)
  pow11 <- pow(design5, n1 = 22, nuisance = 0.2, allocation = "exact",
    recalculation = FALSE)
  pow12 <- pow(design5, n1 = 28, nuisance = 0.3, allocation = "exact",
    recalculation = FALSE)
  pow13 <- pow(design5, n1 = 32, nuisance = 0.4, allocation = "exact",
    recalculation = FALSE)
  pow14 <- pow(design5, n1 = 36, nuisance = 0.5, allocation = "exact",
    recalculation = FALSE)

  expect_equal(pow1, 0.8005762, tolerance = 1e-7)
  expect_equal(pow2, 0.7948942, tolerance = 1e-7)
  expect_equal(pow3, 0.8005885, tolerance = 1e-7)
  expect_equal(pow4, 0.8000787, tolerance = 1e-7)
  expect_equal(pow5, 0.8009629, tolerance = 1e-7)
  expect_equal(pow6, 0.8026954, tolerance = 1e-7)
  expect_equal(pow7, 0.5107838, tolerance = 1e-7)
  expect_equal(pow8, 0.5123605, tolerance = 1e-7)
  expect_equal(pow9, 0.5083387, tolerance = 1e-7)
  expect_equal(pow10, 0.5026442, tolerance = 1e-7)
  expect_equal(pow11, 0.2192882, tolerance = 1e-7)
  expect_equal(pow12, 0.2280977, tolerance = 1e-7)
  expect_equal(pow13, 0.2130656, tolerance = 1e-7)
  expect_equal(pow14, 0.2123263, tolerance = 1e-7)
})

test_that("power of the recalculation design is close to target power", {
  design1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2,
    r = 1, delta = 0, delta_NI = 0.3)
  pow1 <- pow(design1, n1 = 50, nuisance = c(0.2, 0.3, 0.4, 0.5),
    recalculation = TRUE)

  design2 <- setupFarringtonManning(alpha = 0.025, beta = 0.4,
    r = 1, delta = 0, delta_NI = 0.2)
  pow2 <- pow(design2, n1 = 30, nuisance = c(0.2, 0.3, 0.4, 0.5),
    recalculation = TRUE)

  design3 <- setupFarringtonManning(alpha = 0.025, beta = 0.6,
    r = 1, delta = 0, delta_NI = 0.2)
  pow3 <- pow(design3, n1 = 20, nuisance = c(0.2, 0.3, 0.4, 0.5),
    recalculation = TRUE)

  expect_equal(pow1, rep(0.8, 4), tolerance = 0.025)
  expect_equal(pow2, rep(0.6, 4), tolerance = 0.025)
  expect_equal(pow3, rep(0.4, 4), tolerance = 0.025)
})

test_that("errors are defined correctly", {
  d1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1)
  expect_error(pow(d1, n1 = 30, nuisance = 1.1, TRUE))

  d2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1, n_max = 301)
  expect_error(pow(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(pow(d2, 20, 0.5, TRUE, "exact"))

  expect_error(pow(d2, d2@n_max + 1, 0.7, TRUE, "approximate"))

  expect_error(pow(d2, c(20, 30), c(0.6, 0.7), TRUE, "approximate"))

})




test_that("vectorization in n1 works", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.2)
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
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.3)
  n <- round(n_fix(d, 0.4), -1)

  expect_equal(
    pow(d, n1 = n, nuisance = c(0.3, 0.4), recalculation = TRUE, allocation = "approximate"),
    sapply(c(0.3, 0.4), function(p) {
      pow(d, n1 = n, nuisance = p, recalculation = TRUE, allocation = "approximate")
    })
  )

  expect_equal(
    pow(d, n1 = n, nuisance = c(0.3, 0.4), recalculation = FALSE, allocation = "exact"),
    sapply(c(0.3, 0.4), function(p) {
      pow(d, n1 = n, nuisance = p, recalculation = FALSE, allocation = "exact")
    })
  )

})

