context("test pow for the Farrington-Manning test")

test_that("pow works for the fixed design (Friede et al. 2007)", {
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


  expect_equal(pow1, 0.8005762, tolerance = 1e-7)
  expect_equal(pow2, 0.7948942, tolerance = 1e-7)
  expect_equal(pow3, 0.8005885, tolerance = 1e-7)
  expect_equal(pow4, 0.8000787, tolerance = 1e-7)
  expect_equal(pow5, 0.8009629, tolerance = 1e-7)
  expect_equal(pow6, 0.8026954, tolerance = 1e-7)
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

