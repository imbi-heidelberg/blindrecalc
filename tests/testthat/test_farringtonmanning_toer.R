context("test toer for the Farrington-Manning test")

test_that("toer gives same values as pow for delta = 0", {
  design1 <- suppressWarnings(setupFarringtonManning(alpha = 0.025, beta = 0.2,
                                                     r = 1, delta = 0,
                                                     delta_NI = 0, n_max = 250))
  toer1 <- toer(design1, n1 = 80, nuisance = 0.4, recalculation = FALSE)
  power1 <- pow(design1, n1 = 80, nuisance = 0.4, recalculation = FALSE)

  expect_equal(toer1, power1)
})


test_that("errors are defined correctly", {
  d1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1)
  expect_error(toer(d1, n1 = 30, nuisance = 1.1, TRUE))

  d2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1, n_max = 301)
  expect_error(toer(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(toer(d2, 20, 0.5, TRUE, "exact"))

  expect_error(toer(d2, d2@n_max + 1, 0.7, TRUE, "approximate"))

  expect_error(toer(d2, c(20, 30), c(0.6, 0.7), TRUE, "approximate"))

})




test_that("vectorization in n1 works", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.25)
  n <- n_fix(d, 0.25)

  expect_equal(
    toer(d, n1 = c(n/2, n+10), nuisance = 0.25, recalculation = TRUE, allocation = "approximate"),
    sapply(c(n/2, n+10), function(n) {
      toer(d, n1 = n, nuisance = 0.25, recalculation = TRUE, allocation = "approximate")
    })
  )

  expect_equal(
    toer(d, n1 = c(n/2, n+10), nuisance = 0.25, recalculation = FALSE, allocation = "approximate"),
    sapply(c(n/2, n+10), function(n) {
      toer(d, n1 = n, nuisance = 0.25, recalculation = FALSE, allocation = "approximate")
    })
  )

})


test_that("vectorization in nuisance works", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.2)
  expect_equal(
    toer(d, n1 = 20, nuisance = c(0.4, 0.6), recalculation = TRUE, allocation = "approximate"),
    sapply(c(0.4, 0.6), function(p) {
      toer(d, n1 = 20, nuisance = p, recalculation = TRUE, allocation = "approximate")
    })
  )

  expect_equal(
    toer(d, n1 = 20, nuisance = c(0.4, 0.6), recalculation = FALSE, allocation = "exact"),
    sapply(c(0.4, 0.6), function(p) {
      toer(d, n1 = 20, nuisance = p, recalculation = FALSE, allocation = "exact")
    })
  )

})


