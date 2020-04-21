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
  expect_error(toer(d1, 1.1, TRUE))

  d2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1, n_max = 301)
  expect_error(toer(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(toer(d2, 20, 0.5, TRUE, "exact"))

  expect_error(toer(d2, d2@n_max + 1, 0.7, TRUE, "approximate"))

  expect_error(toer(d2, c(20, 30), c(0.6, 0.7), TRUE, "approximate"))

})
