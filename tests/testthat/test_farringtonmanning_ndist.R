context("test n_dist for Farrington-Manning test")

test_that("errors are defined correctly", {
  d1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1)
  expect_error(n_dist(d1, 1.01, TRUE))

  d2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1, n_max = 301)
  expect_error(n_dist(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(n_dist(d2, 20, 0.5, TRUE, "exact"))

  expect_error(n_dist(d2, c(20, 30), c(0.6, 0.7), TRUE, "approximate"))

})


test_that("summary option", {
c
  n <- n_dist(d, n1 = 21, nuisance = 0.23, summary = FALSE)

  n_table <- n_dist(d, n1 = 21, nuisance = 0.23, summary = TRUE)

  expect_equal(as.vector(n_table), as.vector(summary(n$`p = 0.23`)))

})



test_that("plot option", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                              delta = 0, delta_NI = 0.1, n_max = 210)

  n <- n_dist(d, n1 = 21, nuisance = 0.8, allocation = "approximate",
              summary = FALSE, plot = TRUE)

  expect_equal(class(n), "list")

})
