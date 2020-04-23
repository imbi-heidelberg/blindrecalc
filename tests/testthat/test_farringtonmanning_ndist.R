context("test n_dist for Farrington-Manning test")

test_that("errors are defined correctly", {
  d1 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1)
  expect_error(n_dist(d1, n1 = 30, nuisance = -0.01, TRUE))

  d2 <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                               delta = 0, delta_NI = 0.1, n_max = 301)
  expect_error(n_dist(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(n_dist(d2, 20, 0.5, TRUE, "exact"))

  expect_error(n_dist(d2, c(21, 30), c(0.6, 0.7), TRUE, "approximate"))

})


test_that("n_dist works when n1 and nuisance are of length 1", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                              delta = 0, delta_NI = 0.1, n_max = 210)

  n <- n_dist(d, n1 = 21, nuisance = 0.23, summary = FALSE)
  n_table <- n_dist(d, n1 = 21, nuisance = 0.23, summary = TRUE)
  expect_equal(as.vector(n_table), as.vector(summary(n$`p = 0.23`)))

  n_plot <- n_dist(d, n1 = 21, nuisance = 0.8, allocation = "exact",
                   summary = FALSE, plot = TRUE)
  expect_equal(class(n_plot), "list")

})



test_that("n_dist works for multiple n1 values", {
  d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 2,
                              delta = 0, delta_NI = 0.1, n_max = 210)

  n <- n_dist(d, n1 = c(21, 30), nuisance = 0.23, summary = FALSE)
  n_table <- n_dist(d, n1 = c(21, 30), nuisance = 0.23, summary = TRUE)
  expect_equal(as.vector(n_table), as.vector(c(summary(n$`n1 = 21`), summary(n$`n1 = 30`))))

  n_plot <- n_dist(d, n1 = c(21, 30), nuisance = 0.8, allocation = "approximate",
                   summary = FALSE, plot = TRUE)
  expect_equal(class(n_plot), "list")

})

