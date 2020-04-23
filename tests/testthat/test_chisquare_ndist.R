context("test n_dist for ChiSquare test")

test_that("error messages are thrown correctly", {
  d1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.1)
  expect_error(n_dist(d1, n1 = 20, nuisance = 1.1, TRUE))

  d2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 2, delta = 0.1, n_max = 301)
  expect_error(n_dist(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(n_dist(d2, 20, 0.5, TRUE, "exact"))

  expect_error(n_dist(d2, c(21, 30), c(0.6, 0.7), TRUE, "approximate"))

})


test_that("n_dist works for multiple n1 values", {
  d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)

  n <- n_dist(d, n1 = c(20, 40), nuisance = 0.4, summary = FALSE)
  n_table <- n_dist(d, n1 = c(20, 40), nuisance = 0.4, summary = TRUE)
  expect_equal(as.vector(n_table), as.vector(c(summary(n$`n1 = 20`), summary(n$`n1 = 40`))))

  n_plot <- n_dist(d, n1 = c(20, 40), nuisance = 0.4, allocation = "approximate",
                   summary = FALSE, plot = TRUE)
  expect_equal(class(n_plot), "list")

})


test_that("n_dist works for multiple nuisance values", {
  d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)

  n <- n_dist(d, n1 = 30, nuisance = c(0.35, 0.4), summary = FALSE)
  n_table <- n_dist(d, n1 = 30, nuisance = c(0.35, 0.4), summary = TRUE)
  expect_equal(as.vector(n_table), as.vector(c(summary(n$`p = 0.35`), summary(n$`p = 0.4`))))

  n_plot <- n_dist(d, n1 = 30, nuisance = c(0.35, 0.4), allocation = "approximate",
                   summary = FALSE, plot = TRUE)
  expect_equal(class(n_plot), "list")
})
