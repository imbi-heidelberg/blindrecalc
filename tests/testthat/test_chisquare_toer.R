context("test toer for the ChiSquare-test")

test_that("toer gives same values as pow for delta = 0", {
  design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0)
  size1 <- toer(design1, n1 = 100, nuisance = 0.2, recalculation = FALSE)
  pow1 <- pow(design1, n1 = 100, nuisance = 0.2, recalculation = FALSE)
  expect_equal(size1, pow1)

  skip_on_cran()
  design2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0, n_max = 150)
  size2 <- toer(design2, n1 = 75, nuisance = 0.4, recalculation = TRUE,
                allocation = "approximate")
  pow2 <- pow(design2, n1 = 75, nuisance = 0.4, recalculation = TRUE,
               allocation = "approximate")
  expect_equal(size2, pow2)
})


test_that("errors are thrown correctly", {
  d1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.1)
  expect_error(toer(d1, n1 = 20, nuisance = 1.1, TRUE))

  d2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 2, delta = 0.1, n_max = 301)
  expect_error(toer(d2, 21, 0.5, TRUE, "exact"))
  d2@n_max <- 300
  expect_error(toer(d2, 20, 0.5, TRUE, "exact"))

  expect_error(toer(d2, d2@n_max + 1, 0.7, TRUE, "approximate"))

  expect_error(toer(d2, c(20, 30), c(0.6, 0.7), TRUE, "approximate"))

})


test_that("vectorization in n1 works", {
  d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
  expect_equal(
    toer(d, n1 = c(10, 20), nuisance = 0.25, recalculation = TRUE, allocation = "approximate"),
    sapply(c(10, 20), function(n) {
      toer(d, n1 = n, nuisance = 0.25, recalculation = TRUE, allocation = "approximate")
    })
  )

  expect_equal(
    toer(d, n1 = c(10, 20), nuisance = 0.25, recalculation = FALSE, allocation = "approximate"),
    sapply(c(10, 20), function(n) {
      toer(d, n1 = n, nuisance = 0.25, recalculation = FALSE, allocation = "approximate")
    })
  )

})



test_that("alternative can be 'smaller'", {
  d1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2, alternative = "greater")
  d2 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2, alternative = "smaller")

  expect_equal(
    toer(d1, n1 = 20, nuisance = 0.5, recalculation = TRUE),
    toer(d2, n1 = 20, nuisance = 0.5, recalculation = TRUE)
  )

  expect_equal(
    toer(d1, n1 = 20, nuisance = 0.5, recalculation = FALSE),
    toer(d2, n1 = 20, nuisance = 0.5, recalculation = FALSE)
  )
})
