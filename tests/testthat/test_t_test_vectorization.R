context("t-test vectorization")

its <- 1e4

design <- setupStudent(alpha = .025, beta = .2, r = 2, delta = 3.5,
                       delta_NI = 1.5, n_max = Inf)


test_that("error messages when n1 and nuisance are vectors", {
  expect_error(
    n_dist(design, c(10, 20), c(4, 5), T, T, its)
  )

  expect_error(
    toer(design, c(20, 30), c(2, 3, 5), F, F, its)
  )

  expect_error(
    pow(design, c(20, 30), c(2, 3, 5), F, F, its)
  )
})


test_that("Vectorization in n1 works", {
  expect_equal(
    as.numeric(unlist(n_dist(design, c(10, 20), 5, FALSE, FALSE, its, 2020))),
    as.numeric(unlist(sapply(c(10, 20), function(n1) n_dist(design, n1, 5, FALSE, FALSE, its, 2020))))
  )

  expect_equal(
    toer(design, c(10, 20), 5, TRUE, its, 2020),
    sapply(c(10, 20), function(n1) toer(design, n1, 5, TRUE, its, 2020))
  )

  expect_equal(
    pow(design, c(10, 20), 5, TRUE, its, 2703),
    sapply(c(10, 20), function(n1) pow(design, n1, 5, TRUE, its, 2703))
  )

})



test_that("Vectorization in nuisance works", {
  expect_equal(
    as.numeric(unlist(n_dist(design, 17, c(2, 7), FALSE, FALSE, its, 2020))),
    as.numeric(unlist(sapply(c(2, 7), function(sigma) n_dist(design, 17, sigma, FALSE, FALSE, its, 2020))))
  )

  expect_equal(
    toer(design, 23, c(1, 8.5), TRUE, its, 2020),
    sapply(c(1, 8.5), function(sigma) toer(design, 23, sigma, TRUE, its, 2020))
  )

  expect_equal(
    pow(design, 3, c(2, 4, 6), TRUE, its, 2703),
    sapply(c(2, 4, 6), function(sigma) pow(design, 3, sigma, TRUE, its, 2703))
  )

})

