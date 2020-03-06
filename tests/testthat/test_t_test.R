context("test t-test")

its <- 1e5

test_that("Example of Lu can be reproduced", {
  design <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5,
                         delta_NI = 0, n_max = Inf)

  expect_equal(
    toer(design, 40, 5.5, TRUE, its),
    design@alpha,
    tolerance = 2 / sqrt(its), scale = 1)


  expect_equal(
    pow(design, 40, 5.5, TRUE, its),
    1 - design@beta,
    tolerance = 2 / sqrt(its), scale = 1)

})


test_that("Alpha can be adjusted in non-inferiority case", {
  design <- setupStudent(alpha = .025, beta = .2, r = 2, delta = 3.5,
                         delta_NI = 1.5, n_max = Inf)

  expect_gte(
    toer(design, 10, 5, TRUE, its, seed = 2020),
    design@alpha
  )

  alpha_adj <- adjusted_alpha(design, 10, 5, 1e-4, its)

  expect_lte(
    alpha_adj,
    design@alpha
  )

  des       <- design
  des@alpha <- alpha_adj
  expect_lte(
    toer(des, 10, 5, TRUE, its, seed = 1702),
    design@alpha
  )

})


test_that("Vectorization works", {
  design <- setupStudent(alpha = .025, beta = .2, r = 2, delta = 3.5,
                         delta_NI = 1.5, n_max = Inf)

  expect_error(
    n_dist(design, c(10,20), c(4, 5), T, T, 1e4)
  )

  expect_equal(
    as.numeric(unlist(n_dist(design, c(10, 20), 5, FALSE, FALSE, 1e4, 2020))),
    as.numeric(unlist(sapply(c(10, 20), function(n1) n_dist(design, n1, 5, FALSE, FALSE, 1e4, 2020))))
  )
})


test_that("alternative equals smaller", {
  #test errors
  expect_error(setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5,
                            delta_NI = 0, alternative = "smaller", n_max = Inf)
  )

  expect_error(setupStudent(alpha = .025, beta = .2, r = 1, delta = - 3.5,
                            delta_NI = 0, alternative = "greater", n_max = Inf)
  )

  expect_error(setupStudent(alpha = .025, beta = .2, r = 1, delta = - 3.5,
                            delta_NI = -1, alternative = "smaller", n_max = Inf)
  )


  # test toer
  design_smaller <- setupStudent(alpha = 0.025, beta = .1, r = 2, delta = -2,
                                 delta_NI = 0, alternative = "smaller", n_max = 100)
  design_greater <- setupStudent(alpha = 0.025, beta = .1, r = 2, delta = 2,
                                 delta_NI = 0, alternative = "greater", n_max = 100)
  expect_equal(toer(design_smaller, 20, 3, TRUE, 1e4, 42),
               toer(design_greater, 20, 3, TRUE, 1e4, 42))

})

