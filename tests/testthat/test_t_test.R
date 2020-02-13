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
  set.seed(2020)

  design <- setupStudent(alpha = .025, beta = .2, r = 2, delta = 3.5,
                         delta_NI = -1.5, n_max = Inf)

  expect_gte(
    toer(design, 10, 5, TRUE, its),
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
    toer(des, 10, 5, TRUE, its),
    design@alpha
  )

})
