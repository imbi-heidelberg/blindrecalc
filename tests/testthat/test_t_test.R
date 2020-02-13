context("test t-test")

design <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5,
                       delta_NI =  0, n_max = Inf)

test_that("Example of Lu can be reproduced", {
  its <- 1e5

  expect_equal(
    toer(design, 40, 5.5, TRUE, its),
    design@alpha,
    tolerance = 2 / sqrt(its), scale = 1)


  expect_equal(
    pow(design, 40, 5.5, TRUE, its),
    1 - design@beta,
    tolerance = 2 / sqrt(its), scale = 1)

})
