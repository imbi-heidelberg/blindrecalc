context("test t-test")

design <- setupStudent(.025, .2, 1, 3.5, 0, Inf)

test_that("Example of Lu can be reproduced", {
  its <- 1e5

  expect_equal(
    toer(design, 20, 5.5, TRUE, its),
    design@alpha,
    tolerance = 2 / sqrt(its), scale = 1)


  expect_equal(
    pow(design, 20, 5.5, TRUE, its),
    1 - design@beta,
    tolerance = 2 / sqrt(its), scale = 1)

})
