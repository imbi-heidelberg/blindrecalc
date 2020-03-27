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



test_that("fixed sample size", {
  des <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 2,
                      delta_NI = 0, n_max = Inf)
  n <- n_fix(des, 4)

  type_one <- simulation(design = des, n1 = n, nuisance = 4, recalculation = FALSE,
                         delta_true = 0, iters = its, seed = 2020,
                         allocation = "approximate")

  expect_equal(type_one$rejection_probability,
               des@alpha,
               tolerance = 2 / sqrt(its), scale = 1)

  expect_equal(type_one$sample_sizes,
               rep(n, its))
})



test_that("exact allocation", {
  des <- setupStudent(alpha = 0.05, beta = 0.1, r = 1, delta = 3, delta_NI = 1,
                      alternative = "greater", n_max = 299)

  expect_error(simulation(des, 21, 5, TRUE, 2, 1e4, 2020, "exact"))
  expect_error(simulation(des, 20, 5, TRUE, 2, 1e4, 2020, "exact"))

  des@n_max <- 300

  expect_lte(
    simulation(des, 20, 5, TRUE, 1, 1, 2020, "approximate")$sample_sizes,
    simulation(des, 20, 5, TRUE, 1, 1, 2020, "exact")$sample_sizes
  )

})


test_that("summary function", {
  des <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 2,
                      delta_NI = 0, n_max = Inf)

  n <- data.frame(n_dist(design, 20, 3.5, FALSE, FALSE, its, 2020))

  n_table <- n_dist(design, 20, 3.5, TRUE, FALSE, its, 2020)

  expect_equal(as.vector(n_table), as.vector(summary(n)))

})
