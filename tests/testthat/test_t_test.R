context("test t-test")

its <- 1e5

test_that("Examples of Lu (2019) can be reproduced", {
  set.seed(123)
  design1 <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5,
                         delta_NI = 0, n_max = Inf)

  expect_equal(
    toer(design1, 40, 5.5, TRUE, its),
    design1@alpha,
    tolerance = 2 / sqrt(its), scale = 1)

  expect_equal(
    pow(design1, 40, 5.5, TRUE, its),
    1 - design1@beta,
    tolerance = 2 / sqrt(its), scale = 1)


  design2 <- setupStudent(alpha = .05, beta = 1-.903, r = 1, delta = 0,
                          delta_NI = 0.5, n_max = Inf)

  expect_equal(
    toer(design2, 60, 1, TRUE, its),
    design2@alpha,
    tolerance = 2 / sqrt(its), scale = 1)

  expect_equal(
    pow(design2, 60, 1, TRUE, its),
    1 - design2@beta,
    tolerance = 2 / sqrt(its), scale = 1)

})




test_that("Different values of alpha and beta work", {
  # try design with low power and small alpha
  design1 <- setupStudent(alpha = .025, beta = .8, r = 1, delta = 1,
                         delta_NI = 0, n_max = Inf)
  expect_equal(
    toer(design1, 10, 2, TRUE, its),
    design1@alpha,
    tolerance = 2 / sqrt(its), scale = 1)
  expect_equal(
    pow(design1, 10, 2, TRUE, its),
    1 - design1@beta,
    tolerance = 2 / sqrt(its), scale = 1)


  # try design with low power and large alpha
  design2 <- setupStudent(alpha = .2, beta = .5, r = 1, delta = 1,
                         delta_NI = 0, n_max = Inf)
  expect_equal(
    toer(design2, 20, 5, TRUE, its),
    design2@alpha,
    tolerance = 2 / sqrt(its), scale = 1)
  expect_equal(
    pow(design2, 20, 5, TRUE, its),
    1 - design2@beta,
    tolerance = 2 / sqrt(its), scale = 1)


  # try design with large power and large alpha
  design3 <- setupStudent(alpha = .3, beta = .05, r = 1, delta = 1,
                         delta_NI = 0, n_max = Inf)
  expect_equal(
    toer(design3, 50, 3, TRUE, its),
    design3@alpha,
    tolerance = 2 / sqrt(its), scale = 1)
  expect_equal(
    pow(design3, 50, 3, TRUE, its),
    1 - design3@beta,
    tolerance = 2 / sqrt(its), scale = 1)

})


test_that("Smaller power requires smaller sample sizes", {
  design1 <- setupStudent(alpha = .05, beta = .1, r = 1, delta = 1,
                          delta_NI = 0, n_max = Inf)
  design2 <- setupStudent(alpha = .05, beta = .5, r = 1, delta = 1,
                          delta_NI = 0, n_max = Inf)
  design3 <- setupStudent(alpha = .05, beta = .8, r = 1, delta = 1,
                          delta_NI = 0, n_max = Inf)
  design4 <- setupStudent(alpha = .05, beta = .99, r = 1, delta = 1,
                          delta_NI = 0, n_max = Inf)

  expect_lte(
    mean(unlist(n_dist(design4, 10, 5, FALSE, FALSE))),
    mean(unlist(n_dist(design3, 10, 5, FALSE, FALSE)))
  )

  expect_lte(
    mean(unlist(n_dist(design3, 10, 5, FALSE, FALSE))),
    mean(unlist(n_dist(design2, 10, 5, FALSE, FALSE)))
  )

  expect_lte(
    mean(unlist(n_dist(design2, 10, 5, FALSE, FALSE))),
    mean(unlist(n_dist(design1, 10, 5, FALSE, FALSE)))
  )

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

  expect_error(setupStudent(alpha = .025, beta = .2, r = 1, delta = - 3.5,
                            delta_NI = 1, alternative = "smaller", n_max = Inf)
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


test_that("summary option", {
  des <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 2,
                      delta_NI = 0, n_max = Inf)

  n <- data.frame(n_dist(des, 20, 3.5, FALSE, FALSE, its, 2020))

  n_table <- n_dist(des, 20, 3.5, TRUE, FALSE, its, 2020)

  expect_equal(as.vector(n_table), as.vector(summary(n)))

})



test_that("plot option", {
  des <- setupStudent(alpha = .03, beta = .15, r = 1.5, delta = 2.5,
                      delta_NI = 0, n_max = 350)

  n <- n_dist(des, 20, 1, FALSE, TRUE, its, 42)

  expect_equal(class(n), "data.frame")

})
