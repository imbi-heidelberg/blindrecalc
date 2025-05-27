context("test t-test sample size distribution")

its <- 1e5

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


test_that("ndist respects n_max", {
  des <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 2,
                      delta_NI = 0, n_max = 300)
  distri <- n_dist(des, 20, 3.5, TRUE, FALSE, its, 2020)
  expect_equal(distri[nrow(distri), ], "Max.   :300  ")

})
