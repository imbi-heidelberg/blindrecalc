context("test the restricted ml function")

test_that("p_rml works correctly", {
  reml1 <- diff(p_rml(0.4, 0.6, 1, 0.2))
  reml2 <- diff(p_rml(0.1, 0.2, 3, 0.4))
  reml3 <- diff(p_rml(0.7, 0.7, 2, 0.3))

  expect_equal(reml1, -0.2)
  expect_equal(reml2, -0.4)
  expect_equal(reml3, -0.3)
})
