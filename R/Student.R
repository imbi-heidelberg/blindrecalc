#' Simulate Rejection Probability and Sample Size for Student's t-Test
#'
#' This function simulates the probability that a test defined by
#' \code{\link{setupStudent}} rejects the null hypothesis.
#' Note that here the nuisance parameter \code{nuisance} is the variance
#' of the outcome variable sigma^2.
#'
#' @template methods
#' @param Delta_star effect measure under which the rejection probabilities are computed
#' @template iters
#' @template dotdotdot
#'
#' @details The implementation follows the algorithm in Lu (2019):
#' Distribution of the two-sample t-test statistic following blinded
#' sample size re-estimation.
#' Pharmaceutical Statistics 15: 208-215.
#'
#' @export
simulation <- function(design, n1, nuisance, recalculation = TRUE, Delta_star, iters = 1000, ...) {
            r <- design@r

            # Step 1
            z1 <- stats::rnorm(n = iters, mean = 0, sd = 1)
            v1 <- stats::rchisq(n = iters, df = n1 - 2)

            # Step 2
            var_hat <- nuisance^2 / (n1 - 1) * (v1 + (z1 + sqrt(r * n1 / (1 + r)^2) * Delta_star / nuisance)^2)

            # Step 3
            if (recalculation == FALSE) {
              n <- n1
              } else {
                n_recalc <- (1 + r)^2 / r * (stats::qnorm(1 - design@alpha) + stats::qnorm(1 - design@beta))^2 /
                  (design@delta - design@delta_NI)^2 * var_hat
                n        <- sapply(n_recalc, function(m) min(design@n_max, max(ceiling(m), n1)))
              }
            n2 <- n - n1

            f <- function(i) {
              # Step 4
              if(n2[i] == 0) {
                test_statistic <- (z1[i] + sqrt(n1 * r / (1 + r)^2) * (Delta_star - design@delta_NI) / nuisance) / sqrt(v1[i] / (n1 - 2))
                } else {
                  # Step 5
                  z2 <- stats::rnorm(n = 1, mean = 0, sd = 1)
                  w2 <- stats::rchisq(n = 1, df = n2[i] - 1)
                  v2 <- w2 + (sqrt(n2[i] / n[i]) * z1[i] - sqrt(n1 / n[i]) * z2)^2

                  test_statistic <-
                    (sqrt(n1 / n[i]) * z1[i] + sqrt(n2[i] / n[i]) * z2 + sqrt(n[i] * r / (1 + r)^2) *
                       (Delta_star - design@delta_NI) / nuisance) / sqrt((v1[i] + v2) / (n[i] - 2))
                  }
              return(test_statistic)
            }

            test_statistic <- sapply(seq(1, iters, 1), f)

            critical_value <- stats::qt(1 - design@alpha, df = n - 2)
            reject <- ifelse(test_statistic >= critical_value, 1, 0)

            return(list(
              "rejection_probability" = mean(reject),
              "sample_sizes" = n
              ))
}


#' @template iters
#' @rdname toer
#' @export
setMethod("toer", signature("Student"),
          function(design, n1, nuisance, recalculation = TRUE, iters, ...) {
            sapply(nuisance, function(sigma)
              simulation(design, n1, sigma, recalculation, design@delta_NI, iters, ...)$rejection_probability)
          })


#' @template iters
#' @rdname pow
#' @export
setMethod("pow", signature("Student"),
          function(design, n1, nuisance, recalculation = TRUE, iters, ...) {
            sapply(nuisance, function(sigma)
              simulation(design, n1, sigma, recalculation, design@delta, iters, ...)$rejection_probability)
          })


#' @rdname n_fix
#' @export
setMethod("n_fix", signature("Student"),
          function(design, nuisance, ...) {
            (1 + design@r) * 2 * (stats::qnorm(1 - design@alpha) + stats::qnorm(1 - design@beta))^2 /
              (design@delta - design@delta_NI)^2 * nuisance^2
          })

