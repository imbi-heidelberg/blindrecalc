#' Simulate Rejection Probability and Sample Size for Student's t-Test
#'
#' This function simulates the probability that a test defined by
#' \code{\link{setupStudent}} rejects the null hypothesis.
#' Note that here the nuisance parameter \code{nuisance} is the variance
#' of the outcome variable sigma^2.
#'
#' @template methods
#' @template recalculation
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
simulation <- function(design, n1, nuisance, recalculation = TRUE, Delta_star, iters = 1000, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  alloc <- design@r / (1 + design@r)^2

  # Step 1
  z1 <- stats::rnorm(n = iters, mean = 0, sd = 1)
  v1 <- stats::rchisq(n = iters, df = n1 - 2)

  # Step 2
  var_hat <- nuisance^2 / (n1 - 1) * (v1 + (z1 + sqrt(n1 * alloc) * Delta_star / nuisance)^2)

  # Step 3
  if (recalculation == FALSE) {
    n <- n1
    } else {
      n_recalc <- 1 / alloc * (stats::qnorm(1 - design@alpha) + stats::qnorm(1 - design@beta))^2 /
        (design@delta - design@delta_NI)^2 * var_hat
      n        <- sapply(n_recalc, function(m) min(design@n_max, max(ceiling(m), n1)))
    }

  n2 <- n - n1

  f <- function(i) {
    # Step 4
    if(n2[i] == 0) {
      test_statistic <- (z1[i] + sqrt(n1 * alloc) * (Delta_star - design@delta_NI) / nuisance) / sqrt(v1[i] / (n1 - 2))
      } else {
        # Step 5
        z2 <- stats::rnorm(n = 1, mean = 0, sd = 1)
        w2 <- stats::rchisq(n = 1, df = n2[i] - 1)
        v2 <- w2 + (sqrt(n2[i] / n[i]) * z1[i] - sqrt(n1 / n[i]) * z2)^2

        test_statistic <-
          (sqrt(n1 / n[i]) * z1[i] + sqrt(n2[i] / n[i]) * z2 + sqrt(n[i] * alloc) *
             (Delta_star - design@delta_NI) / nuisance) / sqrt((v1[i] + v2) / (n[i] - 2))
        }
    return(test_statistic)
    }

  test_statistic <- sapply(seq(1, iters, 1), f)
  critical_value <- stats::qt(1 - design@alpha, df = n - 2)
  reject         <- ifelse(test_statistic >= critical_value, 1, 0)

  return(list("rejection_probability" = mean(reject),
              "sample_sizes" = n
  ))
}




#' @template iters
#' @rdname toer
#' @export
setMethod("toer", signature("Student"),
          function(design, n1, nuisance, recalculation = TRUE, iters = 1e4, seed = NULL, ...) {
            if (length(nuisance) > 1 && length(n1) > 1) {
              stop("Either the nuisance parameter or the internal pilot study sample size must be of length 1!")
            }

            if (length(n1) == 1) {
              return(sapply(nuisance, function(sigma)
                simulation(design, n1, sigma, recalculation, design@delta_NI, iters, seed, ...)$rejection_probability))
            }

            if (length(nuisance) == 1) {
              return(sapply(n1, function(n1)
                simulation(design, n1, nuisance, recalculation, design@delta_NI, iters, seed, ...)$rejection_probability))
            }
          })




#' @template iters
#' @rdname pow
#' @export
setMethod("pow", signature("Student"),
          function(design, n1, nuisance, recalculation = TRUE, iters = 1e4, seed = NULL, ...) {
            if (length(nuisance) > 1 && length(n1) > 1) {
              stop("Either the nuisance parameter or the internal pilot study sample size must be of length 1!")
            }

            if (length(n1) == 1) {
              return(sapply(nuisance, function(sigma)
                simulation(design, n1, sigma, recalculation, design@delta, iters, seed, ...)$rejection_probability))
            }

            if (length(nuisance) == 1) {
              return(sapply(n1, function(n1)
                simulation(design, n1, nuisance, recalculation, design@delta, iters, seed, ...)$rejection_probability))
            }

          })




#' @template iters
#' @rdname sample_size_dist
#' @export
setMethod("sample_size_dist", signature("Student"),
          function(design, n1, nuisance, summary = TRUE, plot = FALSE, iters = 1e4, seed = NULL, ...) {
            if (length(nuisance) > 1 && length(n1) > 1) {
              stop("Either the nuisance parameter or the internal pilot study sample size must be of length 1!")
            }

            if (length(n1) == 1) {
              n <- sapply(nuisance, function(sigma)
                simulation(design, n1, sigma, recalculation = TRUE, design@delta, iters, seed, ...)$sample_sizes)

              if (plot == TRUE) {
                graphics::par(c(list(mfrow = c(1, length(nuisance)))))
                for (i in 1:length(nuisance)) {
                  graphics::boxplot(n[, i], range = 0, xlab = paste(expression(sigma),"=",nuisance[i]),
                                    ylim = c(min(n), max(n)))
                }
              }

              n <- data.frame(n)
              for (i in 1:ncol(n))
                colnames(n)[i] <- paste(expression(sigma),"=",nuisance[i])
            }

            if (length(nuisance) == 1) {
              n <- sapply(n1, function(n1)
                simulation(design, n1, nuisance, recalculation = TRUE, design@delta, iters, seed, ...)$sample_sizes)

              if (plot == TRUE) {
                graphics::par(c(list(mfrow = c(1, length(n1)))))
                for (i in 1:length(n1)) {
                  graphics::boxplot(n[, i], range = 0, xlab = paste(expression(n_1),"=",n1[i]),
                                    ylim = c(min(n), max(n)))
                }
              }

              n <- data.frame(n)
              for (i in 1:ncol(n))
                colnames(n)[i] <- paste(expression(n_1),"=",n1[i])
            }

            if (summary == TRUE) return(summary(n))
            else return(n)

          })




#' @rdname n_fix
#' @export
setMethod("n_fix", signature("Student"),
          function(design, nuisance, ...) {
            sapply(nuisance, function(sigma)
              (1 + design@r)^2 / design@r * (stats::qnorm(1 - design@alpha) + stats::qnorm(1 - design@beta))^2 /
                (design@delta - design@delta_NI)^2 * sigma^2
            )
          })




#' @param tol desired absolute tolerance
#' @template iters
#' @rdname adjusted_alpha
#'
#' @details In the case of the Student's t-test, the adjusted alpha is calculated
#' using the algorithm by Kieser and Friede (2000):
#' "Re-calculating the sample size in internal pilot study designs
#' with control of the type I error rate"
#'
#' @export
setMethod("adjusted_alpha", signature("Student"),
          function(design, n1, nuisance, tol, iters = 1e4, seed = NULL, ...) {
            alpha_max <- function(alp) {
              d       <- design
              d@alpha <- alp
              return(max(toer(d, n1, nuisance, TRUE, iters, seed)))
            }

            alpha_adj <- design@alpha
            alpha_act <- alpha_max(alpha_adj)

            while(alpha_act - design@alpha > tol) {
              alpha_adj <- alpha_adj * design@alpha / alpha_act
              alpha_act <- alpha_max(alpha_adj)
            }

            return(alpha_adj)

          })
