#' Simulate Rejection Probability and Sample Size for Student's t-Test
#'
#' This function simulates the probability that a test defined by
#' \code{\link{setupStudent}} rejects the null hypothesis.
#' Note that here the nuisance parameter \code{nuisance} is the variance
#' of the outcome variable sigma^2.
#'
#' @template methods_student
#' @template recalculation
#' @param delta_true effect measure under which the rejection probabilities are computed
#' @template iters
#' @template allocation
#' @template dotdotdot
#'
#' @return Simulated rejection probabilities and sample sizes for
#'    each nuisance parameter.
#'
#' @details The implementation follows the algorithm in Lu (2019):
#' Distribution of the two-sample t-test statistic following blinded
#' sample size re-estimation.
#' Pharmaceutical Statistics 15: 208-215.
#' Since Lu (2019) assumes negative non-inferiority margins, the non-inferiority
#' margin of \code{design} is multiplied with -1 internally.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' simulation(d, n1 = 20, nuisance = 5.5, recalculation = TRUE, delta_true = 3.5)
#'
#' @export
simulation <- function(design, n1, nuisance, recalculation = TRUE, delta_true,
                       iters = 1000, seed = NULL, allocation = c("approximate", "exact"), ...) {

  if (!is.null(seed)) set.seed(seed)

  if (design@alternative == "smaller") {
    design@delta <- -design@delta
    delta_true   <- -delta_true
    design@r     <- 1 / design@r
  }

  # check if allocation can be done exactly
  allocation <- match.arg(allocation)
  if (allocation == "exact") {
    if (sum(n1 %% (design@r + 1) != 0) > 0) {
      stop("n1 cannot be allocated exactly!")
    }
    if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
      stop("n_max cannot be allocated exactly!")
    }
  }
  alloc <- design@r / (1 + design@r)^2

  # the following implements the 5 steps of the algorithm by Lu (2019), p.210
  ## Step 1
  z1 <- stats::rnorm(n = iters, mean = 0, sd = 1)
  v1 <- stats::rchisq(n = iters, df = n1 - 2)

  ## Step 2
  var_hat <- nuisance^2 / (n1 - 1) * (v1 + (z1 + sqrt(n1 * alloc) * delta_true / nuisance)^2)

  ## Step 3
  if (recalculation == FALSE) {
    n <- rep(n1, iters)
  } else {
    n_recalc <- ceiling(1 / alloc * (stats::qnorm(1 - design@alpha) + stats::qnorm(1 - design@beta))^2 /
      (design@delta - design@delta_NI)^2 * var_hat)
    if (allocation == "exact") {
      while ((sum(n_recalc %% (design@r + 1) != 0) > 0)) n_recalc <- n_recalc + 1
    }
    n <- sapply(n_recalc, function(m) min(design@n_max, max(m, n1)))
  }

  n2 <- n - n1

  f <- function(i) {
    ## Step 4
    if (n2[i] == 0) {
      test_statistic <- (z1[i] + sqrt(n1 * alloc) * (delta_true - design@delta_NI) / nuisance) / sqrt(v1[i] / (n1 - 2))
      } else {
        ## Step 5
        z2 <- stats::rnorm(n = 1, mean = 0, sd = 1)
        w2 <- stats::rchisq(n = 1, df = n2[i] - 1)
        v2 <- w2 + (sqrt(n2[i] / n[i]) * z1[i] - sqrt(n1 / n[i]) * z2)^2

        test_statistic <-
          (sqrt(n1 / n[i]) * z1[i] + sqrt(n2[i] / n[i]) * z2 + sqrt(n[i] * alloc) *
             (delta_true - design@delta_NI) / nuisance) / sqrt((v1[i] + v2) / (n[i] - 2))
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


#' Type I Error Rate
#'
#' Computes the type I error rate of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods_student
#' @template recalculation
#' @template iters
#' @template allocation
#' @template dotdotdot
#'
#' @return One type I error rate value for every nuisance parameter
#'  and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' toer(d, n1 = 20, nuisance = 5.5, recalculation = TRUE)
#'
#' @rdname toer.Student
#' @export
setMethod("toer", signature("Student"),
          function(design, n1, nuisance, recalculation = TRUE, iters = 1e4, seed = NULL,
                   allocation = c("approximate", "exact"), ...) {

            if (length(nuisance) > 1 && length(n1) > 1) {
              stop("Either the nuisance parameter or the internal pilot study sample size must be of length 1!")
            }

            # apply simulation function at the non-inferiority boundary (i.e., the null hypothesis)
            if (length(n1) == 1) {
              return(sapply(nuisance, function(sigma)
                simulation(design, n1, sigma, recalculation, design@delta_NI, iters, seed, allocation, ...)$rejection_probability))
            } else if (length(nuisance) == 1) {
              return(sapply(n1, function(n1)
                simulation(design, n1, nuisance, recalculation, design@delta_NI, iters, seed, allocation, ...)$rejection_probability))
            }
          })



#' Power
#'
#' Calculates the power of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods_student
#' @template recalculation
#' @template iters
#' @template allocation
#' @template dotdotdot
#'
#' @return One power value for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' pow(d, n1 = 20, nuisance = 5.5, recalculation = TRUE)
#'
#' @rdname pow.Student
#' @export
setMethod("pow", signature("Student"),
          function(design, n1, nuisance, recalculation = TRUE, iters = 1e4, seed = NULL,
                   allocation = c("approximate", "exact"), ...) {
            if (length(nuisance) > 1 && length(n1) > 1) {
              stop("Either the nuisance parameter or the internal pilot study sample size must be of length 1!")
            }

            # apply simulation function at the specified effect size (i.e., the alternative hypothesis)
            if (length(n1) == 1) {
              return(sapply(nuisance, function(sigma)
                simulation(design, n1, sigma, recalculation, design@delta, iters, seed, allocation, ...)$rejection_probability))
            } else if (length(nuisance) == 1) {
              return(sapply(n1, function(n1)
                simulation(design, n1, nuisance, recalculation, design@delta, iters, seed, allocation, ...)$rejection_probability))
            }

          })




#' Distribution of the Sample Size
#'
#' Calculates the distribution of the total sample sizes of designs
#' with blinded sample size recalculation for different values of the
#' nuisance parameter or of n1.
#'
#' @template methods_student
#' @template summary
#' @template plot
#' @template iters
#' @template allocation
#' @param range this determines how far the plot whiskers extend out from the box.
#'    If range is positive, the whiskers extend to the most extreme data point
#'    which is no more than range times the interquartile range from the box.
#'    A value of zero causes the whiskers to extend to the data extremes.
#' @template dotdotdot
#'
#' @return Summary and/or plot of the sample size distribution for
#'   every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' n_dist(d, n1 = 20, nuisance = 5.5, summary = TRUE, plot = FALSE, seed = 2020)
#'
#' @rdname n_dist.Student
#' @export
setMethod("n_dist", signature("Student"),
          function(design, n1, nuisance, summary = TRUE, plot = FALSE, iters = 1e4,
                   seed = NULL, range = 0, allocation = c("approximate", "exact"), ...) {
            if (length(nuisance) > 1 && length(n1) > 1) {
              stop("Only one of n1 and nuisance can have length > 1.")
            }

            # create data frame that includes the simulated sample sizes
            if (length(n1) == 1) {
              n <- sapply(nuisance, function(sigma)
                simulation(design, n1, sigma, recalculation = TRUE, design@delta, iters, seed, allocation, ...)$sample_sizes)
              n <- data.frame(n)
              for (i in 1:ncol(n))
                colnames(n)[i] <- paste(expression(sigma),"=",nuisance[i])
            } else if (length(nuisance) == 1) {
              n <- sapply(n1, function(n1)
                simulation(design, n1, nuisance, recalculation = TRUE, design@delta, iters, seed, allocation, ...)$sample_sizes)
              n <- data.frame(n)
              for (i in 1:ncol(n))
                colnames(n)[i] <- paste(expression(n_1),"=",n1[i])
            }

            if (plot == TRUE) graphics::boxplot(n, range = range, ...)

            if (summary == TRUE) return(summary(n))
            else return(n)

          })



#' Adjusted level of significance
#'
#' This method returns an adjusted significance level that can be used
#' such that the actual type I error rate is preserved.
#'
#' @template methods_student
#' @param tol desired absolute tolerance
#' @template iters
#' @template dotdotdot
#'
#' @return Value of the adjusted significance level for every nuisance
#'  parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @details In the case of the Student's t-test, the adjusted alpha is calculated
#' using the algorithm by Kieser and Friede (2000):
#' "Re-calculating the sample size in internal pilot study designs
#' with control of the type I error rate". Statistics in Medicine 19: 901-911.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 0, delta_NI = 1.5, n_max = 848)
#' sigma <- c(2, 5.5, 9)
#' adjusted_alpha(design = d, n1 = 20, nuisance = sigma, tol = 1e-4, iters = 1e3)
#'
#' @rdname adjusted_alpha.Student
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

            # iteratively reduce the significance level until it is sufficiently small
            while(alpha_act - design@alpha > tol) {
              alpha_adj <- alpha_adj * design@alpha / alpha_act
              alpha_act <- alpha_max(alpha_adj)
            }

            return(alpha_adj)
        })



#' Fixed Sample Size
#'
#' Returns the sample size of a fixed design without sample size recalculation.
#'
#' @param design test statistic object
#' @param nuisance nuisance parameter
#' @template dotdotdot
#'
#' @return One value of the fixed sample size for every nuisance parameter
#'  and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' n_fix(design = d, nuisance = 5.5)
#'
#' @rdname n_fix.Student
#' @export
setMethod("n_fix", signature("Student"),
          function(design, nuisance, ...) {
            # apply known formula for fixed sample size
            sapply(nuisance, function(sigma)
              (1 + design@r)^2 / design@r * (stats::qnorm(1 - design@alpha) + stats::qnorm(1 - design@beta))^2 /
                (design@delta - design@delta_NI)^2 * sigma^2
            )
          })
