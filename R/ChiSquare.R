#' Type I Error Rate
#'
#' Computes the type I error rate of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods_chisquare
#' @template recalculation
#' @template allocation_chisquare
#' @template dotdotdot
#'
#' @return One type I error rate value for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#'   d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
#'   toer(d, n1 = c(10, 20), nuisance = 0.25, recalculation = TRUE)
#'
#' @rdname toer.ChiSquare
#' @export
setMethod("toer", signature("ChiSquare"),
  function(design, n1, nuisance, recalculation,
           allocation = c("exact", "approximate", "kf_approx"), ...) {
    allocation <- match.arg(allocation)
    # Check if input is valid
    if (allocation == "exact") {
      if (sum(n1 %% (design@r + 1) != 0) > 0) {
        stop("No integer sample sizes.")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("No integer sample sizes for n_max.")
      }
    }
    if (sum(nuisance < 0) + sum(nuisance > 1) != 0) {
      stop("Nuisance has to be within [0, 1].")
    }
    if (sum(design@n_max < n1) > 0) {
      stop("n_max is smaller than n1.")
    }
    # Check whether n1 or nuisance is a vector
    if ((length(n1) > 1) & (length(nuisance) > 1)) {
      stop("only one of n1 and nuisance can have length > 1")
    } else if (length(n1) > 1) {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- lapply(n1, function(x) get_nmat_chisq(design, x, allocation, ...))
        mapply(chisq_recalc_reject, n1 = n1, nmat = nmat,
          MoreArgs = list(design = design, nuisance = nuisance, type = "size"))
      } else {
        sapply(n1, function(x) chisq_fix_reject(design, x, nuisance, "size"))
      }
    } else if (length(nuisance) > 1) {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        sapply(nuisance, function(x) chisq_recalc_reject(design, n1, x, "size", nmat))
      } else {
        sapply(nuisance, function(x) chisq_fix_reject(design, n1, x, "size"))
      }
    } else {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        chisq_recalc_reject(design, n1, nuisance, "size", nmat)
      } else {
        chisq_fix_reject(design, n1, nuisance, "size")
      }
    }
  })



#' Power
#'
#' Calculates the power of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods_chisquare
#' @template recalculation
#' @template allocation_chisquare
#' @template dotdotdot
#'
#' @return One power value for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#'   d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
#'   pow(d, n1 = 20, nuisance = c(0.2, 0.4), recalculation = TRUE)
#'
#' @rdname pow.ChiSquare
#' @export
setMethod("pow", signature("ChiSquare"),
  function(design, n1, nuisance, recalculation,
           allocation = c("exact", "approximate", "kf_approx"), ...) {
    allocation <- match.arg(allocation)
    # Check if input is valid
    if (allocation == "exact") {
      if (n1 %% (design@r + 1) != 0) {
        stop("No integer sample sizes for first stage.")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("No integer sample sizes for n_max.")
      }
    }
    if (sum(nuisance < 0) + sum(nuisance > 1) > 0) {
      stop("Nuisance has to be within [0, 1].")
    }
    if (sum(design@n_max < n1) > 0) {
      stop("n_max is smaller than n1.")
    }

    # Check whether n1 or nuisance is a vector
    if ((length(n1) > 1) & (length(nuisance) > 1)) {
      stop("Only one of n1 and nuisance can have length > 1")
    } else if (length(n1) > 1) {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- lapply(n1, function(x) get_nmat_chisq(design, x, allocation, ...))
        mapply(chisq_recalc_reject, n1 = n1, nmat = nmat,
          MoreArgs = list(design = design, nuisance = nuisance, type = "power"))
      } else {
        sapply(n1, function(x) chisq_fix_reject(design, x, nuisance, "power"))
      }
    } else if (length(nuisance) > 1) {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        sapply(nuisance, function(x) chisq_recalc_reject(design, n1, x, "power", nmat))
      } else {
        sapply(nuisance, function(x) chisq_fix_reject(design, n1, x, "power"))
      }
    } else {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        chisq_recalc_reject(design, n1, nuisance, "power", nmat)
      } else {
        chisq_fix_reject(design, n1, nuisance, "power")
      }
    }
  })



#' Distribution of the Sample Size
#'
#' Calculates the distribution of the total sample sizes of designs
#' with blinded sample size recalculation for different values of the
#' nuisance parameter or of n1.
#'
#' @template methods_chisquare
#' @template summary
#' @template plot
#' @template allocation_chisquare
#' @template dotdotdot
#'
#' @details Only sample sizes that occur with a probability of at least 0.01% are
#' considered.
#'
#' @return Summary and/or plot of the sample size distribution for
#'   every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#'   d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
#'   n_dist(d, n1 = 20, nuisance = 0.25, summary = TRUE, plot = FALSE)
#'
#' @rdname n_dist.ChiSquare
#' @export
setMethod("n_dist", signature("ChiSquare"),
          function(design, n1, nuisance, summary = TRUE, plot = FALSE,
                   allocation = c("exact", "approximate"), ...) {
            allocation <- match.arg(allocation)
            # Check if input is valid
            if (allocation == "exact") {
              if (sum(n1 %% (design@r + 1) != 0) > 0) {
                stop("No integer sample sizes for first stage.")
              }
              if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
                stop("No integer sample sizes for n_max.")
              }
            }
            if (sum(nuisance < 0) + sum(nuisance > 1) > 0) {
              stop("Nuisance has to be within [0, 1].")
            }

            # Check whether n1 or nuisance is a vector
            if ((length(n1) > 1) & (length(nuisance) > 1)) {
              stop("Only one of n1 and nuisance can have length > 1.")
            } else if (length(n1) == 1) {
              # Calculate possible total sample sizes and probabilities
              out <- lapply(nuisance, function(x) n_distrib_chisq(design, n1, x, allocation, ...))
              out <- Map(cbind, out, nuisance = nuisance)
              out <- do.call("rbind", out)
              out <- with(out, data.frame(n = rep(n, prob * 10000), p = rep(nuisance, prob * 10000)))
              out.list <- split(out$n, paste0("p = ", out$p))

              if (plot) {
                graphics::boxplot(out.list, ...)
              }
              if (summary) {
                sapply(out.list, summary)
              } else {
                out.list
              }
            } else {
              out <- lapply(n1, function(x) n_distrib_chisq(design, x, nuisance, allocation, ...))
              out <- Map(cbind, out, n1 = n1)
              out <- do.call("rbind", out)
              out <- with(out, data.frame(n = rep(n, prob * 10000), n1 = rep(n1, prob * 10000)))
              out.list <- split(out$n, paste0("n1 = ", out$n1))

              if (plot) {
                graphics::boxplot(out.list, ...)
              }
              if (summary) {
                sapply(out.list, summary)
              } else {
                out.list
              }
            }
          })




#' Adjusted level of significance
#'
#' This method returns an adjusted significance level that can be used
#' such that the actual type I error rate is preserved.
#'
#' @template methods_chisquare
#' @template adjalpha_binary
#' @template recalculation
#' @template allocation_chisquare
#' @template dotdotdot
#'
#' @return Value of the adjusted significance level for every
#'  nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#'   d <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
#'   adjusted_alpha(d, n1 = 10, nuisance = 0.3, gamma = 0.001,
#'      nuis_ass = 0.3, precision = 0.001, recalculation = TRUE)
#'
#' @rdname adjusted_alpha.ChiSquare
#' @export
setMethod("adjusted_alpha", signature("ChiSquare"),
  function(design, n1, nuisance, nuis_ass, precision = 0.001, gamma = 0,
           recalculation, allocation = c("exact", "approximate"), ...) {
    allocation <- match.arg(allocation)
    # Check if input is valid
    if (allocation == "exact") {
      if (n1 %% (design@r + 1) != 0) {
        stop("No integer sample sizes for first stage.")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("No integer sample sizes for n_max.")
      }
    }
    if (sum(nuisance < 0) + sum(nuisance >1) > 0) {
      stop("Nuisance has to be within [0, 1].")
    }

    alpha_nom <- design@alpha - gamma
    if (recalculation) {
      # iteratively reduce the significance level until it is sufficiently small
      repeat {
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        alpha_max <- max(sapply(nuisance,
          function(x) chisq_recalc_reject(design, n1, x, "size", nmat)))
        if (alpha_max <= alpha_nom) break
        design@alpha <- design@alpha - precision
      }
    } else {
      # iteratively reduce the significance level until it is sufficiently small
      repeat {
        alpha_max <- max(sapply(nuisance,
          function(x) chisq_fix_reject(design, n1, x, "size")))
        if (alpha_max <= alpha_nom) break
        design@alpha <- design@alpha - precision
        if (allocation == "exact") {
          n1 <- n_fix(design, nuis_ass, ...)
        } else {
          n1 <- n_fix(design, nuis_ass, rounded = FALSE, ...)
        }
      }
    }
    return(design@alpha)
  })



#' Fixed Sample Size
#'
#' Returns the sample size of a fixed design without sample size recalculation.
#'
#' @param design Object of class \code{ChiSquare} created by \code{setupChiSquare}.
#' @param nuisance Value of the nuisance parameter. For the
#'   Chi-Squared test this is the overall response rate.
#' @param variance A character string indicating whether the "\code{heterogenous}" (default)
#'   or the "\code{homogeneous}" variance formula should be used.
#' @param rounded Whether the calculated sample size should be rounded up such that
#'   the allocation ratio is preserved.
#' @template dotdotdot
#'
#' @return One value of the fixed sample size for every nuisance parameter
#'  and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#'   design1 <- setupChiSquare(alpha = 0.025, beta = 0.2, r = 1, delta = 0.2)
#'   n_fix(design1, nuisance = c(0.2, 0.3))
#'
#' @rdname n_fix.ChiSquare
#' @export
setMethod("n_fix", signature("ChiSquare"),
          function(design, nuisance, variance = c("heterogeneous", "homogeneous"),
                   rounded = TRUE, ...) {
            variance <- match.arg(variance)
            # Check if input is valid
            if (sum(nuisance < 0) + sum(nuisance > 1) > 0) {
              stop("Nuisance has to be within [0, 1].")
            }
            # Use recursion if nuisance is a vector
            if (length(nuisance) > 1) {
              sapply(nuisance, function(x) n_fix(design = design, nuisance = x,
                                                 variance = variance, rounded = rounded, ...))
            } else {
              p_e <- nuisance + design@delta / (1 + design@r)
              p_c <- p_e - design@delta
              z_a <- stats::qnorm(1 - design@alpha)
              z_b <- stats::qnorm(1 - design@beta)

              if (p_e < 0 | p_c < 0 | p_e > 1 | p_c > 1) {
                return(NA)
              } else {
                if (variance == "heterogeneous") {
                  n <- (1 + design@r) / design@r * (z_a * sqrt((1 + design@r) *
                                                                 nuisance * (1 - nuisance)) +
                                                      z_b * sqrt(design@r * p_c *
                                                                   (1 - p_c) + p_e * (1 - p_e)))^2 /design@delta^2
                } else if (variance == "homogeneous") {
                  n <- (1 + design@r)^2 / design@r * (z_a + z_b)^2 / design@delta^2 *
                    nuisance * (1 - nuisance)
                }
                if (rounded) {
                  n <- ceiling(n)
                  if (n %% (design@r + 1) == 0) {
                    return(n)
                  } else {
                    n <- n + design@r + 1 - n %% (design@r + 1)
                    return(n)
                  }
                } else {
                  return(n)
                }
              }
            }
          })
