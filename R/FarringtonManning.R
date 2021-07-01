#' Type I Error Rate
#'
#' Computes the type I error rate of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods_fm
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @return One type I error rate value for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.2)
#' toer(d, n1 = 20, nuisance = 0.25, recalculation = TRUE, allocation = "approximate")
#'
#' @rdname toer.FarringtonManning
#' @export
setMethod("toer", signature("FarringtonManning"),
  function(design, n1, nuisance, recalculation,
           allocation = c("exact", "approximate"), ...) {
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
    if (sum(nuisance < 0) + sum(nuisance > 1) > 0) {
      stop("Nuisance has to be within [0, 1].")
    }
    if (sum(design@n_max < n1) > 0) {
      stop("n_max is smaller than n1.")
    }

    # Check whether n1 or nuisance is a vector
    if ((length(n1) > 1) & (length(nuisance) > 1)) {
      stop("Only one of n1 and nuisance can have length > 1.")
    } else if (length(n1) > 1) {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- lapply(n1, function(x) get_nmat_fm(design, x, allocation, ...))
        mapply(fm_recalc_reject, n1 = n1, nmat = nmat,
          MoreArgs = list(design = design, nuisance = nuisance, type = "size"))
      } else {
        sapply(n1, function(x) fm_fix_reject(design, x, nuisance, "size"))
      }
    } else if (length(nuisance) > 1) {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- get_nmat_fm(design, n1, allocation, ...)
        sapply(nuisance, function(x) fm_recalc_reject(design, n1, x, "size", nmat))
      } else {
        sapply(nuisance, function(x) fm_fix_reject(design, n1, x, "size"))
      }
    } else {
      if (recalculation) {
        # Create matrix with total sample sizes
        nmat <- get_nmat_fm(design, n1, allocation, ...)
        fm_recalc_reject(design, n1, nuisance, "size", nmat)
      } else {
        fm_fix_reject(design, n1, nuisance, "size")
      }
    }
  })



#' Power
#'
#' Calculates the power of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods_fm
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @return One power value for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.25)
#' pow(d, n1 = 30, nuisance = 0.4, allocation = "approximate", recalculation = TRUE)
#'
#' @rdname pow.FarringtonManning
#' @export
setMethod("pow", signature("FarringtonManning"),
function(design, n1, nuisance, recalculation,
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
  if (sum(design@n_max < n1) > 0) {
    stop("n_max is smaller than n1.")
  }

  # Check whether n1 or nuisance is a vector
  if ((length(n1) > 1) & (length(nuisance) > 1)) {
    stop("only one of n1 and nuisance can have length > 1")
  } else if (length(n1) > 1) {
    if (recalculation) {
      # Create matrix with total sample sizes
      nmat <- lapply(n1, function(x) get_nmat_fm(design, x, allocation, ...))
      mapply(fm_recalc_reject, n1 = n1, nmat = nmat,
        MoreArgs = list(design = design, nuisance = nuisance, type = "power"))
    } else {
      sapply(n1, function(x) fm_fix_reject(design, x, nuisance, "power"))
    }
  } else if (length(nuisance) > 1) {
    if (recalculation) {
      # Create matrix with total sample sizes
      nmat <- get_nmat_fm(design, n1, allocation, ...)
      sapply(nuisance, function(x) fm_recalc_reject(design, n1, x, "power", nmat))
    } else {
      sapply(nuisance, function(x) fm_fix_reject(design, n1, x, "power"))
    }
  } else {
    if (recalculation) {
      # Create matrix with total sample sizes
      nmat <- get_nmat_fm(design, n1, allocation, ...)
      fm_recalc_reject(design, n1, nuisance, "power", nmat)
    } else {
      fm_fix_reject(design, n1, nuisance, "power")
    }
  }
})



#' Distribution of the Sample Size
#'
#' Calculates the distribution of the total sample sizes of designs
#' with blinded sample size recalculation for different values of the
#' nuisance parameter or of n1.
#'
#' @template methods_fm
#' @template summary
#' @template plot
#' @template allocation
#' @template dotdotdot
#'
#' @details Only sample sizes that occur with a probability of at least 0.01% are
#' considered.
#'
#' @return Summary and/or plot of the sample size distribution for
#'   each nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.25)
#' n_dist(d, n1 = 30, nuisance = 0.2, summary = TRUE, plot = FALSE)
#'
#' @rdname n_dist.FarringtonManning
#' @export
setMethod("n_dist", signature("FarringtonManning"),
          function(design, n1, nuisance, summary, plot,
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
              # Calculate possible sample sizes and probabilities
              out <- lapply(nuisance, function(x) n_distrib_fm(design, n1, x, allocation, ...))
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
              # Calculate possible total sample sizes and probabilities
              out <- lapply(n1, function(x) n_distrib_fm(design, x, nuisance, allocation, ...))
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
#' @template methods_fm
#' @template adjalpha_binary
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @return Value of the adjusted significance level for every nuisance
#'  parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @examples
#' d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.25)
#' adjusted_alpha(d, n1 = 20, nuisance = 0.5, recalculation = TRUE)
#'
#' @rdname adjusted_alpha.FarringtonManning
#' @export
setMethod("adjusted_alpha", signature("FarringtonManning"),
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
        nmat <- get_nmat_fm(design, n1, allocation, ...)
        alpha_max <- max(sapply(nuisance,
          function(x) fm_recalc_reject(design, n1, x, "size", nmat)))
        if (alpha_max <= alpha_nom) break
        design@alpha <- design@alpha - precision
      }
    } else {
      repeat {
        # iteratively reduce the significance level until it is sufficiently small
        alpha_max <- max(sapply(nuisance,
          function(x) fm_fix_reject(design, n1, x, "size")))
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
#' @param design Object of class \code{FarringtonManning} created
#'   by \code{setupFarringtonManning}.
#' @param nuisance Value of the nuisance parameter. For the
#'   Farrington-Manning test this is the overall response rate.
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
#' d <- setupFarringtonManning(alpha = 0.025, beta = 0.2, r = 1, delta = 0, delta_NI = 0.25)
#' n_fix(d, nuisance = 0.3)
#'
#' @rdname n_fix.FarringtonManning
#' @export
setMethod("n_fix", signature("FarringtonManning"),
          function(design, nuisance, rounded = TRUE, ...) {
            # Check if input is valid
            if (design@delta_NI <= 0) {
              stop("delta_NI has to be positive.")
            }
            if (sum(nuisance < 0) + sum(nuisance > 1) > 0) {
              stop("Nuisance has to be within [0, 1].")
            }
            # Use recursion if nuisance is a vector
            if (length(nuisance) > 1) {
              sapply(nuisance, function(x) n_fix(design = design, nuisance = x,
                                                 rounded = rounded, ...))
            } else {
              p_e <- nuisance + design@delta / (1 + design@r)
              p_c <- p_e - design@delta

              pt <- p_rml(p_c, p_e, design@r, design@delta_NI)
              pt_c <- pt[1]
              pt_e <- pt[2]

              v0 <- design@r * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)
              v1 <- design@r * p_c * (1 - p_c) + p_e * (1 - p_e)
              z_a <- stats::qnorm(1 - design@alpha)
              z_b <- stats::qnorm(1 - design@beta)

              n <- ((1 + design@r) / design@r) * (z_a * sqrt(v0) + z_b *
                sqrt(v1))^2 / (design@delta + design@delta_NI)^2

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
          })
