#' @param design Object of class \code{ChiSquare} created by \code{setupChiSquare}.
#' @param nuisance Value of the nuisance parameter. For the
#'   Chi-Squared test this is the overall response rate.
#' @param variance A character string indicating whether the "\code{heterogenous}" (default)
#'   or the "\code{homogeneous}" variance formula should be used.
#' @param rounded Whether the calculated sample size should be rounded up such that
#'   the allocation ratio is preserved.
#' @template dotdotdot
#'
#' @rdname ChiSquare
#' @export
setMethod("n_fix", signature("ChiSquare"),
  function(design, nuisance, variance = c("heterogeneous", "homogeneous"),
           rounded = TRUE, ...) {
    variance <- match.arg(variance)
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
            nuisance * (1 - nuisance)) + z_b * sqrt(design@r * p_c *
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

#' @template methods_chisquare
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @rdname ChiSquare
#' @export
setMethod("toer", signature("ChiSquare"),
  function(design, n1, nuisance, recalculation,
           allocation = c("exact", "approximate"), ...) {
    allocation <- match.arg(allocation)
    if (allocation == "exact") {
      if (sum(n1 %% (design@r + 1) != 0) > 0) {
        stop("no integer sample sizes")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("no integer sample sizes for n_max")
      }
    }
    if (sum(nuisance < 0) + sum(nuisance > 1) != 0) {
      stop("nuisance has to be within [0, 1]")
    }
    if (sum(design@n_max < n1) > 0) {
      stop("n_max is smaller than n1")
    }

    if ((length(n1) > 1) & (length(nuisance) > 1)) {
      stop("only one of n1 and nuisance can have length > 1")
    } else if (length(n1) > 1) {
      if (recalculation) {
        nmat <- lapply(n1, function(x) get_nmat_chisq(design, x, allocation, ...))
        mapply(chisq_recalc_reject, n1 = n1, nmat = nmat,
          MoreArgs = list(design = design, nuisance = nuisance, type = "size"))
      } else {
        sapply(n1, function(x) chisq_fix_reject(design, x, nuisance, "size"))
      }
    } else if (length(nuisance) > 1) {
      if (recalculation) {
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        sapply(nuisance, function(x) chisq_recalc_reject(design, n1, x, "size", nmat))
      } else {
        sapply(nuisance, function(x) chisq_fix_reject(design, n1, x, "size"))
      }
    } else {
      if (recalculation) {
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        chisq_recalc_reject(design, n1, nuisance, "size", nmat)
      } else {
        chisq_fix_reject(design, n1, nuisance, "size")
      }
    }
  })

#' @template methods_chisquare
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @rdname ChiSquare
#' @export
setMethod("pow", signature("ChiSquare"),
  function(design, n1, nuisance, recalculation,
           allocation = c("exact", "approximate"), ...) {
    allocation <- match.arg(allocation)
    if (allocation == "exact") {
      if (n1 %% (design@r + 1) != 0) {
        stop("no integer sample sizes for first stage")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("no integer sample sizes for n_max")
      }
    }
    if (sum(nuisance < 0) + sum(nuisance > 1) > 0) {
      stop("nuisance has to be within [0, 1]")
    }
    if (sum(design@n_max < n1) > 0) {
      stop("n_max is smaller than n1")
    }

    if ((length(n1) > 1) & (length(nuisance) > 1)) {
      stop("only one of n1 and nuisance can have length > 1")
    } else if (length(n1) > 1) {
      if (recalculation) {
        nmat <- lapply(n1, function(x) get_nmat_chisq(design, x, allocation, ...))
        mapply(chisq_recalc_reject, n1 = n1, nmat = nmat,
          MoreArgs = list(design = design, nuisance = nuisance, type = "power"))
      } else {
        sapply(n1, function(x) chisq_fix_reject(design, x, nuisance, "power"))
      }
    } else if (length(nuisance) > 1) {
      if (recalculation) {
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        sapply(nuisance, function(x) chisq_recalc_reject(design, n1, x, "power", nmat))
      } else {
        sapply(nuisance, function(x) chisq_fix_reject(design, n1, x, "power"))
      }
    } else {
      if (recalculation) {
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        chisq_recalc_reject(design, n1, nuisance, "power", nmat)
      } else {
        chisq_fix_reject(design, n1, nuisance, "power")
      }
    }


  })

#' @template methods_chisquare
#' @template adjalpha_binary
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @rdname ChiSquare
#' @export
setMethod("adjusted_alpha", signature("ChiSquare"),
  function(design, n1, nuisance, nuis_ass, precision = 0.001, gamma = 0,
           recalculation, allocation = c("exact", "approximate"), ...) {
    allocation <- match.arg(allocation)
    if (allocation == "exact") {
      if (n1 %% (design@r + 1) != 0) {
        stop("no integer sample sizes for first stage")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("no integer sample sizes for n_max")
      }
    }
    if (sum(nuisance < 0) + sum(nuisance >1) > 0) {
      stop("nuisance has to be within [0, 1]")
    }

    alpha_nom <- design@alpha - gamma
    if (recalculation) {
      repeat {
        nmat <- get_nmat_chisq(design, n1, allocation, ...)
        alpha_max <- max(sapply(nuisance,
          function(x) chisq_recalc_reject(design, n1, x, "size", nmat)))
        if (alpha_max <= alpha_nom) break
        design@alpha <- design@alpha - precision
      }
    } else {
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


#' @template methods_chisquare
#' @param summary Logical. If \code{TRUE} (default) a summary of the sample
#'   size distribution is printed. If \code{FALSE} all sample sizes are
#'   printed.
#' @template plot
#' @template allocation
#' @template dotdotdot
#'
#' @rdname ChiSquare
#' @export
setMethod("n_dist", signature("ChiSquare"),
  function(design, n1, nuisance, summary = TRUE, plot = FALSE,
           allocation = c("exact", "approximate"), ...) {
    allocation <- match.arg(allocation)
    if (allocation == "exact") {
      if (sum(n1 %% (design@r + 1) != 0) > 0) {
        stop("no integer sample sizes for first stage")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("no integer sample sizes for n_max")
      }
    }
    if (sum(nuisance < 0) + sum(nuisance > 1) > 0) {
      stop("nuisance has to be within [0, 1]")
    }

    if (((length(n1) > 1) & (length(nuisance) > 1)) |
        ((length(n1) == 1) & (length(nuisance) == 1))) {
      stop("one of n1 and nuisance must have length > 1")
    } else if (length(nuisance) > 1) {
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
