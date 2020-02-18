#' Sample Size Calculation for the Farrington-Manning Test
#'
#' Calculates the sample size for the Farrington-Manning test for the corresponding
#' one-stage design without sample size recalculation.
#'
#' @param design an object of class \code{FarringtonManning} created
#' by \code{setupFarringtonManning()}.
#' @param nuisance the overall response rate.
#' @param rounded

#' @export
#'
#' @examples
setMethod("n_fix", signature("FarringtonManning"),
  function(design, nuisance, rounded = TRUE, ...) {
    if (design@delta_NI <= 0) stop("delta_NI has to be positive")
    p_e <- nuisance + design@delta / (1 + design@r)
    p_c <- p_e - design@delta

    pt <- p_rml(p_c, p_e, design@r, design@delta_NI)
    pt_c <- pt[1]
    pt_e <- pt[2]

    v0 <- design@r * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)
    v1 <- design@r * p_c * (1 - p_c) + p_e * (1 - p_e)
    z_a <- stats::qnorm(1 - design@alpha / 2)
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
  })

#' Calculation of the Actual Level of the Chi-Squared Test
#'
#' Calculation of the actual level of the chi-squared test for the fixed sample design and the
#' internal pilot study design.
#'
#' @param design an object of class \code{FarringtonManning} created
#'    by \code{setupFarringtonManning()}.
#' @param n1 Either the total sample size (if \code{design} is \code{"fixed"}) or
#'    sample size of the first stage (if \code{design} is \code{"ips"})
#' @param nuisance the overall response rate.
#' @template recalculation
#' @tempalte allocation
#' @template dotdotdot
#'
#' @export
#'
#' @examples
setMethod("toer", signature("FarringtonManning"),
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
        nmat <- lapply(n1, function(x) get_nmat_fm(design, x, allocation, ...))
        mapply(fm_recalc_reject, n1 = n1, nmat = nmat,
          MoreArgs = list(design = design, nuisance = nuisance, type = "size"))
      } else {
        sapply(n1, function(x) fm_fix_reject(design, x, nuisance, "size"))
      }
    } else if (length(nuisance) > 1) {
      if (recalculation) {
        nmat <- get_nmat_fm(design, n1, allocation, ...)
        sapply(nuisance, function(x) fm_recalc_reject(design, n1, x, "size", nmat))
      } else {
        sapply(nuisance, function(x) fm_fix_reject(design, n1, x, "size"))
      }
    } else {
      if (recalculation) {
        nmat <- get_nmat_fm(design, n1, allocation, ...)
        fm_recalc_reject(design, n1, nuisance, "size", nmat)
      } else {
        fm_fix_reject(design, n1, nuisance, "size")
      }
    }
  })

#' Calculation of the Power of the Farrington-Manning test
#'
#' Calculation of the power of the Farrington-Manning test for the fixed sample design
#' and the internal pilot study design.
#'
#' @param design
#' @param n1 Either the total sample size (if \code{design} is \code{"fixed"}) or
#' sample size of the first stage (if \code{design} is \code{"ips"})
#' @param nuisance the overall response rate.
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @export
#'
#' @examples
setMethod("pow", signature("FarringtonManning"),
function(design, n1, nuisance, recalculation,
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
  if (sum(design@n_max < n1) > 0) {
    stop("n_max is smaller than n1")
  }

  if ((length(n1) > 1) & (length(nuisance) > 1)) {
    stop("only one of n1 and nuisance can have length > 1")
  } else if (length(n1) > 1) {
    if (recalculation) {
      nmat <- lapply(n1, function(x) get_nmat_fm(design, x, allocation, ...))
      mapply(fm_recalc_reject, n1 = n1, nmat = nmat,
        MoreArgs = list(design = design, nuisance = nuisance, type = "power"))
    } else {
      sapply(n1, function(x) fm_fix_reject(design, x, nuisance, "power"))
    }
  } else if (length(nuisance) > 1) {
    if (recalculation) {
      nmat <- get_nmat_fm(design, n1, allocation, ...)
      sapply(nuisance, function(x) fm_recalc_reject(design, n1, x, "power", nmat))
    } else {
      sapply(nuisance, function(x) fm_fix_reject(design, n1, x, "power"))
    }
  } else {
    if (recalculation) {
      nmat <- get_nmat_fm(design, n1, allocation, ...)
      fm_recalc_reject(design, n1, nuisance, "power", nmat)
    } else {
      fm_fix_reject(design, n1, nuisance, "power")
    }
  }
})

#' Calculation of the Adjusted Alpha for the Farrington-Manning
#'
#' Calculates the adjusted alpha that is necessary to maintain the nominimal type 1
#' error rate for the fixed sample design and the internal pilot study design in the
#' Farrington-Manning test.
#'
#' @param design
#' @param n1 Either the total sample size or sample size of the first stage
#' @param nuisance A vector of nuisance parameters
#' @param precision
#' @template recalculation
#' @template allocation
#' @template dotdotdot
#'
#' @export
setMethod("adjusted_alpha", signature("FarringtonManning"),
  function(design, n1, nuisance, precision = 0.001, recalculation,
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
    if (sum(nuisance < 0) + sum(nuisance >1) > 0) {
      stop("nuisance has to be within [0, 1]")
    }

    alpha_nom <- design@alpha / 2
    if (recalculation) {
      repeat {
        nmat <- get_nmat_fm(design, n1, allocation, ...)
        alpha_max <- max(sapply(nuisance,
          function(x) fm_recalc_reject(design, n1, x, "size", nmat)))
        if (alpha_max <= alpha_nom) break
        design@alpha <- design@alpha - precision
      }
    } else {
      repeat {
        alpha_max <- max(sapply(nuisance,
          function(x) fm_fix_reject(design, n1, x, "size")))
        if (alpha_max <= alpha_nom) break
        design@alpha <- design@alpha - precision
      }
    }
    return(design@alpha)
  })
