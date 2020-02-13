#' Sample Size Calculation for the Chi-Squared Test
#'
#' Calculates the sample size for the chi-squared test for the corresponding
#' one-stage design without sample size recalculation.
#'
#' @param design an object of class \code{ChiSquare} created by \code{setupChiSquare()}.
#' @param nuisance the overall response rate.
#' @param variance a character string indicating whether the "\code{heterogenous}" (default)
#' or the "\code{homogeneous}" variance formula should be used.
#'
#' @return
#' @export
#'
#' @examples
setMethod("n_fix", signature("ChiSquare"),
  function(design, nuisance, variance = c("heterogeneous", "homogeneous"),
           rounded = TRUE, ...) {
    variance <- match.arg(variance)
    if (length(nuisance) > 1) {
      sapply(nuisance, function(x) n_fix(design = design, nuisance = x,
        variance = variance))
    } else {
      p_e <- nuisance + design@delta / (1 + design@r)
      p_c <- p_e - design@delta
      z_a <- stats::qnorm(1 - design@alpha / 2)
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

#' Calculation of the Actual Level of the Chi-Squared Test
#'
#' Calculation of the actual level of the chi-squared test for the fixed sample design and the
#' internal pilot study design.
#'
#' @param n1 Either the total sample size (if \code{design} is \code{"fixed"}) or
#' sample size of the first stage (if \code{design} is \code{"ips"})
#' @param nuisance the overall response rate.
#' @param recalculation
#' @param allocation
#'
#' @return
#' @export
#'
#' @examples
setMethod("toer", signature("ChiSquare"),
  function(design, n1, nuisance, recalculation,
           allocation = c("exact", "approximate"), ...) {
    allocation <- match.arg(allocation)
    if (allocation == "exact") {
      if (n1 %% (design@r + 1) != 0) {
        stop("no integer sample sizes")
      }
      if (is.finite(design@n_max) & design@n_max %% (design@r + 1) != 0) {
        stop("no integer sample sizes for n_max")
      }
    }
    if (nuisance < 0 | nuisance > 1) {
      stop("nuisance has to be within [0, 1]")
    }
    if (design@n_max < n1) {
      stop("n_max is smaller than n1")
    }

    if (recalculation) {
      nmat <- get_nmat_chisq(design, n1, allocation, ...)
      chisq_recalc_reject(design, n1, nuisance, "size", nmat)
    } else {
      chisq_fix_reject(design, n1, nuisance, "size")
    }
  })

#' Calculation of the Power of the Chi-Squared Test
#'
#' Calculation of the power of the chi-squared test for the fixed sample design and the
#' internal pilot study design.
#'
#' @param design
#' @param n1 Either the total sample size (if \code{design} is \code{"fixed"}) or
#' sample size of the first stage (if \code{design} is \code{"ips"})
#' @param nuisance the overall response rate.
#' @param recalculation
#' @param allocation
#'
#' @return
#' @export
#'
#' @examples
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
    if (nuisance < 0 | nuisance > 1) {
      stop("nuisance has to be within [0, 1]")
    }

    if (recalculation) {
      nmat <- get_nmat_chisq(design, n1, allocation, ...)
      chisq_recalc_reject(design, n1, nuisance, "power", nmat)
    } else {
      chisq_fix_reject(design, n1, nuisance, "power")
    }
  })
