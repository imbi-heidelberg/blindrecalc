#' Type I Error Rate
#'
#' Computes the type I error rate of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods
#' @template recalculation
#' @template dotdotdot
#'
#' @return One type I error rate value for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @details The method is implemented for the classes \code{\link{Student}},
#' \code{\link{ChiSquare}}, and \code{\link{FarringtonManning}}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' toer(d, n1 = 20, nuisance = 5.5, recalculation = TRUE)
#'
#' @export
setGeneric("toer", function(design, n1, nuisance, recalculation, ...) {
  standardGeneric("toer")
  })


#' Power
#'
#' Calculates the power of designs with blinded sample size recalculation
#' or of fixed designs for one or several values of the nuisance parameter.
#'
#' @template methods
#' @template recalculation
#' @template dotdotdot
#'
#' @return One power value for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @details The method is implemented for the classes \code{\link{Student}},
#' \code{\link{ChiSquare}}, and \code{\link{FarringtonManning}}.
#'
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' pow(d, n1 = 20, nuisance = 5.5, recalculation = TRUE)
#'
#'
#' @export
setGeneric("pow", function(design, n1, nuisance, recalculation, ...) {
  standardGeneric("pow")
})


#' Distribution of the Sample Size
#'
#' Calculates the distribution of the total sample sizes of designs
#' with blinded sample size recalculation for different values of the
#' nuisance parameter or of n1.
#'
#' @template methods
#' @param summary logical - is a summary of the sample size distribution desired?
#'    Otherwise, a vector with sample sizes is returned.
#' @template plot
#' @template dotdotdot
#'
#' @return Summary and/or plot of the sample size distribution for
#'   every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @details The method is implemented for the classes \code{\link{Student}},
#' \code{\link{ChiSquare}}, and \code{\link{FarringtonManning}}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' n_dist(d, n1 = 20, nuisance = 5.5, summary = TRUE, plot = FALSE, seed = 2020)
#'
#'
#' @export
setGeneric("n_dist", function(design, n1, nuisance, summary = TRUE, plot = FALSE, ...) {
  standardGeneric("n_dist")
})


#' Adjusted level of significance
#'
#' This method returns an adjusted significance level that can be used
#' such that the actual type I error rate is preserved.
#'
#' @template methods
#' @template dotdotdot
#'
#' @return Value of the adjusted significance level
#'  for every nuisance parameter and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @details The method is implemented for the classes \code{\link{Student}},
#' \code{\link{ChiSquare}}, and \code{\link{FarringtonManning}}.
#' Check the class-specific documentation for further parameters that have
#' to be specified.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 0, delta_NI = 1.5, n_max = 848)
#' sigma <- c(2, 5.5, 9)
#' adjusted_alpha(design = d, n1 = 20, nuisance = sigma, tol = 1e-4, iters = 1e3)
#'
#' @export
setGeneric("adjusted_alpha", function(design, n1, nuisance, ...) {
  standardGeneric("adjusted_alpha")
})


#' Fixed Sample Size
#'
#' Returns the total sample size of a fixed design without sample size recalculation.
#'
#' @param design test statistic object created by \code{setup}
#' @param nuisance nuisance parameter for the respective test problem
#' @template dotdotdot
#'
#' @return One value of the fixed sample size for every nuisance parameter
#'  and every value of n1.
#'
#' @details The method is only vectorized in either \code{nuisance}
#'   or \code{n1}.
#'
#' @details The method is implemented for the classes \code{\link{Student}},
#' \code{\link{ChiSquare}}, and \code{\link{FarringtonManning}}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                   alternative = "greater", n_max = 156)
#' n_fix(design = d, nuisance = 5.5)
#'
#' @export
setGeneric("n_fix", function(design, nuisance, ...) {
  standardGeneric("n_fix")
})
