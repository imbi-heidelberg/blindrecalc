#' Test statistics
#'
#' In \code{blindrecalc}, different test statistics are implemented.
#' Currently, those are Student's t-test for superiority and non-inferiorty
#' for continuous outcomes, the Chi^2-test for superiority tests for binary
#' outcomes and the Farrington Manning test for non-inferiority tests for
#' binary outcomes.
#'
#' @aliases TestStatistic
#' @exportClass TestStatistic
setClass("TestStatistic", representation(
  alpha    = "numeric",
  beta     = "numeric",
  r        = "numeric",
  delta    = "numeric",
  delta_NI = "numeric",
  n_max    = "numeric"
))


#' Student's t test
#'
#' TODO
#'
#' @aliases Student
#' @exportClass Student
setClass("Student", contains = "TestStatistic")

setClass("ChiSquare", contains = "TestStatistic")

setClass("FarringtonManning", contains = "TestStatistic")


#' Setup Student's t-test
#'
#' This function creates an object of class \code{Student} that can
#' be used for sample size recalculation.
#'
#' @template setup
#' @template NI
#' @template dotdotdot
#'
#' @export
setupStudent <- function(alpha, beta, r = 1, delta, delta_NI = 0, n_max = Inf, ...) {
  new("Student", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = delta_NI, n_max = n_max)
}



#' Setup a Chi^2 test
#'
#' This function creates an object of class \code{ChiSquare} that can
#' be used for sample size recalculation.
#'
#' @template setup
#' @template dotdotdot
#'
#' @details For non-inferiority trials use the function \code{\link{setupFarringtionManning}}.
#'
#' @export
setupChiSquare <- function(alpha, beta, r = 1, delta, n_max = Inf, ...) {
  new("ChiSquare", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = delta_NI, n_max = n_max)
}




#' Setup a Farrington Manning test
#'
#' This function creates an object of class \code{FarringtonManning} that can
#' be used for sample size recalculation in non-inferiority trials with
#' binary endpoints.
#'
#' @template setup
#' @template NI
#' @template dotdotdot
#'
#' @export
setupFarringtonManning <- function(alpha, beta, r = 1, delta, delta_NI, n_max = Inf, ...) {
  if (delta_NI == 0) warning("The non-inferiority margin equals 0! Do you want to conduct a chi square test?")
  new("FarringtonManning", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = delta_NI, n_max = n_max)
}
