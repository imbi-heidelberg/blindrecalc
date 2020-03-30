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
setClass("TestStatistic", slots = c(
  alpha       = "numeric",
  beta        = "numeric",
  r           = "numeric",
  delta       = "numeric",
  delta_NI    = "numeric",
  alternative = "character",
  n_max       = "numeric"
))


#' Student's t test
#'
#' This class implements Student's t-test for superiority and non-inferiority
#' tests.
#' A trial with continuous outcomes of the two groups \code{T} and \code{C}
#' is assumed.
#' If \code{alternative == "greater"} the null hypothesis for the
#' mean difference
#' \ifelse{html}{\out{&Delta; = &mu;<sub>T</sub> - &mu;<sub>C</sub>}}{\eqn{\Delta = \mu_T - \mu_C}}
#' is
#' \ifelse{html}{\out{<p>H<sub>0</sub>: &Delta; &le; -&delta;<sub>NI</sub>  vs.  H<sub>1</sub>: &Delta; > -&delta;<sub>NI</sub>.</p>}}{\deqn{H_0: \Delta \leq -\delta_{NI}  vs.  H_1: \Delta > -\delta_{NI}.}}
#' Here, \ifelse{html}{\out{&delta;<sub>NI</sub> >0}}{\eqn{\delta_{NI} > 0}} denotes the non-inferiority margin.
#' If \code{alternative=="smaller"}, the direction of the effect is changed.
#'
#' @details The notation is based on the paper of Lu (2019):
#' Distribution of the two-sample t-test statistic following blinded
#' sample size re-estimation. Pharmaceutical Statistics 15: 208-215.
#'
#'
#' @aliases Student
#' @rdname Student
#' @exportClass Student
setClass("Student", contains = "TestStatistic")



#' Chi^2 test
#'
#' TODO
#'
#'
#' @aliases ChiSquare
#' @exportClass ChiSquare
setClass("ChiSquare", contains = "TestStatistic")



#' Farrington Manning test
#'
#' TODO
#'
#' @aliases FarringtonManning
#' @exportClass FarringtonManning
setClass("FarringtonManning", contains = "TestStatistic")


#' Student's t-test
#'
#' The function \code{setupStudent} creates an object of class
#' \code{\link{Student}} that can be used for sample size recalculation.
#'
#' @template setup
#' @template NI
#' @template alternative
#' @template dotdotdot
#'
#' @rdname Student
#' @export
setupStudent <- function(alpha, beta, r = 1, delta, delta_NI = 0,
                         alternative = c("greater", "smaller"), n_max = Inf, ...) {

  if (delta_NI < 0)
    stop("the non-inferiority margin must be non-negative!")

  if (all(alternative == "smaller", delta_NI != 0))
    stop("smaller alternatives are not possible for non-inferiority tests!")

  if (all(alternative == "smaller", delta > 0))
    stop("use negative effect sizes for power calculations if alternative == 'smaller'!")

  if (all(alternative == "greater", delta < 0))
    stop("use positive effect sizes for power calculations if alternative == 'greater'!")

  new("Student", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = -delta_NI, alternative = match.arg(alternative), n_max = n_max)
}



#' Setup a Chi^2 test
#'
#' This function creates an object of class \code{\link{ChiSquare-class}} that
#' can be used for sample size recalculation.
#'
#' @template setup
#' @template alternative
#' @template dotdotdot
#'
#' @details For non-inferiority trials use the function \code{\link{setupFarringtonManning}}.
#'
#' @rdname ChiSquare
#' @export
setupChiSquare <- function(alpha, beta, r = 1, delta,
                           alternative = c("greater", "smaller"), n_max = Inf, ...) {
  new("ChiSquare", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = 0, alternative = match.arg(alternative), n_max = n_max)
}




#' Setup a Farrington Manning test
#'
#' This function creates an object of class \code{\link{FarringtonManning-class}}
#' that can be used for sample size recalculation in non-inferiority trials with
#' binary endpoints.
#'
#' @template setup
#' @template NI
#' @template dotdotdot
#'
#' @rdname FarringtonManning
#' @export
setupFarringtonManning <- function(alpha, beta, r = 1, delta, delta_NI, n_max = Inf, ...) {
  if (delta_NI == 0) warning("The non-inferiority margin equals 0! Do you want to conduct a chi square test?")
  new("FarringtonManning", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = delta_NI, alternative = "greater", n_max = n_max)
}
