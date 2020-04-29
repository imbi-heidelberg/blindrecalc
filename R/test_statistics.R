# abstract class 'TestStatistics' for internal use
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
#' Here, \ifelse{html}{\out{&delta;<sub>NI</sub> &ge; 0}}{\eqn{\delta_{NI} \geq 0}} denotes the non-inferiority margin.
#' For superiority trials, \ifelse{html}{\out{&delta;<sub>NI</sub>}}{\eqn{\delta_{NI}}}
#' can be set to zero (default).
#' If \code{alternative=="smaller"}, the direction of the effect is changed.
#'
#' @details The notation is based on the paper of Lu (2019):
#' Distribution of the two-sample t-test statistic following blinded
#' sample size re-estimation. Pharmaceutical Statistics 15: 208-215.
#'
#' @details The following methods are available for this class:
#' \code{\link{toer}}, \code{\link{pow}}, \code{\link{n_dist}},
#' \code{\link{adjusted_alpha}}, and \code{\link{n_fix}}.
#' Check the design specific documentation for details.
#'
#' @aliases Student
#' @rdname Student
#' @exportClass Student
setClass("Student", contains = "TestStatistic")



#' Chi-squared test
#'
#' This class implements a chi-squared test for superiority trials. A trial
#' with binary outcomes in two groups \code{T} and \code{C} is assumed.
#'
#' @details The following methods are available for this class:
#' \code{\link{toer}}, \code{\link{pow}}, \code{\link{n_dist}},
#' \code{\link{adjusted_alpha}}, and \code{\link{n_fix}}.
#' Check the design specific documentation for details.
#'
#' @aliases ChiSquare
#' @rdname ChiSquare
#' @exportClass ChiSquare
setClass("ChiSquare", contains = "TestStatistic")



#' Farrington Manning test
#'
#' This class implements a Farrington-Manning test for non-inferiority
#' trials. A trial with binary outcomes in two groups \code{T} and
#' \code{C} is assumed.
#'
#' @details The following methods are available for this class:
#' \code{\link{toer}}, \code{\link{pow}}, \code{\link{n_dist}},
#' \code{\link{adjusted_alpha}}, and \code{\link{n_fix}}.
#' Check the design specific documentation for details.
#'
#' @aliases FarringtonManning
#' @rdname FarringtonManning
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
#' @return An object of class \code{\link{Student}}.
#'
#' @examples
#' d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
#'                    alternative = "greater", n_max = 156)
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



#' Setup a chi-squared test
#'
#' The function \code{setupChiSquare} creates an object of class
#' \code{\link{ChiSquare}}.
#'
#' @template setup
#' @template alternative
#' @template dotdotdot
#'
#' @details For non-inferiority trials use the function \code{\link{setupFarringtonManning}}.
#'
#' @return An object of class \code{\link{ChiSquare}}.
#'
#' @examples
#' design <- setupChiSquare(alpha = .025, beta = .2, r = 1, delta = 0.2,
#' alternative = "greater")
#'
#' @rdname ChiSquare
#' @export
setupChiSquare <- function(alpha, beta, r = 1, delta,
                           alternative = c("greater", "smaller"), n_max = Inf, ...) {
  new("ChiSquare", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = 0, alternative = match.arg(alternative), n_max = n_max)
}



#' Setup a Farrington-Manning test
#'
#' The function \code{\link{setupFarringtonManning}} creates an object of
#' \code{\link{FarringtonManning}}.
#'
#' @template setup
#' @template NI
#' @template dotdotdot
#'
#' @return An object of class \code{\link{FarringtonManning}}.
#'
#' @examples
#' design <- setupFarringtonManning(alpha = .025, beta = .2, r = 1, delta = 0,
#' delta_NI = .15)
#'
#' @rdname FarringtonManning
#' @export
setupFarringtonManning <- function(alpha, beta, r = 1, delta, delta_NI, n_max = Inf, ...) {
  if (delta_NI == 0) warning("The non-inferiority margin equals 0! Do you want to conduct a chi square test?")
  new("FarringtonManning", alpha = alpha, beta = beta, r = r, delta = delta,
      delta_NI = delta_NI, alternative = "greater", n_max = n_max)
}
