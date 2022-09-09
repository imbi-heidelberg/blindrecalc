#' Blinded Sample Size Recalculation
#'
#' The package \pkg{blindrecalc} provides characteristics and plots
#' of trial designs with blinded sample size recalculation where
#' a nuisance parameter is estimated at an blinded interim analysis.
#'
#' Currently, for continuous outcomes, a t-test is implemented for superiority
#' and non-inferiority trials.
#' For superiority trials with binary endpoint, the chi^2-test is implemented.
#' The Farrington Manning test covers non-inferiority trials with binary endpoint.
#'
#' @useDynLib blindrecalc
#' @import methods
#' @importFrom Rcpp sourceCpp
#' @docType package
#' @name blindrecalc
NULL
