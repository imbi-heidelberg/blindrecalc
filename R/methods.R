#' Type I error rate
#'
#' Computes the type I error rate of a design with blinded sample size recalculation
#' by simulation (in the continuous case) or exactly (in the binary case)
#' for one or several values of the nuisance parameter.
#'
#' @template methods
#'
#' @return one type I error rate value for every nuisance parameter
#'
#' @export
setGeneric("toer", function(s, n1, nuisance, ...) standardGeneric("toer"))




#' Power
#'
#' @template methods
#'
#' @export
setGeneric("pow", function(s, n1, nuisance, ...) standardGeneric("pow"))


#' Distribution of the sample size
#'
#' @template methods
#'
#' @export
setGeneric("sample_size_dist", function(s, n1, nuisance, ...) standardGeneric("sample_size_dist"))


#' Fixed sample size
#'
#' Returns the sample size for the corresponding one-stage design without
#' sample size recalculation.
#'
#' @param s test statistic object
#' @param nuisance nuisance parameter
#'
#' @export
setGeneric("n_fix", function(s, nuisance, ...) standardGeneric("n_fix"))
