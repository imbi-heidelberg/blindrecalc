#' Type I Error Rate
#'
#' Computes the type I error rate of a design with blinded sample size recalculation
#' by simulation (in the continuous case) or exactly (in the binary case)
#' for one or several values of the nuisance parameter.
#'
#' @template methods
#' @template recalculation
#' @template dotdotdot
#'
#' @return one type I error rate value for every nuisance parameter
#'
#' @export
setGeneric("toer", function(design, n1, nuisance, recalculation, ...) {
  standardGeneric("toer")
  })


#' Power
#'
#' @template methods
#' @template recalculation
#' @template dotdotdot
#'
#' @export
setGeneric("pow", function(design, n1, nuisance, recalculation, ...) {
  standardGeneric("pow")
})


#' Distribution of the Sample Size
#'
#' @template methods
#' @template recalculation
#' @template dotdotdot
#'
#' @export
setGeneric("sample_size_dist", function(design, n1, nuisance, recalculation,  ...) {
  standardGeneric("sample_size_dist")
})



#' Adjusted level of significance
#'
#' This method returns an adjusted significance level that can be used
#' such that the actual type I error rate is protected.
#'
#' @template methods
#' @param tol desired absolute tolerance
#'
#' @export
setGeneric("adjusted_alpha", function(design, n1, nuisance, ...) {
  standardGeneric("adjusted_alpha")
})




#' Fixed Sample Size
#'
#' Returns the sample size for the corresponding one-stage design without
#' sample size recalculation.
#'
#' @param s test statistic object
#' @param nuisance nuisance parameter
#' @template dotdotdot
#'
#' @export
setGeneric("n_fix", function(design, nuisance, ...) {
  standardGeneric("n_fix")
})
