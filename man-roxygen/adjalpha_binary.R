#' @param nuis_ass If \code{recalculation = FALSE} this is the value for
#'   the overall response rate that is used to calculate the sample size
#'   for the adjusted significance level.
#' @param precision Value by which the nominal type 1 error rate is
#'   reduced in each iteration until the nominal type 1 error rate is
#'   preserved.
#' @param gamma If \code{gamma > 0} then the significance level is adjusted
#'   such that the actual level is at most \code{alpha - gamma}. This is
#'   necessary to maintain the nomininal significance level if a confidence
#'   interval approach proposed by Friede & Kieser (2011) is used.
