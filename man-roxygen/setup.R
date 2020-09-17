#' @param alpha One-sided type I error rate.
#' @param beta Type II error rate.
#' @param r Allocation ratio between experimental and control group.
#' @param delta Difference of effect size between alternative and null hypothesis.
#' @param n_max Maximal overall sample size. If the recalculated sample size
#'   is greater than \code{n_max} it is set to \code{n_max}.
#'   The default value equals \code{n_max = Inf}.
