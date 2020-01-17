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
  function(design, nuisance, variance = c("heterogeneous", "homogeneous"), ...) {
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
          return(n)
        } else if (variance == "homogeneous") {
          n <- (1 + design@r)^2 / design@r * (z_a + z_b)^2 / design@delta^2 *
            nuisance * (1 - nuisance)
          return(n)
        }
      }
    }
  })


setMethod("toer", signature("ChiSquare"),
  function(design, n1, p, recalculation = TRUE, ...) {
    size <- 0
    krit <- stats::qnorm(1 - design@alpha / 2)

    if (recalculation) {
      n_c1 <- ceiling(n1 / (design@r + 1))
      n_e1 <- ceiling(n1 * design@r / (design@r + 1))
      for (i in 0:n_c1) {
        for (j in 0:n_e1) {
          p_c <- i / n_c1
          p_e <- j / n_e1
          p_hat <- (p_c + design@r * p_e) / (1 + design@r)
          n_new <- getn_chisq(p_hat, design@delta, design@r, design@alpha, 1 - design@beta)
          n_new <- min(n_new, design@n_max)

          if ((n_new > 0 & n_new <= n) | is.na(n_new)) {
            if (p_hat == 0 | p_hat == 1) next
            ts <- sqrt(n_c1 * n_e1 / (n_c1 + n_e1)) * (p_e - p_c) /
              sqrt(p_hat * (1 - p_hat))
            ind <- ts > krit
            if (ind) {
              size <- size + choose(n_c1, i) * choose(n_e1, j) * p^(i + j) *
                (1 - p)^(n - (i + j))
            }
          } else {
            n2 <- n_new - (n_c1 + n_e1)
            n_c2 <- ceiling(n2 / (r + 1))
            n_e2 <- ceiling(n2 * r / (r + 1))

            for (k in 0:n_c2) {
              for (l in 0:n_e2) {
                n_cdot <- n_c1 + n_c2
                n_edot <- n_e1 + n_e2
                p_c <- (i + k) / n_cdot
                p_e <- (j + l) / n_edot
                p_hat <- (p_c + r * p_e) / (1 + r)
                ts <- sqrt(n_cdot * n_edot / (n_cdot + n_edot)) * (p_e - p_c) /
                  sqrt(p_hat * (1 - p_hat))
                ind <- ts > krit
                if (ind) {
                  x <- i + j + k + l
                  size <- size + choose(n_c1, i) * choose(n_e1, j) * choose(n_c2, k) *
                    choose(n_e2, l) * p^x * (1 - p)^(n_new - x)
                }
              }
            }
          }
        }
      }
      return(size)
  }
})
