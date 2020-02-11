chisq_fix_reject <- function(design, n, nuisance, type = c("size", "power"),
                             ...) {
  reject_prob <- 0
  krit <- stats::qnorm(1 - design@alpha / 2)

  if (type == "power") {
    p1_e <- nuisance + design@delta / (1 + design@r)
    p1_c <- p1_e - design@delta

    if (p1_e > 1 | p1_e < 0 | p1_c < 0 | p1_c > 1) {
      stop("response rates outside [0, 1]")
    }
  }

  n_c <- ceiling(n / (design@r + 1))
  n_e <- ceiling(n * design@r / (design@r + 1))

  for (i in 0:n_c) {
    for (j in 0:n_e) {
      if(i + j == 0 | i + j == n_c + n_e) next
      p_c <- i / n_c
      p_e <- j / n_e
      p_hat <- (i + j) / (n_c + n_e)
      ts <- sqrt(n_c * n_e / (n_c + n_e)) * (p_e - p_c) /
        sqrt(p_hat * (1 - p_hat))

      if (design@alternative == "greater") {
        ind <- ts > krit
      } else {
        ind <- ts < -krit
      }

      if (ind & type == "size") {
        n_tot <- n_c + n_e
        reject_prob <- reject_prob + choose(n_c, i) * choose(n_e, j) *
          nuisance^(i + j) * (1 - nuisance)^(n_tot - (i + j))
      } else if (ind & type == "power") {
        reject_prob <- reject_prob + choose(n_c, i) * choose(n_e, j) *
          p1_c^i *(1 - p1_c)^(n_c - i) * p1_e^j * (1 - p1_e)^(n_e - j)
      }
    }
  }
  return(reject_prob)
}

chisq_recalc_reject <- function(design, n1, nuisance, type = c("size", "power"),
                                allocation = c("exact", "approximate"), ...) {
  reject_prob <- 0
  krit <- stats::qnorm(1 - design@alpha / 2)

  if (type == "power") {
    p1_e <- nuisance + design@delta / (1 + design@r)
    p1_c <- p1_e - design@delta
    if (p1_e > 1 | p1_e < 0 | p1_c < 0 | p1_c > 1) {
      stop("response rates outside [0, 1]")
    }
  }

  n_c1 <- ceiling(n1 / (design@r + 1))
  n_e1 <- ceiling(n1 * design@r / (design@r + 1))

  for (i in 0:n_c1) {
    for (j in 0:n_e1) {
      p_c <- i / n_c1
      p_e <- j / n_e1
      p_hat <- (i + j) / (n_c1 + n_e1)

      if (allocation == "exact") {
        n_new <- n_fix(design, p_hat, ...)
      } else {
        n_new <- ceiling(n_fix(design, p_hat,
          rounded = FALSE, ...))
      }
      n_new <- min(n_new, design@n_max)

      if ((n_new > 0 & n_new <= n1) | is.na(n_new)) {
        if (p_hat == 0 | p_hat == 1) next
        ts <- sqrt(n_c1 * n_e1 / (n_c1 + n_e1)) * (p_e - p_c) /
          sqrt(p_hat * (1 - p_hat))

        if (design@alternative == "greater") {
          ind <- ts > krit
        } else {
          ind <- ts < -krit
        }

        if (ind & type == "size") {
          reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
            nuisance^(i + j) * (1 - nuisance)^(n1 - (i + j))
        } else if (ind & type == "power") {
          reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
            p1_c^i * (1 - p1_c)^(n_c1 - i) * p1_e^j * (1 - p1_e)^(n_e1 - j)
        }
      } else {
        n2 <- n_new - (n_c1 + n_e1)
        n_c2 <- ceiling(n2 / (design@r + 1))
        n_e2 <- ceiling(n2 * design@r / (design@r + 1))

        for (k in 0:n_c2) {
          for (l in 0:n_e2) {
            n_cdot <- n_c1 + n_c2
            n_edot <- n_e1 + n_e2
            p_c <- (i + k) / n_cdot
            p_e <- (j + l) / n_edot
            p_hat <- (i + j + k + l) / (n_cdot + n_edot)
            ts <- sqrt(n_cdot * n_edot / (n_cdot + n_edot)) * (p_e - p_c) /
              sqrt(p_hat * (1 - p_hat))

            if (design@alternative == "greater") {
              ind <- ts > krit
            } else {
              ind <- ts < -krit
            }

            if (ind & type == "size") {
              x <- i + j + k + l
              n_tot <- n_c1 + n_e1 + n_c2 + n_e2
              reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
                choose(n_c2, k) * choose(n_e2, l) * nuisance^x *
                (1 - nuisance)^(n_tot - x)
            } else if (ind & type == "power") {
              x_c <- i + k
              x_e <- j + l
              reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
                choose(n_c2, k) * choose(n_e2, l) * p1_c^x_c *
                (1 - p1_c)^(n_cdot - x_c) * p1_e^x_e * (1 - p1_e)^(n_edot - x_e)
            }
          }
        }
      }
    }
  }
  return(reject_prob)
}
