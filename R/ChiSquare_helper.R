# Helper function that creates a matrix with the total sample sizes
# for all possible outcomes of the internal pilot study
get_nmat_chisq <- function(design, n1, allocation, ...) {
  n_c1 <- ceiling(n1 / (design@r + 1))
  n_e1 <- ceiling(n1 * design@r / (design@r + 1))
  p_hatfun <- function(i, j, n_c1, n_e1) {
    (i + j) / (n_c1 + n_e1)
  }
  out.mat <- expand.grid(i = 0:n_c1, j = 0:n_e1)

  out.mat$p_hat <- mapply(p_hatfun, i = out.mat$i, j = out.mat$j,
    MoreArgs = list(n_c1 = n_c1, n_e1 = n_e1))

  if (allocation == "exact") {
    out.mat$n <- n_fix(design, nuisance = out.mat$p_hat, ...)
  } else if (allocation == "kf_approx") {
    out.mat$n <- ceiling(n_fix(design, nuisance = out.mat$p_hat, rounded = FALSE, ...))
  } else {
    out.mat$n <- n_fix(design, nuisance = out.mat$p_hat, rounded = FALSE, ...)
  }
  out.mat$n <- pmin(out.mat$n, design@n_max)
  out.mat$n <- ifelse(is.na(out.mat$n), -99, out.mat$n)
  return(as.matrix(out.mat))
}

# Helper function that calculates the probability for all possible
# total sample sizes based on n1 and nuisance
n_distrib_chisq <- function(design, n1, nuisance, allocation, ...) {
  p_e <- nuisance + design@delta / (1 + design@r)
  p_c <- p_e - design@delta
  n_c1 <- ceiling(n1 / (design@r + 1))
  n_e1 <- ceiling(n1 - n_c1)

  n_new <- prob <- numeric()

  for (i in 0:n_c1) {
    for (j in 0:n_e1) {
      if (i + j == 0 | i +j == n1) next
      p1 <- i / n_c1
      p2 <- j / n_e1
      p_hat <- (i + j) / (n_c1 + n_e1)
      if (allocation == "exact") {
        n_new <- c(n_new, n_fix(design, nuisance = p_hat, ...))
      } else {
        n_new <- c(n_new, ceiling(n_fix(design, nuisance = p_hat, rounded = FALSE, ...)))
      }
      prob <- c(prob, choose(n_c1, i) * choose(n_e1, j) * p_c^i * (1 - p_c)^(n_c1 - i) *
          p_e^j * (1 - p_e)^(n_e1 - j))
    }
  }

  n_new[which(is.na(n_new))] <- n1
  n_new[which(n_new < n1)] <- n1
  out <- stats::aggregate(prob, list(n_new), sum)
  colnames(out) <- c("n", "prob")
  return(out)
}
