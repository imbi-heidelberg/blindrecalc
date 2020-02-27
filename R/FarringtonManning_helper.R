get_nmat_fm <- function(design, n1, allocation, ...) {
  n_c1 <- ceiling(n1 / (design@r + 1))
  n_e1 <- ceiling(n1 * design@r / (design@r + 1))
  p_hatfun <- function(i, j, n_c1, n_e1) {
    (i + j) / (n_c1 + n_e1)
  }
  out.mat <- expand.grid(i = 0:n_c1, j = 0:n_e1)

  out.mat$p_hat <- mapply(p_hatfun, i = out.mat$i, j = out.mat$j,
    MoreArgs = list(n_c1 = n_c1, n_e1 = n_e1))

  if (allocation == "exact") {
    out.mat$n <- sapply(out.mat$p_hat, function(x) n_fix(design, nuisance = x, ...))
  } else {
    out.mat$n <- sapply(out.mat$p_hat, function(x) n_fix(design, nuisance = x,
      rounded = FALSE, ...))
  }
  out.mat$n <- pmin(out.mat$n, design@n_max)
  out.mat$n <- ifelse(is.na(out.mat$n), -99, out.mat$n)
  return(as.matrix(out.mat))
}

n_distrib_fm <- function(design, n1, nuisance, allocation, ...) {
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

# p_rml1 <- function(p_c, p_e, r, margin) {
#   a <- 1 + (1 / r)
#   b <- -(1 + (1 / r) + p_e + (1 / r) * p_c - margin * ((1 / r) + 2))
#   c <- margin^2 - margin * (2 * p_e + (1 / r) + 1) + p_e + (1 / r) * p_c
#   d <- p_e * margin * (1 - margin)
#
#   v <- b^3 / (3 * a)^3 - (b * c) / (6 * a^2) + d / (2 * a)
#   u <- sign(v) * sqrt(b^2 / (3 * a)^2 - c / (3 * a))
#   x <- ifelse(v == 0 & u == 0, 0,
#     ifelse(v / u^3 > 1, 1, v / u^3))
#   w <- (1 / 3) * (pi + acos(x))
#
#   pt_e <- max(0, 2 * u * cos(w) - b / (3 * a))
#   pt_c <- min(1, pt_e + margin)
#
#   return(c(pt_c, pt_e))
# }

# fm_fix_reject1 <- function(design, n, nuisance, type = c("size", "power"),
#                           ...) {
#
#   p0_e <- nuisance - design@delta_NI / (1 + design@r)
#   p0_c <- p0_e + design@delta_NI
#   if (p0_e < 0 | p0_c > 1) {
#     stop("combination of nuisance and delta_NI implies probabilities outside [0, 1]")
#   }
#
#   if (type == "power") {
#     p1_e <- nuisance - design@delta / (1 + design@r)
#     p1_c <- p1_e + design@delta
#     if (p1_e < 0 | p1_c > 1) {
#       stop("combination of nuisance and delta implies probabilities outside [0, 1]")
#     }
#   }
#
#   reject_prob <- 0
#   krit <- stats::qnorm(1 - design@alpha / 2)
#   n_c <- ceiling(n / (design@r + 1))
#   n_e <- ceiling(n * design@r / (design@r + 1))
#   for (i in 0:n_c) {
#     for (j in 0:n_e) {
#       if(i + j == 0 | i + j == n_c + n_e) next
#       p_c <- i / n_c
#       p_e <- j / n_e
#       pt <- p_rml(p_c, p_e, design@r, design@delta_NI)
#       pt_c <- pt[1]
#       pt_e <- pt[2]
#
#       se <- sqrt((design@r * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)) / n_e)
#       ts <- (p_e - p_c + design@delta_NI) / se
#
#       ind <- ts > krit
#       if (ind & type == "size") {
#         reject_prob <- reject_prob + choose(n_c, i) * choose(n_e, j) * p0_c^i *
#           (1 - p0_c)^(n_c - i) * p0_e^j * (1 - p0_e)^(n_e - j)
#       } else if (ind & type == "power") {
#         reject_prob <- reject_prob + choose(n_c, i) * choose(n_e, j) * p1_c^i *
#           (1 - p1_c)^(n_c - i) * p1_e^j * (1 - p1_e)^(n_e - j)
#       }
#     }
#   }
#   return(reject_prob)
# }

# fm_recalc_reject <- function(design, n1, nuisance, type = c("size", "power"),
#                              allocation = c("exact", "approximate"), ...) {
#   p0_e <- nuisance - design@delta_NI / (1 + design@r)
#   p0_c <- p0_e + design@delta_NI
#   if (p0_e < 0 | p0_c > 1) {
#     stop("combination of nuisance and delta_NI implies probabilities outside [0, 1]")
#   }
#
#   if (type == "power") {
#     p1_e <- nuisance - design@delta / (1 + design@r)
#     p1_c <- p1_e + design@delta
#     if (p1_e < 0 | p1_c > 1) {
#       stop("combination of nuisance and delta implies probabilities outside [0, 1]")
#     }
#   }
#
#   reject_prob <- 0
#   krit <- stats::qnorm(1 - design@alpha / 2)
#   n_c1 <- ceiling(n1 / (design@r + 1))
#   n_e1 <- ceiling(n1 * design@r / (design@r + 1))
#
#   for (i in 0:n_c1) {
#     for (j in 0:n_e1) {
#       p_c <- i / n_c1
#       p_e <- j / n_e1
#       p_hat <- (i + j) / (n_c1 + n_e1)
#
#       if (allocation == "exact") {
#         n_new <- n_fix(design, p_hat, ...)
#       } else {
#         n_new <- ceiling(n_fix(design, p_hat,
#           rounded = FALSE, ...))
#       }
#       n_new <- min(n_new, design@n_max)
#
#       if ((n_new > 0 & n_new <= n1) | is.na(n_new)) {
#         pt <- p_rml(p_c, p_e, design@r, design@delta_NI)
#         pt_c <- pt[1]
#         pt_e <- pt[2]
#         se <- sqrt((design@r * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)) / n_e1)
#         ts <- (p_e - p_c + design@delta_NI) / se
#         ind <- ts > krit
#         if (ind & type == "size") {
#           reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
#             p0_c^i * (1 - p0_c)^(n_c1 - i) * p0_e^j * (1 - p0_e)^(n_e1 - j)
#         } else if (ind & type == "power") {
# 		  reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
# 			p1_c^i * (1 - p1_c)^(n_c1 - i) * p1_e^j * (1 - p1_e)^(n_e1 - j)
#         }
#       } else {
#         n2 <- n_new - (n_c1 + n_e1)
#         n_c2 <- ceiling(n2 / (design@r + 1))
#         n_e2 <- ceiling(n2 * design@r / (design@r + 1))
#
#         for (k in 0:n_c2) {
#           for (l in 0:n_e2) {
#             n_cdot <- n_c1 + n_c2
#             n_edot <- n_e1 + n_e2
#             p_c <- (i + k) / n_cdot
#             p_e <- (j + l) / n_edot
#             pt <- p_rml(p_c, p_e, design@r, design@delta_NI)
#             pt_c <- pt[1]
#             pt_e <- pt[2]
#             se <- sqrt((design@r * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)) / n_edot)
#             ts <- (p_e - p_c + design@delta_NI) / se
#             ind <- ts > krit
#             if (ind & type == "size") {
#               reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
#                 choose(n_c2, k) * choose(n_e2, l) * p0_c^(i + k) *
#                 (1 - p0_c)^(n_cdot - (i + k)) * p0_e^(j + l) *
#                 (1 - p0_e)^(n_edot - (j + l))
#             } else if (ind & type == "power") {
#               reject_prob <- reject_prob + choose(n_c1, i) * choose(n_e1, j) *
#                 choose(n_c2, k) * choose(n_e2, l) * p1_c^(i + k) *
#                 (1 - p1_c)^(n_cdot - (i + k)) * p1_e^(j + l) *
#                 (1 - p1_e)^(n_edot - (j + l))
#             }
#           }
#         }
#       }
#     }
#   }
#   return(reject_prob)
# }
