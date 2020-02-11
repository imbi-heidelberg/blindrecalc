p_rml <- function(p_c, p_e, r, margin) {
  a <- 1 + (1 / r)
  b <- -(1 + (1 / r) + p_e + (1 / r) * p_c - margin * ((1 / r) + 2))
  c <- margin^2 - margin * (2 * p_e + (1 / r) + 1) + p_e + (1 / r) * p_c
  d <- p_e * margin * (1 - margin)

  v <- b^3 / (3 * a)^3 - (b * c) / (6 * a^2) + d / (2 * a)
  u <- sign(v) * sqrt(b^2 / (3 * a)^2 - c / (3 * a))
  x <- ifelse(v == 0 & u == 0, 0,
    ifelse(v / u^3 > 1, 1, v / u^3))
  w <- (1 / 3) * (pi + acos(x))

  pt_e <- max(0, 2 * u * cos(w) - b / (3 * a))
  pt_c <- min(1, pt_e + margin)

  return(c(pt_c, pt_e))
}
