setMethod("n_fix", signature("FarringtonManning"),
  function(design, nuisance, rounded = TRUE, ...) {
    if (design@delta_NI <= 0) stop("delta_NI has to be positive")
    p_e <- nuisance + design@delta / (1 + design@r)
    p_c <- p_e - design@delta

    pt <- p_rml(p_c, p_e, design@r, design@delta_NI)
    pt_c <- pt[1]
    pt_e <- pt[2]

    v0 <- design@r * pt_c * (1 - pt_c) + pt_e * (1 - pt_e)
    v1 <- design@r * p_c * (1 - p_c) + p_e * (1 - p_e)
    z_a <- stats::qnorm(1 - design@alpha / 2)
    z_b <- stats::qnorm(1 - design@beta)

    n <- ((1 + design@r) / design@r) * (z_a * sqrt(v0) + z_b *
        sqrt(v1))^2 / (design@delta + design@delta_NI)^2

    if (rounded) {
      n <- ceiling(n)
      if (n %% (r + 1) == 0) {
        return(n)
      } else {
        n <- n + r + 1 - n %% (r + 1)
        return(n)
      }
    } else {
      return(n)
    }
  })
