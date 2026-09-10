# Frequency-cycling moving periodogram: no padding, whole blocks only.
# x: finite real vector; m: half-window length; thin: positive block spacing.
# Returns u (normalized window-center time), omega (radians/sample), mi (power),
# re/im (normalized complex Fourier coefficient parts); one row per ordinate.
moving_periodogram <- function(x, m = 16L, thin = 2L) {
  stopifnot(is.numeric(x), is.null(dim(x)), all(is.finite(x)),
            length(m) == 1L, m >= 1L, m == as.integer(m),
            length(thin) == 1L, thin >= 1L, thin == as.integer(thin))
  n <- length(x)
  blocks <- (n - 2L * m) %/% (thin * m)
  if (blocks < 1L) stop("Series too short for m and thin")
  start <- rep(thin * m * (seq_len(blocks) - 1L), each = m) +
    rep(0:(m - 1L), blocks)
  j <- rep(seq_len(m), blocks)
  omega <- 2 * pi * seq_len(m) / (2 * m + 1)
  phase <- exp(-1i * outer(0:(2 * m), omega))
  coeff <- vapply(seq_along(start), function(r)
    sum(x[start[r] + seq_len(2 * m + 1L)] * phase[, j[r]]), complex(1)) /
    sqrt(2 * pi * (2 * m + 1))
  data.frame(u = (start + m + 1) / n, omega = omega[j],
             mi = Mod(coeff)^2, re = Re(coeff), im = Im(coeff))
}
