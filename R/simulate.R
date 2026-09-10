# LS2 simulator with zero-based rescaled time u = 0, ..., (n-1)/n.
# n: sample count; optional innovations: at least n+1 standard-normal values.
# Returns a numeric vector of length n. Set the random seed before calling.
simulate_ls2 <- function(n, innovations = rnorm(n + 2L)) {
  stopifnot(n >= 2L, n == as.integer(n), length(innovations) >= n + 1L)
  u <- (seq_len(n) - 1) / n
  a <- 1.1 * cos(1.5 - cos(4 * pi * u))
  innovations[2:(n + 1)] + a * innovations[seq_len(n)]
}

# Two-sided density per radian: white noise has S = 1/(2*pi).
# u and omega: paired coordinates (or one scalar); omega is in radians/sample.
# Returns a density vector, not an outer-product grid; use expand.grid for a grid.
true_psd_ls2 <- function(u, omega) {
  a <- 1.1 * cos(1.5 - cos(4 * pi * u))
  (1 + a^2 + 2 * a * cos(omega)) / (2 * pi)
}
