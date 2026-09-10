# Base-R numerical checks: no Python or Stan installation needed.
source("R/simulate.R")
source("R/moving_periodogram.R")
source("R/pspline.R")
set.seed(21)
x <- rnorm(200)
m <- 8L
ord <- moving_periodogram(x, m, thin = 2L)
# Independent FFT definition checks phase, normalization and window location.
for (r in seq_len(nrow(ord))) {
  center <- round(ord$u[r] * length(x))
  j <- round(ord$omega[r] * (2*m + 1)/(2*pi))
  value <- fft(x[(center-m):(center+m)])[j+1] / sqrt(2*pi*(2*m+1))
  stopifnot(abs(value - complex(real = ord$re[r], imaginary = ord$im[r])) < 1e-12)
}
stopifnot(inherits(try(moving_periodogram(1:10, 8), silent = TRUE), "try-error"))
# White-noise expectation and LS2 variance normalization.
means <- replicate(100, mean(moving_periodogram(rnorm(1024), 16)$mi))
stopifnot(abs(mean(means) - 1/(2*pi)) < 0.015)
u <- 0.37; om <- seq(0, pi, length.out = 10001)
a <- 1.1*cos(1.5-cos(4*pi*u))
s <- true_psd_ls2(u, om)
stopifnot(abs(2 * sum((head(s,-1)+tail(s,-1))/2) * pi/10000 - (1+a^2)) < 1e-10)
sp <- spline_setup(8)
B <- basis(sp, seq(0, 1, length.out = 101))
stopifnot(max(abs(rowSums(B)-1)) < 1e-12, min(B) > -1e-12,
          sum(sp$values == 0) == 2L, abs(sum(diag(sp$penalty))-1) < 1e-12)
# Independent numerical integration of the roughness quadratic form.
w <- rnorm(ncol(B)); xx <- seq(0, 1, length.out = 100001)
dd <- splines::splineDesign(sp$knots, xx, ord = 4, derivs = 2)
rawP <- crossprod(dd) / 100000
rawP <- rawP - (tcrossprod(dd[1,])+tcrossprod(dd[nrow(dd),]))/200000
stopifnot(max(abs(rawP/sum(diag(rawP))-sp$penalty)) < 1e-7)
# Eigenbasis and scattered surface use the same column-major coefficient order.
p <- prepare_model(ord, 3, 4)
C <- matrix(rnorm(p$data$Kt*p$data$Kf), p$data$Kt)
W <- p$time$vectors %*% C %*% t(p$frequency$vectors)
y1 <- rowSums((basis(p$time, ord$u) %*% W) * basis(p$frequency, ord$omega/pi))
y2 <- vapply(seq_len(nrow(ord)), function(r)
  sum(p$data$Bt[r,] * (C %*% p$data$Bf[p$data$rung[r],])), numeric(1))
stopifnot(max(abs(y1-y2)) < 1e-12)
cat("All R numerical checks passed.\n")
