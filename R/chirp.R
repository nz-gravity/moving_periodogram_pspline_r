# Unit-amplitude linear chirp: dt=1 second, known phase zero at t=0.
# theta=(f0 [Hz], fdot [Hz/s]); instantaneous frequency=f0+fdot*t.
chirp_time <- function(n, theta) {
  t <- 0:(n-1L)
  cos(2*pi*(theta[1]*t + 0.5*theta[2]*t^2))
}

# Cached local Fourier operator at arbitrary paired centers/frequencies.
# Each row uses a full 2m+1 window; normalization matches moving_periodogram().
fourier_operator <- function(n, m, centers, omega) {
  stopifnot(length(centers)==length(omega), all(centers-m>=1), all(centers+m<=n))
  index <- outer(centers-m-1L, seq_len(2*m+1L), "+")
  weights <- exp(-1i*outer(omega,0:(2*m))) / sqrt(2*pi*(2*m+1))
  function(x) {
    stopifnot(length(x)==n)
    rowSums(weights * matrix(x[index],nrow=length(centers)))
  }
}

# Exact whitening for LS2's tridiagonal covariance, including the first sample's
# unobserved preceding innovation. Thus sum(ls2_whiten(h)^2) = h' C^{-1} h.
ls2_whiten <- function(x) {
  n <- length(x); u <- (0:(n-1))/n
  a <- 1.1*cos(1.5-cos(4*pi*u))
  diagonal <- sqrt(1+a[1]^2)
  y <- numeric(n); y[1] <- x[1]/diagonal
  for (i in 2:n) {
    sub <- a[i]/diagonal
    diagonal <- sqrt(1+a[i]^2-sub^2)
    y[i] <- (x[i]-sub*y[i-1])/diagonal
  }
  y
}

# Data-only initialization: scan chirp rates and FFT the dechirped time series.
# This uses the full prior search box, not the injection values; it is not a
# calibrated detection statistic. Local optimization later includes known phase.
chirp_initialize <- function(x, bounds) {
  n <- length(x); t <- 0:(n-1)
  nfft <- 8*n
  bins <- which((0:(nfft-1))/nfft > bounds[1,1] & (0:(nfft-1))/nfft < bounds[1,2])
  best <- -Inf; theta <- numeric(2)
  rates <- seq(bounds[2,1],bounds[2,2],length.out=ceiling(diff(bounds[2,])*4*n^2)+1)
  for (rate in rates) {
    z <- fft(c(x*exp(-1i*pi*rate*t^2),rep(0,nfft-n)))
    j <- bins[which.max(Mod(z[bins]))]
    if (Mod(z[j]) > best) { best <- Mod(z[j]); theta <- c((j-1)/nfft,rate) }
  }
  theta
}
