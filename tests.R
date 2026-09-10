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

# Blocked sampler: compare complex likelihood and log-phi prior ratios with
# independent normalized R densities; verify the amplitude conditional by quadrature.
source("R/blocked.R")
z <- c(1+2i,-0.3+0.2i); eta <- c(-0.5,0.4)
expected <- sum(dnorm(Re(z),0,sqrt(exp(eta)/2),log=TRUE) +
                dnorm(Im(z),0,sqrt(exp(eta)/2),log=TRUE))
stopifnot(abs(tf_loglik(z,eta) - expected - length(z)*log(pi)) < 1e-12)
lt <- c(0,0,1,1); lf <- c(0,1,0,1); cc <- c(.2,-.5,.1,.8)
prior_reference <- function(lp) {
  q <- ifelse(lt+lf==0,.01,exp(lp[1])*lt+exp(lp[2])*lf+1e-6)
  sum(dnorm(cc,0,1/sqrt(q),log=TRUE)) + sum(dgamma(exp(lp),2,1,log=TRUE)+lp)
}
p1 <- log(c(2,3)); p2 <- log(c(4,1))
stopifnot(abs((blocked_prior(cc,p1,lt,lf)-blocked_prior(cc,p2,lt,lf)) -
              (prior_reference(p1)-prior_reference(p2))) < 1e-12)
g <- c(.4+.3i,.8-.2i); ap <- amplitude_conditional(z,g,eta)
peak <- tf_loglik(z-ap[1]*g,eta)+dnorm(ap[1],0,5,log=TRUE)
f <- function(a) vapply(a,function(v) exp(tf_loglik(z-v*g,eta)+dnorm(v,0,5,log=TRUE)-peak),numeric(1))
norm <- integrate(f,-Inf,Inf)$value
mu <- integrate(function(a) a*f(a),-Inf,Inf)$value/norm
variance <- integrate(function(a) (a-mu)^2*f(a),-Inf,Inf)$value/norm
stopifnot(abs(mu-ap[1])<1e-6,abs(variance-ap[2]^2)<1e-6)
cat("Blocked likelihood, prior and Gibbs conditional checks passed.\n")


# The cached transform is the same linear operator as moving_periodogram().
source("R/chirp.R")
set.seed(17); xx <- rnorm(128); oo <- moving_periodogram(xx,m=8)
op <- fourier_operator(128,8,round(oo$u*128),oo$omega)
stopifnot(max(Mod(op(xx)-complex(real=oo$re,imaginary=oo$im))) < 1e-12)
hh <- chirp_time(128,c(.12,.0002))
stopifnot(max(Mod(op(xx+hh)-op(xx)-op(hh))) < 1e-12)
# Independent dense covariance check of the exact LS2 matched-filter norm.
nn <- 32; aa <- 1.1*cos(1.5-cos(4*pi*(0:(nn-1))/nn))
covariance <- diag(1+aa^2)
for (i in 2:nn) covariance[i,i-1] <- covariance[i-1,i] <- aa[i]
v <- rnorm(nn)
stopifnot(abs(sum(ls2_whiten(v)^2)-sum(v*solve(covariance,v))) < 1e-9)
cat("Chirp transform and exact LS2 SNR checks passed.\n")
# The joint non-centered phi/coefficient move includes the change-of-variables
# Jacobian, cancelling the Gaussian normalization in the simplified MH ratio.
q1 <- blocked_precision(p1,lt,lf); q2 <- blocked_precision(p2,lt,lf)
cp <- cc*sqrt(q1/q2)
full <- blocked_prior(cp,p2,lt,lf)-blocked_prior(cc,p1,lt,lf) + .5*sum(log(q1/q2))
stopifnot(abs(full-(sum(2*p2-exp(p2))-sum(2*p1-exp(p1)))) < 1e-12)
cat("Joint precision/coefficient MH Jacobian check passed.\n")
