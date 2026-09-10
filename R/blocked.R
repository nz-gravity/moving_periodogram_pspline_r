# Residual complex-normal log likelihood, dropping the constant -N*log(pi).
tf_loglik <- function(residual, eta) -sum(eta + Mod(residual)^2 * exp(-eta))

# Diagonal tensor precision; the null directions have a proper N(0, 10^2) prior.
blocked_precision <- function(log_phi, lt, lf) {
  ifelse(lt + lf == 0, 0.01, exp(log_phi[1])*lt + exp(log_phi[2])*lf + 1e-6)
}

# Log prior in (c, log phi), including the determinant AND log-transform Jacobian.
blocked_prior <- function(c, log_phi, lt, lf) {
  q <- blocked_precision(log_phi, lt, lf)
  0.5*sum(log(q) - q*c^2) + sum(2*log_phi - exp(log_phi))
}

# A|rest ~ Normal: real and imaginary noise components each have variance S/2.
amplitude_conditional <- function(data, g, eta, prior_sd = 5) {
  precision <- 1/prior_sd^2 + 2*sum(Mod(g)^2 * exp(-eta))
  c(mean = 2*sum(Re(Conj(g)*data) * exp(-eta))/precision, sd = 1/sqrt(precision))
}

# Returns draws and POST-warmup acceptance rates for one chain. Gaussian random
# walk steps are symmetric in their chosen coordinates; the joint rescaling
# includes its Jacobian. Scales adapt only during warmup.
blocked_chain <- function(data, template, B, lt, lf, initial_theta, bounds, theta_units,
                          seed, warmup = 2000L, samples = 20000L) {
  set.seed(seed)
  K <- ncol(B)
  inside <- function(th) all(th > bounds[,1]) && all(th < bounds[,2])
  theta <- initial_theta + rnorm(2, 0, 0.005)
  if (!inside(theta)) theta <- initial_theta
  g <- template(theta)
  A <- sum(Re(Conj(g)*data))/sum(Mod(g)^2)
  log_phi <- rnorm(2, log(2), 0.1)
  q0 <- blocked_precision(log(c(2,2)), lt, lf)
  # Fixed approximate Fisher preconditioner. It changes proposals, not the target.
  covariance <- solve(crossprod(B) + diag(q0))
  L <- t(chol(covariance))
  c <- as.vector(covariance %*% crossprod(B, log(pmax(Mod(data-A*g)^2,1e-12)) + 0.5772156649))
  eta <- as.vector(B %*% c)
  Z <- B %*% L
  blocks <- split(seq_len(K), ceiling(seq_len(K)/12))
  theta_L <- diag(c(1,2.6)) %*% t(chol(matrix(c(1,-0.97,-0.97,1),2)))
  scales <- c(theta = 0.02, psd = 2.38/sqrt(12), phi = 0.5, phi_rescale = 0.2)
  accepted <- window <- numeric(4)
  names_out <- c("A","f0","fdot","phi[1]","phi[2]",sprintf("c[%d]",1:K),"lp__")
  draws <- matrix(NA_real_,samples,length(names_out),dimnames=list(NULL,names_out))
  for (i in seq_len(warmup + samples)) {
    # Block 1: exact amplitude Gibbs update, conditional on the current PSD/track.
    ap <- amplitude_conditional(data,g,eta)
    A <- rnorm(1,ap[1],ap[2])
    residual <- data-A*g
    ll <- tf_loglik(residual,eta)
    hit <- numeric(4)
    # Block 2: joint frequency/intercept and slope MH update, uniform box prior.
    candidate <- theta + scales[1] * as.vector(theta_L %*% rnorm(2))
    if (inside(candidate)) {
      gp <- template(candidate)
      lp <- tf_loglik(data-A*gp,eta)
      if (is.finite(lp) && log(runif(1)) < lp-ll) {
        theta <- candidate; g <- gp; residual <- data-A*g; ll <- lp; hit[1] <- 1
      }
    }
    # Block 3: cycle through small blocks of preconditioned PSD coordinates.
    # L and Z are fixed; each proposal is symmetric and has a full prior ratio.
    for (jj in blocks) {
      step <- scales[2]*rnorm(length(jj))
      cp <- c + as.vector(L[,jj,drop=FALSE] %*% step)
      etap <- eta + as.vector(Z[,jj,drop=FALSE] %*% step)
      lp <- tf_loglik(residual,etap)
      ratio <- lp-ll + blocked_prior(cp,log_phi,lt,lf)-blocked_prior(c,log_phi,lt,lf)
      if (is.finite(ratio) && log(runif(1)) < ratio) {
        c <- cp; eta <- etap; ll <- lp; hit[2] <- hit[2]+1/length(blocks)
      }
    }
    # Block 4: log smoothing precisions MH; additive tensor penalties are not
    # conjugate Gamma updates. The prior above includes their normalizer.
    pp <- log_phi + scales[3]*rnorm(2)
    ratio <- blocked_prior(c,pp,lt,lf)-blocked_prior(c,log_phi,lt,lf)
    if (is.finite(ratio) && log(runif(1)) < ratio) { log_phi <- pp; hit[3] <- 1 }
    # Also update phi in non-centered coordinates, holding z=sqrt(q)*c fixed.
    # Rescaling c jointly avoids slow smoothing/coefficient coupling. The
    # Gaussian-prior determinant cancels the rescaling Jacobian, leaving the
    # likelihood ratio and Gamma priors (with log-phi Jacobians).
    pp <- log_phi + scales[4]*rnorm(2)
    cp <- c*sqrt(blocked_precision(log_phi,lt,lf)/blocked_precision(pp,lt,lf))
    etap <- as.vector(B %*% cp)
    lp <- tf_loglik(residual,etap)
    ratio <- lp-ll + sum(2*pp-exp(pp))-sum(2*log_phi-exp(log_phi))
    if (is.finite(ratio) && log(runif(1)) < ratio) {
      log_phi <- pp; c <- cp; eta <- etap; ll <- lp; hit[4] <- 1
    }
    if (i <= warmup) {
      window <- window+hit
      if (i %% 50L == 0L) {
        scales <- scales * exp(0.5*(window/50-c(0.3,0.234,0.3,0.3)))
        window[] <- 0
      }
    } else {
      accepted <- accepted+hit
      draws[i-warmup,] <- c(A,theta/theta_units,exp(log_phi),c,
        ll+blocked_prior(c,log_phi,lt,lf)+dnorm(A,0,5,log=TRUE))
    }
  }
  list(draws=draws, acceptance=setNames(accepted/samples,names(scales)), scales=scales)
}
