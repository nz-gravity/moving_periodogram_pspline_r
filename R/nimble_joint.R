# Built-in NUTS for the PSD, then built-in block Metropolis for the signal.
# Input: prepare_joint() result. Output: draws, sampler diagnostics and timings.
fit_nimble_joint <- function(d, warmup=1000L, samples=3000L, seed=501L) {
  if (!requireNamespace("nimbleHMC",quietly=TRUE))
    stop('Install dependencies: install.packages(c("nimbleHMC", "posterior"))')
  library(nimbleHMC)
  K <- ncol(d$B); N <- length(d$data); L <- 2*d$m+1L
  centers <- round(d$ordinates$u*d$n)
  weights <- exp(-1i*outer(d$ordinates$omega,0:(L-1)))/sqrt(2*pi*L)
  constants <- list(T=d$n,N=N,K=K,L=L,t=(0:(d$n-1))/d$n,
    start=centers-d$m,end=centers+d$m,wr=Re(weights),wi=Im(weights),
    B=d$B,lt=d$lt,lf=d$lf,null=as.numeric(d$lt+d$lf==0),
    lower=d$bounds[,1],upper=d$bounds[,2])
  code <- nimble::nimbleCode({
    A ~ dnorm(0,sd=5)
    theta[1] ~ dunif(lower[1],upper[1])
    theta[2] ~ dunif(lower[2],upper[2])
    f0 <- theta[1]/T
    fdot <- theta[2]/(T*T)
    for(j in 1:2) phi[j] ~ dgamma(2,rate=1)
    for(k in 1:K) {
      q[k] <- null[k]*0.01+(1-null[k])*(phi[1]*lt[k]+phi[2]*lf[k]+1e-6)
      c[k] ~ dnorm(0,tau=q[k])
    }
    wave[1:T] <- A*cos(2*3.141592653589793*(theta[1]*t[1:T]+0.5*theta[2]*t[1:T]^2))
    for(i in 1:N) {
      eta[i] <- inprod(B[i,1:K],c[1:K])
      sigma[i] <- exp(eta[i]/2)/sqrt(2)
      mur[i] <- inprod(wr[i,1:L],wave[start[i]:end[i]])
      mui[i] <- inprod(wi[i,1:L],wave[start[i]:end[i]])
      dre[i] ~ dnorm(mur[i],sd=sigma[i])
      dimag[i] ~ dnorm(mui[i],sd=sigma[i])
    }
  })
  g <- d$template(d$initial)
  A0 <- sum(Re(Conj(g)*d$data))/sum(Mod(g)^2)
  q0 <- ifelse(d$lt+d$lf==0,.01,2*d$lt+2*d$lf+1e-6)
  c0 <- as.vector(solve(crossprod(d$B)+diag(q0),
    crossprod(d$B,log(pmax(Mod(d$data-A0*g)^2,1e-12))+.5772156649)))
  inits <- list(A=A0,theta=d$initial,c=c0,phi=c(2,2))
  started <- proc.time()[[3]]
  model <- nimble::nimbleModel(code,constants=constants,
    data=list(dre=Re(d$data),dimag=Im(d$data)),inits=inits,buildDerivs=TRUE)
  # Independent likelihood check, including the complex-normal constant.
  expected <- -sum(d$B %*% c0+Mod(d$data-A0*g)^2*exp(-d$B %*% c0))-N*log(pi)
  actual <- model$calculate(model$getNodeNames(dataOnly=TRUE))
  stopifnot(abs(actual-expected)<1e-7)
  conf <- nimble::configureMCMC(model,nodes=NULL,monitors=c("A","f0","fdot","c","phi"))
  conf$addSampler(target=c("c","phi"),type="NUTS",
    control=list(warmupMode="iterations",warmup=warmup,delta=.9))
  # Scale and covariance are in (amplitude, n*f0, n^2*fdot) coordinates.
  propCov <- diag(c(.05,.03,.08)) %*%
    matrix(c(1,0,0,0,1,-.97,0,-.97,1),3) %*% diag(c(.05,.03,.08))
  conf$addSampler(target=c("A","theta"),type="RW_block",control=list(propCov=propCov))
  mcmc <- nimble::buildMCMC(conf)
  compiled <- nimble::compileNimble(model,mcmc)
  compile_seconds <- proc.time()[[3]]-started
  started <- proc.time()[[3]]
  sampler <- compiled$mcmc$samplerFunctions[[1]]
  signal_sampler <- compiled$mcmc$samplerFunctions[[2]]
  diagnostics <- vector("list",4)
  chains <- lapply(1:4,function(ch) {
    message("NIMBLE chain ",ch,"/4")
    set.seed(seed+ch-1L)
    compiled$model$setInits(list(A=A0+rnorm(1,0,.03),theta=d$initial+rnorm(2,0,.005),
      c=c0+rnorm(K,0,.01),phi=exp(rnorm(2,log(2),.1))))
    signal_sampler$adaptive <- TRUE
    compiled$mcmc$run(warmup,nburnin=warmup)
    divergences_before <- sampler$numDivergences
    signal_sampler$adaptive <- FALSE
    compiled$mcmc$run(samples,reset=FALSE,resetMV=TRUE)
    diagnostics[[ch]] <<- data.frame(chain=ch,
      divergences=sampler$numDivergences-divergences_before)
    as.matrix(compiled$mcmc$mvSamples)
  })
  draws <- posterior::as_draws_array(aperm(simplify2array(chains),c(1,3,2)))
  list(draws=draws,seconds=proc.time()[[3]]-started,compile_seconds=compile_seconds,
       likelihood_error=abs(actual-expected),diagnostics=do.call(rbind,diagnostics))
}
