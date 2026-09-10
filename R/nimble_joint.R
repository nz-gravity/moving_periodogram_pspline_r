# Vectorized complex-normal density. This changes evaluation, not the likelihood.
# A matching generator supports NIMBLE model initialization.
if (!requireNamespace("nimbleHMC",quietly=TRUE))
  stop('Install dependencies: install.packages(c("nimbleHMC", "posterior"))')
dComplex <- nimble::nimbleFunction(
  run=function(x=double(2), re=double(1), im=double(1), eta=double(1),
               log=integer(0,default=0)) {
    ans <- -sum(eta+((x[,1]-re)^2+(x[,2]-im)^2)*exp(-eta))-length(eta)*log(pi)
    returnType(double(0))
    if(log) return(ans)
    return(exp(ans))
  },buildDerivs=TRUE)

rComplex <- nimble::nimbleFunction(
  run=function(n=integer(0),re=double(1),im=double(1),eta=double(1)) {
    returnType(double(2))
    if(n!=1) stop("rComplex supports one matrix draw")
    ans <- matrix(0,nrow=length(eta),ncol=2)
    ans[,1] <- rnorm(length(eta),re,exp(eta/2)/sqrt(2))
    ans[,2] <- rnorm(length(eta),im,exp(eta/2)/sqrt(2))
    return(ans)
  })

# Built-in NUTS for the PSD, then built-in block Metropolis for the signal.
# Input: prepare_joint() result. Output: draws, sampler diagnostics and timings.
fit_nimble_joint <- function(d, warmup=1000L, samples=3000L, seed=501L, chains=4L) {
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
      z[k] ~ dnorm(0,sd=1)
      c[k] <- z[k]/sqrt(q[k])
    }
    wave[1:T] <- A*cos(2*3.141592653589793*(theta[1]*t[1:T]+0.5*theta[2]*t[1:T]^2))
    eta[1:N] <- B[1:N,1:K] %*% c[1:K]
    for(i in 1:N) {
      mur[i] <- inprod(wr[i,1:L],wave[start[i]:end[i]])
      mui[i] <- inprod(wi[i,1:L],wave[start[i]:end[i]])
    }
    obs[1:N,1:2] ~ dComplex(mur[1:N],mui[1:N],eta[1:N])
  })
  g <- d$template(d$initial)
  A0 <- sum(Re(Conj(g)*d$data))/sum(Mod(g)^2)
  q0 <- ifelse(d$lt+d$lf==0,.01,2*d$lt+2*d$lf+1e-6)
  c0 <- as.vector(solve(crossprod(d$B)+diag(q0),
    crossprod(d$B,log(pmax(Mod(d$data-A0*g)^2,1e-12))+.5772156649)))
  inits <- list(A=A0,theta=d$initial,z=c0*sqrt(q0),phi=c(2,2))
  started <- proc.time()[[3]]
  model <- nimble::nimbleModel(code,constants=constants,
    data=list(obs=cbind(Re(d$data),Im(d$data))),inits=inits,buildDerivs=TRUE)
  # Verify both data and normalized priors independently, before and after C++ compilation.
  check_target <- function(m) {
    actual <- m$calculate()
    q <- ifelse(d$lt+d$lf==0,.01,m$phi[1]*d$lt+m$phi[2]*d$lf+1e-6)
    eta <- as.vector(d$B %*% m$c)
    expected <- -sum(eta+Mod(d$data-m$A*d$template(m$theta))^2*exp(-eta))-N*log(pi)+
      sum(dnorm(m$c,0,1/sqrt(q),log=TRUE))+sum(dgamma(m$phi,2,1,log=TRUE))+
      dnorm(m$A,0,5,log=TRUE)+sum(dunif(m$theta,d$bounds[,1],d$bounds[,2],log=TRUE))
    # Jacobian for c = z/sqrt(q); required when checking the transformed target.
    expected <- expected-.5*sum(log(q))
    error <- abs(actual-expected)
    stopifnot(is.finite(error),error<1e-7)
    error
  }
  target_error <- check_target(model)
  conf <- nimble::configureMCMC(model,nodes=NULL,print=FALSE,monitors=c("A","f0","fdot","c","phi","z"))
  conf$addSampler(target=c("z","phi"),type="NUTS",
    control=list(warmupMode="iterations",warmup=warmup,delta=.9))
  # Scale and covariance are in (amplitude, n*f0, n^2*fdot) coordinates.
  propCov <- diag(c(.05,.03,.08)) %*%
    matrix(c(1,0,0,0,1,-.97,0,-.97,1),3) %*% diag(c(.05,.03,.08))
  conf$addSampler(target=c("A","theta"),type="RW_block",control=list(propCov=propCov))
  stopifnot(length(conf$getUnsampledNodes())==0L)
  mcmc <- nimble::buildMCMC(conf)
  # Expose nested sampler fields for adaptation controls and divergence counts.
  old_interfaces <- nimble::getNimbleOption("buildInterfacesForCompiledNestedNimbleFunctions")
  nimble::nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions=TRUE)
  on.exit(nimble::nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions=old_interfaces))
  compiled <- nimble::compileNimble(model,mcmc)
  target_error <- max(target_error,check_target(compiled$model))
  compiled$model$setInits(list(A=A0+.02,theta=d$initial+c(.01,-.02),
    z=(c0+.001)*sqrt(ifelse(d$lt+d$lf==0,.01,1.5*d$lt+2.5*d$lf+1e-6)),phi=c(1.5,2.5)))
  target_error <- max(target_error,check_target(compiled$model))
  compile_seconds <- proc.time()[[3]]-started
  started <- proc.time()[[3]]
  sampler <- compiled$mcmc$samplerFunctions[[1]]
  signal_sampler <- compiled$mcmc$samplerFunctions[[2]]
  stopifnot(length(sampler$numDivergences)==1L,length(signal_sampler$adaptive)==1L)
  diagnostics <- vector("list",chains)
  chain_draws <- lapply(seq_len(chains),function(ch) {
    message("NIMBLE chain ",ch,"/",chains)
    set.seed(seed+ch-1L)
    phi0 <- exp(rnorm(2,log(2),.1))
    qinit <- ifelse(d$lt+d$lf==0,.01,phi0[1]*d$lt+phi0[2]*d$lf+1e-6)
    compiled$model$setInits(list(A=A0+rnorm(1,0,.03),theta=d$initial+rnorm(2,0,.005),
      z=(c0+rnorm(K,0,.01))*sqrt(qinit),phi=phi0))
    signal_sampler$adaptive <- TRUE
    compiled$mcmc$run(warmup,nburnin=warmup)
    divergences_before <- sampler$numDivergences
    signal_sampler$adaptive <- FALSE
    compiled$mcmc$run(samples,reset=FALSE,resetMV=TRUE,time=TRUE)
    draws_ch <- as.matrix(compiled$mcmc$mvSamples)
    diagnostics[[ch]] <<- data.frame(chain=ch,
      divergences=sampler$numDivergences-divergences_before,
      signal_acceptance=mean(diff(draws_ch[,"A"])!=0),
      psd_seconds=compiled$mcmc$getTimes()[1],signal_seconds=compiled$mcmc$getTimes()[2])
    draws_ch
  })
  draws <- posterior::as_draws_array(aperm(simplify2array(chain_draws),c(1,3,2)))
  list(draws=draws,seconds=proc.time()[[3]]-started,compile_seconds=compile_seconds,
       target_error=target_error,diagnostics=do.call(rbind,diagnostics))
}
