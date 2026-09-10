# Same LS2/chirp data and priors as blocked_joint.R, with joint Stan NUTS.
# Run from the repository root after setup.R: Rscript examples/stan_joint.R
source("R/simulate.R"); source("R/moving_periodogram.R")
source("R/pspline.R"); source("R/chirp.R"); source("R/blocked.R")
if (!requireNamespace("cmdstanr",quietly=TRUE)) stop("Run Rscript setup.R first")
source("R/joint_data.R")
with(prepare_joint(), {
# The indexing and weights are the exact same cached transform used by R.
L <- 2*m+1L
centers <- round(ordinates$u*n)
index <- outer(centers-m-1L,seq_len(L),"+")
weights <- exp(-1i*outer(ordinates$omega,0:(L-1))) / sqrt(2*pi*L)
stan_data <- list(T=n,N=length(data),K=ncol(B),L=L,t=(0:(n-1))/n,
  index=as.integer(index),wr=Re(weights),wi=Im(weights),dre=Re(data),dimag=Im(data),
  B=B,lt=lt,lf=lf,f0_lower=bounds[1,1],f0_upper=bounds[1,2],
  fdot_lower=bounds[2,1],fdot_upper=bounds[2,2])
g <- template(initial)
A_initial <- sum(Re(Conj(g)*data))/sum(Mod(g)^2)
q <- blocked_precision(log(c(2,2)),lt,lf)
c_initial <- as.vector(solve(crossprod(B)+diag(q),
  crossprod(B,log(pmax(Mod(data-A_initial*g)^2,1e-12))+.5772156649)))
model <- cmdstanr::cmdstan_model("joint.stan")
warmup <- 1000L; samples <- 1000L
started <- proc.time()[[3]]
fit <- model$sample(data=stan_data,seed=501,chains=4,parallel_chains=4,
  iter_warmup=warmup,iter_sampling=samples,metric="dense_e",adapt_delta=.95,
  max_treedepth=12,step_size=.01,refresh=0,
  init=function() list(A=A_initial,cycles_f0=initial[1]+rnorm(1,0,.001),
    cycles_fdot=initial[2]+rnorm(1,0,.001),c=c_initial+rnorm(ncol(B),0,.01),phi=c(2,2)))
draws <- fit$draws(c("A","f0","fdot","phi","c","lp__"))
summary <- posterior::summarise_draws(draws)
diagnostics <- fit_diagnostics(fit,c("A","cycles_f0","cycles_fdot","c","phi"))
ok <- diagnostics$ok
print(summary[summary$variable %in% c("A","f0","fdot","phi[1]","phi[2]"),])
print(diagnostics)
message("Elapsed: ",round(proc.time()[[3]]-started,1)," s; joint NUTS gate: ",ok)
# Regular display grid: evaluate both waveform transforms at identical centers.
grid <- expand.grid(u=round(seq(.1,.85,length.out=40)*n)/n,
                    f=seq(.06,.44,length.out=40))
display_transform <- fourier_operator(n,m,round(grid$u*n),2*pi*grid$f)
bt <- basis(prep$time,grid$u) %*% prep$time$vectors
bf <- basis(prep$frequency,2*grid$f) %*% prep$frequency$vectors
B_display <- do.call(cbind,lapply(seq_len(ncol(bf)),function(j) bt*bf[,j]))
out <- file.path("results",paste0("stan-ls2-",format(Sys.time(),"%Y%m%d-%H%M%S")))
dir.create(out,recursive=TRUE)
case <- list(grid=grid,truth=truth,S=true_psd_ls2(grid$u-1/n,2*pi*grid$f),
             h=display_transform(signal),data=data,B=B_display,draws=draws,summary=summary,
             diagnostics=diagnostics,ok=ok,signal=signal,noise=noise,x=x,ordinates=ordinates,
             settings=list(n=n,m=m,thin=2,kt=8,kf=6,target_snr=target_snr,
               whittle_snr=whittle_snr,simulation_seed=4821,sampler="Stan NUTS",seed=501,
               warmup=warmup,samples=samples,initial=initial/theta_units,bounds=bounds_hz),
             seconds=proc.time()[[3]]-started)
saveRDS(case,file.path(out,"case.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
write.csv(as.data.frame(summary),file.path(out,"summary.csv"),row.names=FALSE)
saveRDS(stan_data,file.path(out,"stan_data.rds"))
saveRDS(fit$draws("log_lik"),file.path(out,"log_lik.rds"))
source("R/blocked_plot.R")
plot_blocked(case,out)
if(!ok) warning("Some chain diagnostics failed. Inspect the saved results before interpreting recovery.")
message("Saved joint Stan LS2 example to ",out)
})
