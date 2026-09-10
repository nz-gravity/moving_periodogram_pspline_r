# Alternating built-in NUTS (PSD) and block Metropolis (signal).
# Install once: install.packages(c("nimbleHMC", "posterior"))
source("R/simulate.R"); source("R/moving_periodogram.R")
source("R/pspline.R"); source("R/chirp.R")
source("R/joint_data.R"); source("R/nimble_joint.R")
d <- prepare_joint()
warmup <- 1000L; samples <- 3000L
fit <- fit_nimble_joint(d,warmup,samples)
draws <- fit$draws
summary <- posterior::summarise_draws(draws)
ok <- all(is.finite(summary$rhat)) && max(summary$rhat)<1.01 &&
  min(summary$ess_bulk,summary$ess_tail)>100 && sum(fit$diagnostics$divergences)==0
print(summary[summary$variable %in% c("A","f0","fdot","phi[1]","phi[2]"),])
print(fit$diagnostics)
message("Sampling: ",round(fit$seconds,1)," s; compilation: ",round(fit$compile_seconds,1)," s; gate: ",ok)
with(d, {
# Regular display grid: evaluate both waveform transforms at identical centers.
grid <- expand.grid(u=round(seq(.1,.85,length.out=40)*n)/n,
                    f=seq(.06,.44,length.out=40))
display_transform <- fourier_operator(n,m,round(grid$u*n),2*pi*grid$f)
bt <- basis(prep$time,grid$u) %*% prep$time$vectors
bf <- basis(prep$frequency,2*grid$f) %*% prep$frequency$vectors
B_display <- do.call(cbind,lapply(seq_len(ncol(bf)),function(j) bt*bf[,j]))
out <- file.path("results",paste0("nimble-ls2-",format(Sys.time(),"%Y%m%d-%H%M%S")))
dir.create(out,recursive=TRUE)
case <- list(grid=grid,truth=truth,S=true_psd_ls2(grid$u-1/n,2*pi*grid$f),
             h=display_transform(signal),data=data,B=B_display,draws=draws,summary=summary,
             diagnostics=fit$diagnostics,ok=ok,signal=signal,noise=noise,x=x,ordinates=ordinates,
             settings=list(n=n,m=m,thin=2,kt=8,kf=6,target_snr=40,
               whittle_snr=whittle_snr,simulation_seed=seed,chain_seeds=501:504,
               warmup=warmup,samples=samples,initial=initial/c(n,n^2),bounds=bounds/c(n,n^2),sampler="NIMBLE NUTS + RW_block"),
             seconds=fit$seconds,compile_seconds=fit$compile_seconds,target_error=fit$target_error)
saveRDS(case,file.path(out,"case.rds"))
write.csv(as.data.frame(summary),file.path(out,"summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
source("R/blocked_plot.R")
plot_blocked(case,out)
if(!ok) warning("Some diagnostics failed; inspect the saved results before interpreting recovery.")
message("Saved NIMBLE example to ",out)
})
