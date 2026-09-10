# LS2 time series + a high-SNR linear chirp; Gibbs/blocked MH, no Stan needed.
source("R/simulate.R"); source("R/moving_periodogram.R")
source("R/pspline.R"); source("R/chirp.R"); source("R/blocked.R")
if (!requireNamespace("posterior",quietly=TRUE)) stop('Run install.packages("posterior") in R')
source("R/joint_data.R")
with(prepare_joint(), {
warmup <- 2000L; samples <- 20000L
started <- proc.time()[[3]]
chains <- lapply(1:4,function(ch) {
  message("Blocked chain ",ch,"/4")
  blocked_chain(data,template,B,lt,lf,initial,bounds,theta_units,
                seed=500+ch,warmup=warmup,samples=samples)
})
draws <- posterior::as_draws_array(simplify2array(lapply(chains,`[[`,"draws")) |> aperm(c(1,3,2)))
summary <- posterior::summarise_draws(draws)
ok <- all(is.finite(summary$rhat)) && max(summary$rhat)<1.01 &&
  min(summary$ess_bulk,summary$ess_tail)>100
print(summary[summary$variable %in% c("A","f0","fdot","phi[1]","phi[2]"),])
print(do.call(rbind,lapply(chains,`[[`,"acceptance")))
message("Elapsed: ",round(proc.time()[[3]]-started,1)," s; all-parameter diagnostic gate: ",ok)
# Regular display grid: evaluate both waveform transforms at identical centers.
grid <- expand.grid(u=round(seq(.1,.85,length.out=40)*n)/n,
                    f=seq(.06,.44,length.out=40))
display_transform <- fourier_operator(n,m,round(grid$u*n),2*pi*grid$f)
bt <- basis(prep$time,grid$u) %*% prep$time$vectors
bf <- basis(prep$frequency,2*grid$f) %*% prep$frequency$vectors
B_display <- do.call(cbind,lapply(seq_len(ncol(bf)),function(j) bt*bf[,j]))
out <- file.path("results",paste0("blocked-ls2-",format(Sys.time(),"%Y%m%d-%H%M%S")))
dir.create(out,recursive=TRUE)
case <- list(grid=grid,truth=truth,S=true_psd_ls2(grid$u-1/n,2*pi*grid$f),
             h=display_transform(signal),data=data,B=B_display,draws=draws,summary=summary,
             chains=chains,ok=ok,signal=signal,noise=noise,x=x,ordinates=ordinates,
             settings=list(n=n,m=m,thin=2,kt=8,kf=6,target_snr=target_snr,
               whittle_snr=whittle_snr,simulation_seed=4821,chain_seeds=501:504,
               warmup=warmup,samples=samples,initial=initial/theta_units,bounds=bounds_hz),
             seconds=proc.time()[[3]]-started)
saveRDS(case,file.path(out,"case.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
write.csv(as.data.frame(summary),file.path(out,"summary.csv"),row.names=FALSE)
source("R/blocked_plot.R")
plot_blocked(case,out)
if(!ok) warning("Some chain diagnostics failed. Inspect the saved results before interpreting recovery.")
message("Saved blocked LS2 example to ",out)
})
