# Common display grid and summaries; uses all draws for interval/accuracy metrics.
# These grid metrics describe one realization, not repeated-data calibration.
summarize_joint <- function(d,fit,warmup,samples,chain_seed=501L) {
  grid <- expand.grid(u=round(seq(.1,.85,length.out=40)*d$n)/d$n,
                      f=seq(.06,.44,length.out=40))
  bt <- basis(d$prep$time,grid$u) %*% d$prep$time$vectors
  bf <- basis(d$prep$frequency,2*grid$f) %*% d$prep$frequency$vectors
  B <- do.call(cbind,lapply(seq_len(ncol(bf)),function(j) bt*bf[,j]))
  summary <- posterior::summarise_draws(fit$draws)
  ok <- all(is.finite(summary$rhat)) && max(summary$rhat)<1.01 &&
    min(summary$ess_bulk,summary$ess_tail)>100 && sum(fit$diagnostics$divergences)==0
  draws <- unclass(posterior::as_draws_matrix(fit$draws))
  eta <- draws[,sprintf("c[%d]",1:ncol(B)),drop=FALSE] %*% t(B)
  surface <- cbind(grid,mean_log=colMeans(eta),
    lower_log=apply(eta,2,quantile,.05),upper_log=apply(eta,2,quantile,.95),
    median_psd=exp(apply(eta,2,median)))
  truth <- true_psd_ls2(grid$u-1/d$n,2*pi*grid$f)
  signal <- colMeans(draws[,c("A","f0","fdot")])
  error <- signal[1]*chirp_time(d$n,signal[2:3])-d$signal
  metrics <- data.frame(n=d$n,seed=d$seed,rmse_log=sqrt(mean((surface$mean_log-log(truth))^2)),
    coverage=mean(surface$lower_log<=log(truth) & surface$upper_log>=log(truth)),
    width_log=mean(surface$upper_log-surface$lower_log),
    residual_signal_snr=sqrt(sum(ls2_whiten(error)^2)),
    rhat=max(summary$rhat),min_ess=min(summary$ess_bulk,summary$ess_tail),
    seconds=fit$seconds,compile_seconds=fit$compile_seconds,ok=ok)
  transform <- fourier_operator(d$n,d$m,round(grid$u*d$n),2*pi*grid$f)
  list(grid=grid,S=truth,h=transform(d$signal),B=B,surface=surface,
    truth=d$truth,signal=d$signal,noise=d$noise,x=d$x,data=d$data,ordinates=d$ordinates,
    draws=fit$draws,summary=summary,diagnostics=fit$diagnostics,metrics=metrics,ok=ok,
    seconds=fit$seconds,compile_seconds=fit$compile_seconds,target_error=fit$target_error,
    settings=list(n=d$n,m=d$m,thin=2,kt=8,kf=6,target_snr=d$target_snr,
      whittle_snr=d$whittle_snr,simulation_seed=d$seed,chain_seeds=chain_seed+seq_len(dim(fit$draws)[2])-1L,
      warmup=warmup,samples=samples,initial=d$initial/c(d$n,d$n^2),bounds=d$bounds_hz,
      sampler="NIMBLE NUTS + RW_block"))
}
