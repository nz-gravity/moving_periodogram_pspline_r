# LS2 time series + a high-SNR linear chirp; Gibbs/blocked MH, no Stan needed.
source("R/simulate.R"); source("R/moving_periodogram.R")
source("R/pspline.R"); source("R/chirp.R"); source("R/blocked.R")
if (!requireNamespace("posterior",quietly=TRUE)) stop('Run install.packages("posterior") in R')
set.seed(4821)
n <- 1024L; m <- 16L; target_snr <- 40
truth <- c(A=NA_real_,f0=0.12,fdot=0.00012)
unit_signal <- chirp_time(n,truth[2:3])
truth[1] <- target_snr/sqrt(sum(ls2_whiten(unit_signal)^2))
signal <- truth[1]*unit_signal
noise <- simulate_ls2(n)
x <- signal+noise
ordinates <- moving_periodogram(x,m=m,thin=2L)
transform <- fourier_operator(n,m,round(ordinates$u*n),ordinates$omega)
data <- transform(x)
# Transform the WHOLE waveform for each proposal, retaining its complex phase.
# Internal parameters are cycles across the observation, for well-scaled proposals.
theta_units <- c(n,n^2)
template <- function(th) transform(chirp_time(n,th/theta_units))
prep <- prepare_model(ordinates,kt=8L,kf=6L)
B <- do.call(cbind,lapply(seq_len(prep$data$Kf),function(j)
  prep$data$Bt * prep$data$Bf[prep$data$rung,j]))
lt <- prep$data$lt; lf <- prep$data$lf
S_observed <- true_psd_ls2(ordinates$u-1/n,ordinates$omega)
whittle_snr <- sqrt(2*sum(Mod(transform(signal))^2/S_observed))
bounds_hz <- rbind(c(.08,.18),c(0,.0002))
bounds <- bounds_hz*theta_units
initial <- chirp_initialize(x,bounds_hz)*theta_units
# Profile amplitude to refine the data-only initial search using the actual transform.
objective <- function(th) {g<-template(th); -sum(Re(Conj(g)*data))^2/sum(Mod(g)^2)}
initial <- optim(initial,objective,method="L-BFGS-B",lower=bounds[,1],upper=bounds[,2])$par
message("Injected exact LS2 optimal SNR=",target_snr,"; moving-Whittle SNR=",round(whittle_snr,2))
message("Data-only initial f0, fdot: ",paste(signif(initial/theta_units,6),collapse=", "))
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
