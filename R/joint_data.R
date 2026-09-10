# Shared LS2/chirp inputs. Source simulate, moving_periodogram, pspline and chirp first.
# Returns data, spline design, exact Fourier operator and a data-only initialization.
prepare_joint <- function(n=1024L, seed=4821L) {
  set.seed(seed)
  m <- 16L; target_snr <- 40
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
  list(n=n,m=m,truth=truth,signal=signal,noise=noise,x=x,ordinates=ordinates,
       data=data,transform=transform,template=template,prep=prep,B=B,lt=lt,lf=lf,
       initial=initial,bounds=bounds,bounds_hz=bounds_hz,theta_units=theta_units,
       target_snr=target_snr,whittle_snr=whittle_snr,seed=seed)
}
