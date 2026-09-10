# Optional integration check after examples/stan_joint.R.
# Rscript examples/check_stan_joint.R results/stan-ls2-TIMESTAMP
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Supply the results directory from stan_joint.R")
source("R/chirp.R"); source("R/blocked.R")
case <- readRDS(file.path(args[1],"case.rds"))
sd <- readRDS(file.path(args[1],"stan_data.rds"))
d <- posterior::as_draws_matrix(case$draws)
# Evaluate at an actual posterior draw, with full precision output. This checks
# the transform inside Stan independently of the R likelihood implementation.
v <- as.numeric(d[1,]); names(v) <- colnames(d)
c <- v[paste0("c[",seq_len(sd$K),"]")]
init <- list(A=v["A"],cycles_f0=sd$T*v["f0"],
  cycles_fdot=sd$T^2*v["fdot"],c=unname(c),phi=unname(v[c("phi[1]","phi[2]")]))
fit <- cmdstanr::cmdstan_model("joint.stan")$sample(data=sd,chains=1,
  iter_warmup=0,iter_sampling=1,fixed_param=TRUE,init=list(init),seed=812,
  sig_figs=16,refresh=0)
transform <- fourier_operator(sd$T,case$settings$m,
  round(case$ordinates$u*sd$T),case$ordinates$omega)
r_ll <- tf_loglik(case$data-v["A"]*transform(chirp_time(sd$T,v[c("f0","fdot")])),
  as.vector(sd$B %*% c))
stan_ll <- as.numeric(fit$draws("log_lik"))
error <- abs(r_ll-stan_ll)
stopifnot(is.finite(error),error<1e-7)
cat("R versus Stan joint log-likelihood absolute difference:",error,"\n")
