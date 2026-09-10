# Alternating built-in NUTS (PSD) and block Metropolis (signal).
# Rscript examples/nimble_joint.R [n=2048] [simulation_seed=4821]
args <- commandArgs(trailingOnly=TRUE)
if("--help" %in% args) {
  cat("Rscript examples/nimble_joint.R [n=2048] [simulation_seed=4821]\n")
  quit(status=0)
}
n <- if(length(args)>0) as.integer(args[1]) else 2048L
seed <- if(length(args)>1) as.integer(args[2]) else 4821L
if(length(args)>2 || !is.finite(n) || n<512 || n>2048 || !is.finite(seed) || seed<0)
  stop("Use 512 <= n <= 2048 and a nonnegative integer simulation seed.")
source("R/simulate.R"); source("R/moving_periodogram.R")
source("R/pspline.R"); source("R/chirp.R"); source("R/blocked_plot.R")
source("R/joint_data.R"); source("R/nimble_joint.R"); source("R/joint_summary.R")
d <- prepare_joint(n,seed)
warmup <- 1000L; samples <- 3000L
fit <- fit_nimble_joint(d,warmup,samples)
case <- summarize_joint(d,fit,warmup,samples)
out <- file.path("results",paste0("nimble-ls2-",format(Sys.time(),"%Y%m%d-%H%M%S")))
dir.create(out,recursive=TRUE)
saveRDS(case,file.path(out,"case.rds"))
write.csv(as.data.frame(case$summary),file.path(out,"summary.csv"),row.names=FALSE)
write.csv(case$metrics,file.path(out,"metrics.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
plot_blocked(case,out)
plot_joint_slices(case,out)
print(case$metrics); print(fit$diagnostics)
if(!case$ok) warning("Some diagnostics failed; inspect the saved results before interpreting recovery.")
message("Saved NIMBLE example to ",out)
