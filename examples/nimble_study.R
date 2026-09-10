# Small repeated-data check; six fits, not a coverage-calibration campaign.
source("R/simulate.R"); source("R/moving_periodogram.R")
source("R/pspline.R"); source("R/chirp.R"); source("R/blocked_plot.R")
source("R/joint_data.R"); source("R/nimble_joint.R"); source("R/joint_summary.R")
sizes <- c(1024L,2048L); seeds <- 4821:4823
warmup <- 1000L; samples <- 3000L
out <- file.path("results",paste0("nimble-study-",format(Sys.time(),"%Y%m%d-%H%M%S")))
dir.create(out,recursive=TRUE)
metrics <- list()
for(n in sizes) for(seed in seeds) {
  d <- prepare_joint(n,seed)
  fit <- fit_nimble_joint(d,warmup,samples)
  case <- summarize_joint(d,fit,warmup,samples)
  saveRDS(case,file.path(out,sprintf("n%d-seed%d.rds",n,seed)))
  metrics[[length(metrics)+1L]] <- case$metrics
  write.csv(do.call(rbind,metrics),file.path(out,"metrics.csv"),row.names=FALSE)
  print(case$metrics)
  # Select the representative case in advance, including if diagnostics fail.
  if(n==max(sizes) && seed==seeds[1]) {
    plot_blocked(case,out); plot_joint_slices(case,out)
    saveRDS(case,file.path(out,"representative.rds"))
  }
}
results_table <- do.call(rbind,metrics)
study_means <- aggregate(results_table[c("rmse_log","coverage","width_log","seconds")],
  by=list(n=results_table$n),FUN=mean)
study_means$replicates <- length(seeds)
write.csv(study_means,file.path(out,"study_summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"session.txt"))
if(!all(do.call(rbind,metrics)$ok)) warning("Some fits failed diagnostics; all results are retained.")
message("Saved repeated-data check to ",out)
