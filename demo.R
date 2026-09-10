# Run from the repository root: Rscript demo.R [smoke|pilot|study] [replicates] [seed-base]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && args[1] %in% c("--help", "-h")) {
  cat("Usage: Rscript demo.R [smoke|pilot|study] [replicates: 1-99] [seed-base]\n",
      "Defaults: pilot; 8 replicates for study, otherwise 1; seed base 300000.\n",
      "Run from the repository root. See README.md for setup and outputs.\n", sep = "")
  quit(status = 0)
}
if (!file.exists("model.stan")) stop("Run from the repository root containing model.stan")
source("R/simulate.R")
source("R/moving_periodogram.R")
source("R/pspline.R")
source("R/plot.R")
mode <- if (length(args)) args[1] else "pilot"
stopifnot(mode %in% c("smoke", "pilot", "study"))
if (!requireNamespace("cmdstanr", quietly = TRUE) || !requireNamespace("posterior", quietly = TRUE))
  stop("Run Rscript setup.R first (CmdStanR and posterior required)")
replicates <- if (length(args) >= 2) as.integer(args[2]) else if (mode == "study") 8L else 1L
seed_base <- if (length(args) >= 3) as.integer(args[3]) else 300000L
stopifnot(is.finite(replicates), replicates >= 1L, replicates < 100L,
          is.finite(seed_base), seed_base >= 0L, seed_base < 1e9)
sizes <- switch(mode, smoke = 512L, pilot = c(512L, 4096L), study = 2^(9:12))
warmup <- 750L
samples <- 750L
out <- file.path("results", paste0(mode, "-", format(Sys.time(), "%Y%m%d-%H%M%S")))
dir.create(out, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
compile_start <- proc.time()[[3]]
model <- cmdstanr::cmdstan_model("model.stan")
compile_seconds <- proc.time()[[3]] - compile_start
# Frozen starting design, not tuned to force nominal coverage. Window grows
# slowly with n; the fixed interior grid lies inside every retained domain.
grid <- expand.grid(u = seq(0.1, 0.85, length.out = 40),
                    omega = pi * seq(0.12, 0.88, length.out = 40))
# Simulation uses zero-based t/n; moving-periodogram coordinates use t (1-based).
# Subtract 1/n when evaluating truth below to align those conventions exactly.
settings <- list(mode = mode, sizes = sizes, replicates = replicates, seed_base = seed_base,
  warmup = warmup, samples = samples, chains = 4L, thin = 2L,
  m_rule = "round(16*(n/1024)^(1/3))",
  kt_rule = "round(8*(n/1024)^(1/4))", kf = 10L,
  grid = grid, compile_seconds = compile_seconds,
  cmdstan_version = as.character(cmdstanr::cmdstan_version()))
saveRDS(settings, file.path(out, "settings.rds"))
results <- list()
for (n in sizes) for (rep in seq_len(replicates)) {
  seed <- as.integer(seed_base + 100 * match(n, 2^(9:15)) + rep)
  set.seed(seed)
  message("Fitting n=", n, ", replicate=", rep)
  # Includes data generation, ordinates, fitting and posterior summaries;
  # excludes one-time compilation, file writes and figure rendering.
  start <- proc.time()[[3]]
  m <- as.integer(round(16 * (n / 1024)^(1/3)))
  ord <- moving_periodogram(simulate_ls2(n), m = m)
  stopifnot(all(range(grid$u) >= min(ord$u)), all(range(grid$u) <= max(ord$u)),
            min(grid$omega) >= min(ord$omega), max(grid$omega) <= max(ord$omega))
  kt <- as.integer(round(8 * (n / 1024)^(1/4)))
  prep <- prepare_model(ord, kt = kt)
  fit <- fit_pspline(prep, model, seed, warmup, samples)
  diagnostics <- fit_diagnostics(fit)
  attempts <- 1L
  if (!diagnostics$ok) {
    message("Extending the same dataset to twice the warmup and retained draws")
    attempts <- 2L
    fit <- fit_pspline(prep, model, seed, 2L * warmup, 2L * samples)
    diagnostics <- fit_diagnostics(fit)
  }
  surface <- summarize_surface(fit, prep, grid)
  truth <- log(true_psd_ls2(grid$u - 1/n, grid$omega))
  row <- cbind(data.frame(n = n, replicate = rep, seed = seed, m = m, kt = kt, attempts = attempts,
    ordinates = nrow(ord), rmse = sqrt(mean((surface$mean_log - truth)^2)),
    coverage = mean(surface$lower_log <= truth & truth <= surface$upper_log),
    width = mean(surface$upper_log - surface$lower_log), seconds = proc.time()[[3]] - start), diagnostics)
  results[[length(results) + 1L]] <- row
  write.csv(do.call(rbind, results), file.path(out, "metrics.csv"), row.names = FALSE)
  saveRDS(list(surface = surface, truth = truth, diagnostics = diagnostics,
               summary = fit$summary(), draws = fit$draws(c("c", "phi"))),
          file.path(out, sprintf("n%d-rep%d.rds", n, rep)))
  if (n == max(sizes) && rep == 1L) {
    case <- list(surface = surface, truth = truth, diagnostics = diagnostics,
                 summary = fit$summary(), draws = fit$draws(c("c", "phi", "lp__")))
    saveRDS(case, file.path(out, "representative.rds"))
    for (ext in c("pdf", "png")) {
      plot_convergence(case, file.path(out, paste0("convergence.", ext)))
      plot_spectra(case, file.path(out, paste0("spectra.", ext)))
    }
  }
  print(row)
}
results <- do.call(rbind, results)
write.csv(study_summary(results), file.path(out, "summary.csv"), row.names = FALSE)
if (all(results$ok)) {
  plot_study(results, file.path(out, "scaling.pdf"))
  plot_study(results, file.path(out, "scaling.png"))
} else
  warning("Some fits failed diagnostics. Results saved, scaling plot withheld.")
writeLines(c("Points: independent-replicate means; ribbons: approximate 90% t intervals for those means.",
  "Coverage: average pointwise inclusion across the common interior grid, not simultaneous coverage.",
  "Uniform knots: time count grows as n^(1/4), frequency count=10; m grows as n^(1/3), thinning=2.",
  "Time includes generation, transformation, sampling and summaries; compilation excluded.",
  paste("Mode:", mode, "; replicates per size:", replicates),
  "Intervals use independent simulations, not correlated time-frequency grid cells, as replicates.",
  "Eight replicates provide an initial Monte Carlo check, not a general calibration result."),
  file.path(out, "caption.txt"))
message("Saved results to ", out)
