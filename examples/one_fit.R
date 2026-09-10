# Run from the repository root: Rscript examples/one_fit.R
# In RStudio, set the working directory to the repository, then source this file.
if (!file.exists("model.stan")) stop("Set the working directory to the repository root first")
source("R/simulate.R")
source("R/moving_periodogram.R")
source("R/pspline.R")
source("R/plot.R")

# 1. Generate a reproducible time series. n is the number of observations.
set.seed(300201)
n <- 1024L
x <- simulate_ls2(n)

# 2. Turn it into scattered time-frequency power measurements.
# m is the half-window length; thin skips whole frequency-cycling blocks.
ordinates <- moving_periodogram(x, m = 16L, thin = 2L)
print(head(ordinates))

# 3. Construct the uniform cubic spline bases and roughness penalties.
# kt and kf count INTERIOR knots (the basis has kt+4 and kf+4 columns).
prepared <- prepare_model(ordinates, kt = 8L, kf = 10L)
model <- cmdstanr::cmdstan_model("model.stan") # compiles once, then caches

# 4. Draw posterior samples and check convergence BEFORE interpreting the fit.
fit <- fit_pspline(prepared, model, seed = 300201, warmup = 750L, samples = 750L)
diagnostics <- fit_diagnostics(fit)
print(diagnostics)
if (!diagnostics$ok) stop("Diagnostics failed: increase sampling or investigate before continuing")

# 5. Reconstruct the posterior on an interior grid; no extrapolation is required.
grid <- expand.grid(u = seq(0.1, 0.85, length.out = 40),
                    omega = pi * seq(0.12, 0.88, length.out = 40))
surface <- summarize_surface(fit, prepared, grid)
truth <- log(true_psd_ls2(grid$u - 1/n, grid$omega))

# 6. Save the fit and make the two diagnostic figures. New folder each time.
case <- list(surface = surface, truth = truth, diagnostics = diagnostics,
             summary = fit$summary(), draws = fit$draws(c("c", "phi", "lp__")))
out <- file.path("results", paste0("one-fit-", format(Sys.time(), "%Y%m%d-%H%M%S")))
dir.create(out, recursive = TRUE)
saveRDS(case, file.path(out, "representative.rds"))
plot_convergence(case, file.path(out, "convergence.png"))
plot_spectra(case, file.path(out, "spectra.png"))
message("Saved example to ", out)
