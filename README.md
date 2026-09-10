# Moving-periodogram P-spline demo in R

Estimate a time-varying PSD from simulated LS2 data using a moving periodogram,
tensor P-splines and Stan NUTS. A repeated study measures error, pointwise coverage,
interval width and runtime for n=512–4096.

## Run

From the repository root:

```sh
Rscript setup.R             # install missing CmdStanR, posterior and CmdStan
Rscript examples/one_fit.R  # one PSD fit
Rscript demo.R study        # repeated study
Rscript demo.R --help       # modes and options
Rscript tests.R             # numerical checks
```

Stan requires a C++ toolchain; see the [CmdStanR setup guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Example results

Eight realizations per n; bands show 90% confidence intervals for mean metrics.
[Saved metrics and asset provenance](assets/README.md).

<img src="assets/scaling.png" width="450" alt="Repeated-study error, coverage, interval width and runtime">

![True and posterior median spectrum](assets/spectra.png)

<details>
<summary>Sampler diagnostics</summary>

![NUTS traces and autocorrelations](assets/convergence.png)

</details>

## Signal + PSD example

Estimate a chirp's amplitude, `(f0, fdot)` and PSD in 2048 LS2 samples at injected SNR 40,
alternating NUTS updates for the PSD with block Metropolis updates for the signal:

```sh
Rscript -e 'install.packages(c("nimbleHMC", "posterior"), repos="https://cloud.r-project.org")'
Rscript examples/nimble_joint.R
```

[NIMBLE sampler details](examples/NIMBLE.md). Uses built-in samplers; requires a
C++ toolchain. The earlier [R Gibbs/Metropolis example](examples/BLOCKED.md) is also available.

![Joint signal and noise inference](assets/nimble_joint.png)

<details>
<summary>PSD uncertainty and sampler diagnostics</summary>

![PSD slices with pointwise intervals](assets/nimble_psd_slices.png)

![Signal traces and posterior intervals](assets/nimble_signal_diagnostics.png)

</details>

## Load results

Each run prints its output directory under `results/`. Studies save figures,
CSV metrics and summaries, posterior draws, settings and software versions.
Replace the path below with your study directory:

```r
out <- "results/study-YYYYMMDD-HHMMSS"
metrics <- read.csv(file.path(out, "metrics.csv"))
case <- readRDS(file.path(out, "representative.rds"))
head(metrics)
case$diagnostics
head(case$surface)

# Redraw without rerunning MCMC:
source("R/plot.R")
plot_spectra(case, file.path(out, "spectra.png"))
```

## Where the functions live

| File | Main functions / purpose |
|---|---|
| `R/simulate.R` | `simulate_ls2`, `true_psd_ls2`: simulation and analytic PSD |
| `R/moving_periodogram.R` | `moving_periodogram`: scattered time-frequency ordinates |
| `R/pspline.R` | `spline_setup`, `basis`: spline bases and penalties; `prepare_model`, `fit_pspline`: Stan fitting; `summarize_surface`, `fit_diagnostics`: summaries |
| `R/plot.R` | `study_summary`, `plot_study`, `plot_convergence`, `plot_spectra` |
| `model.stan` | PSD likelihood and spline prior |
| `R/chirp.R`, `R/joint_data.R` | Chirp waveform, Fourier transform, SNR and shared LS2/chirp inputs |
| `R/nimble_joint.R`, `examples/nimble_joint.R` | `fit_nimble_joint`: alternating NUTS/Metropolis fit and example driver |
| `R/joint_summary.R` | `summarize_joint`: surfaces, diagnostics and fit-quality metrics |
| `R/blocked.R`, `R/blocked_plot.R` | `blocked_chain`, `plot_blocked`, `plot_joint_slices`: R sampler and shared figures |
| `examples/one_fit.R`, `examples/blocked_joint.R` | Single-fit and signal/PSD examples |
| `joint.stan`, `examples/stan_joint.R` | Experimental all-parameter NUTS example |
| `demo.R`, `examples/nimble_study.R` | PSD-only and signal/PSD repeated studies |
| `setup.R`, `tests.R` | Dependencies and numerical checks |

[Validation notes](VALIDATION.md) ·
[Dynamic Whittle reference](https://doi.org/10.1080/01621459.2025.2594191)
