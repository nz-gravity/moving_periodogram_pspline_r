# Moving-periodogram P-spline demo in R

A small, self-contained LS2 simulation demonstration: generate a time series,
compute a frequency-cycling moving periodogram, fit a smooth Bayesian log-spectrum,
and measure error, pointwise interval coverage, interval width and runtime as n grows.
R handles the numerical work and plots; a 33-line Stan model runs NUTS via CmdStanR.
There is no Python dependency or Python-parity requirement.

## Start here: installation and first run

This is a **small collection of R scripts, not an installable R package**.
Do not use `install.packages("moving_periodogram_pspline_r")` or `install_github()`.
No Python, external dataset or GPU is needed. RStudio is optional.

### 1. Install prerequisites once

Install [R from CRAN](https://cran.r-project.org/) and a C++ build toolchain:

| System | Toolchain |
|---|---|
| macOS | Xcode Command Line Tools: run `xcode-select --install` in Terminal if missing. |
| Windows | Install [Rtools matching your R version](https://cran.r-project.org/bin/windows/Rtools/); follow CmdStan's Windows setup instructions. |
| Linux | A C++ compiler and GNU Make; on Ubuntu/Debian, `sudo apt install g++ make`. |

See the official [CmdStan installation guide](https://mc-stan.org/docs/cmdstan-guide/installation.html)
for platform details. A simple writable project path without spaces is safest for
compilation, especially on Windows. Installing CmdStan can take several minutes
and needs internet access. Later simulations run locally without downloading data.

Tested on macOS with R 4.5.3, CmdStanR 0.9.0, CmdStan 2.39.0 and posterior 1.7.0.
Windows/Linux instructions follow upstream guidance but have not been tested here.
Dependencies are not version-locked; each study saves its actual versions.

### 2. Get the code and set the working directory

Clone the repository, or download and extract its ZIP from GitHub:

```sh
git clone https://github.com/nz-gravity/moving_periodogram_pspline_r.git
cd moving_periodogram_pspline_r
```

**All commands below run from the folder containing `README.md`, `demo.R` and
`model.stan`.** The files in `R/` are loaded with `source()`.

### 3. Install the R dependencies and run a small example

In a terminal (these are shell commands, not R expressions):

```sh
Rscript setup.R             # installs missing cmdstanr, posterior and CmdStan
Rscript tests.R             # expected: All R numerical checks passed.
Rscript examples/one_fit.R  # a commented walkthrough of one 1024-point fit
```

`setup.R` uses your R library (a personal library if the default is not writable),
checks the compiler, and installs CmdStan in its default user directory if absent.
An existing configured CmdStan installation is reused. See
[CmdStanR's getting-started guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
for installation details. Setup reports the resolved CmdStan path.

**RStudio alternative:** choose *Session → Set Working Directory → Choose Directory*
and select this repository. Then run in the R console:

```r
source("setup.R")
source("tests.R")
source("examples/one_fit.R")
```

Read `examples/one_fit.R` from top to bottom: simulate → moving periodogram →
splines → NUTS → diagnostics → posterior surface → figures. It prints convergence
checks and the output directory. Open `convergence.png` and `spectra.png` there.
If a diagnostic fails, the example stops rather than presenting the fit as converged.

### 4. Run the repeated study

In a terminal:

```sh
Rscript demo.R --help
Rscript demo.R smoke               # one n=512 fit
Rscript demo.R pilot               # n=512 and 4096, one realization each
Rscript demo.R study 8 300000       # four sizes, eight realizations per size
```

| Argument | Meaning | Default |
|---|---|---|
| Mode | `smoke`, `pilot` or `study` | `pilot` |
| Replicates | Independent datasets at each n, integer 1–99 | 8 for study, otherwise 1 |
| Seed base | Integer 0–999999999; change it for a fresh set of datasets | 300000 |

Using the same seed base repeats the same simulated datasets; changing only the
replicate count retains the earlier seeds. Numerical posterior draws can vary
across software versions/platforms. Each invocation starts a **new** output folder;
there is no automatic resume or reuse of old fits.

The current study uses n=512,1024,2048,4096. Four chains run in parallel, each with
750 warmup iterations and 750 retained draws. After compilation, the eight-repeat
study took about four minutes here (roughly 5–11 seconds per dataset); other machines
will differ. `smoke` and `pilot` are execution checks, not calibration studies.

From the RStudio console, launch the same script using R's own executable:

```r
system2(file.path(R.home("bin"), "Rscript"),
        c("demo.R", "study", "8", "300000"))
```

## Find and inspect the outputs

Each study prints its folder, for example `results/study-YYYYMMDD-HHMMSS/`:

| Output | Contents |
|---|---|
| `scaling.png`, `scaling.pdf` | Error, coverage, width and runtime, with uncertainty across replicates |
| `convergence.png`, `convergence.pdf` | Traces and autocorrelations for the preselected representative fit |
| `spectra.png`, `spectra.pdf` | True and posterior median spectra, sharing a color scale |
| `metrics.csv` | One row per dataset, including diagnostics and number of sampling attempts |
| `summary.csv` | Mean, standard error and 90% confidence interval for each metric and n |
| `n*-rep*.rds` | Posterior surface, log truth, diagnostics, parameter summaries and draws for each fit |
| `representative.rds` | First replicate at largest n; includes log-posterior draws for diagnostics |
| `settings.rds`, `session.txt`, `caption.txt` | Configuration, software versions and figure interpretation |

The teaching example saves only `representative.rds` and the two PNG figures.
Generated results and compiled executables are Git-ignored: a fresh clone will
not include previous figures. Run the scripts to regenerate them.

To inspect an existing run in R, replace `out` with the folder printed by your run:

```r
out <- "results/study-YYYYMMDD-HHMMSS"
metrics <- read.csv(file.path(out, "metrics.csv"))
case <- readRDS(file.path(out, "representative.rds"))
head(metrics)
case$diagnostics
head(case$surface)   # mean_log, lower_log, upper_log, median_psd

# Redraw a figure without rerunning MCMC:
source("R/plot.R")
plot_spectra(case, file.path(out, "spectra.png"))
```

The scaling figure requires every fit to pass the diagnostic gate. A failed check
triggers one repeat on the same dataset with twice the warmup and retained draws.
`attempts` records this; runtime includes both attempts. If that still fails,
results remain available but the scaling plot is withheld. No datasets are dropped.

## Model and conventions

- LS2 is a locally stationary MA(1): `x[t] = w[t] + a(u)*w[t-1]`, where
  `a(u) = 1.1*cos(1.5-cos(4*pi*u))`. Innovations are independent standard normals.
- Spectra use density per angular frequency: white noise has `S = 1/(2*pi)`.
  The analytic LS2 spectrum is `(1+a^2+2*a*cos(omega))/(2*pi)`.
- Each moving window has length `2*m+1`. Its frequency cycles through the m positive
  Fourier frequencies as the center advances. Keep every second block of m centers.
  No padding: only full windows and full retained blocks; a trailing remainder is dropped.
- The dynamic Whittle log-likelihood is `sum(-eta - I*exp(-eta))`, `eta = log(S)`.
  Its independence assumption is approximate because the windows overlap.
- Cubic B-splines use uniform interior knots on normalized time and frequency.
  Penalties integrate products of second derivatives (exact spanwise quadrature),
  each normalized by its trace. The tensor prior has precision eigenvalues
  `phi_time*lambda_time + phi_freq*lambda_freq + 1e-6` outside the joint null space.
  The four joint-null directions have fixed precision `1e-4`. Both smoothing
  precisions have Gamma(shape=2, rate=1) priors and are sampled jointly with the
  centered eigen-coefficients. The normal density includes its normalizing determinant.
- Frequency-rung reuse avoids a large dense observation-by-tensor-basis matrix.

The design is intentionally a starting configuration: n=512,1024,2048,4096;
`m=round(16*(n/1024)^(1/3))`; `round(8*(n/1024)^(1/4))` time and 10 frequency interior knots.
The growing time basis follows an initial pilot in which fixed 8 time knots
undercovered at n=32768. This is a development choice, not independent calibration.
Fixed spline resolution or local-window bias can produce an error plateau or
undercoverage. No configuration is tuned to force the desired plot shape.

## Reading the figure

All runs use a 40x40 evaluation grid with time in [0.10,0.85] and angular frequency
in [0.12*pi,0.88*pi], checked to lie inside every retained observation domain.
The simulator's zero-based time convention is accounted for when evaluating truth.

1. RMSE of posterior mean **log** spectrum relative to analytic log spectrum, in nats.
2. Fraction of grid locations whose pointwise 5th–95th percentile interval contains
   the analytic spectrum. The dotted line is nominal 90%; this is not simultaneous coverage.
3. Average 90% interval width in log-spectrum units (nats).
4. Elapsed generation, transformation, sampling and summary time. Compilation and
   output/plot writes are excluded and compilation time is recorded separately.

Points are means over independent realizations. Ribbons are approximate 90%
Student-t confidence intervals for those means: mean +/- t(0.95, R-1)*sd/sqrt(R),
with coverage bounds clipped to [0,1]. These describe Monte Carlo uncertainty
in the mean metric, not posterior credible bands or the spread of individual
realizations. Eight replicates are an initial check; the t approximation can be
rough with so few realizations. `summary.csv` records the mean, standard error
and interval for each metric. A single replicate has no ribbon. The requested four aligned panels share the n axis;
error, width and runtime use logarithmic y axes, coverage uses a linear y axis.
The PDF uses one blue series and a dotted black nominal reference.
Pilot/smoke results demonstrate execution only; repeated simulations are needed
before making a coverage claim. Neighboring grid points are not independent replicates.

Four chains must have R-hat <1.01, minimum bulk/tail ESS >=100, no divergences or
maximum-tree-depth hits, and minimum E-BFMI >=0.3 to pass the plotting gate.
These checks assess computation, not statistical calibration.

The convergence figure shows post-warmup traces and autocorrelations for both
smoothing parameters, the log posterior and the coefficient with the largest
R-hat. Four chain colors and line styles are consistent across panels. The
footer reports diagnostics across all coefficients and smoothing parameters.
The case is selected in advance (first replicate at the largest n), not by fit quality.
These plots help assess mixing and do not prove convergence.

The spectrum figure compares analytic S(t,f) and the pointwise posterior median
S(t,f) on the common interior grid, using identical linear color limits in density
per radian. The posterior median is calculated from the retained surface draws;
it is distinct from the posterior-mean log-spectrum used for RMSE.

## Where the functions live

Start with `examples/one_fit.R` for the workflow, then inspect these functions.
Definitions have short input/output comments; this is source code, not generated
package help (`?fit_pspline` is not available).

| File | Functions | Input → output |
|---|---|---|
| `R/simulate.R` | `simulate_ls2`, `true_psd_ls2` | Length → time series; paired time/frequency coordinates → analytic density |
| `R/moving_periodogram.R` | `moving_periodogram` | Time series + window/thinning → dataframe of scattered ordinates |
| `R/pspline.R` | `spline_setup`, `basis` | Interior-knot count → penalty/eigenbasis; coordinates → spline matrix |
| `R/pspline.R` | `prepare_model`, `fit_pspline` | Ordinates → Stan inputs → CmdStanMCMC fit |
| `R/pspline.R` | `summarize_surface`, `fit_diagnostics` | Fit → surface summaries or one-row diagnostic report |
| `R/plot.R` | `study_summary`, `plot_study` | Replicate metrics → summary table or four-panel figure |
| `R/plot.R` | `plot_convergence`, `plot_spectra` | Saved representative case → diagnostic/comparison figures |
| `model.stan` | Stan data, parameters and model blocks | Dynamic Whittle likelihood and Gaussian tensor prior |
| `demo.R` | Experiment driver | Runs all datasets, checks diagnostics, saves results and plots |
| `setup.R`, `tests.R` | Installation and numerical checks | Setup dependencies; verify core numerical conventions |

### Small changes to try

- **More replicates/new datasets:** use the command-line arguments above.
- **Sampler budget:** edit `warmup` and `samples` near the top of `demo.R`.
- **Time-series lengths:** edit `sizes` in `demo.R`; the current seed mapping supports
  powers of two from 512 to 32768. Other sizes also require adjusting that mapping.
- **Knots/window:** edit the `kt` and `m` expressions inside the loop; keep the
  matching descriptions in `settings` consistent. `kf` defaults to 10 in `prepare_model`.
- **Understand one component:** source its `R/` file and call it directly. For example,
  `source("R/simulate.R"); x <- simulate_ls2(1024)` does not start MCMC.

Keep four chains for convergence checks. `fit_pspline` controls `parallel_chains`;
reduce parallelism there on a smaller machine rather than removing diagnostic chains.
Use a fixed configuration across replicates. The functions accept ordinary R vectors,
matrices and dataframes; no package classes are needed until CmdStan returns a fit.

## Troubleshooting

| Symptom | Next step |
|---|---|
| `Rscript` not found | Install R/add it to PATH, or use the RStudio console route above. |
| Cannot open `R/...` or `model.stan` | Check `getwd()`; change to the repository root. |
| Missing `cmdstanr` or `posterior` | Run `setup.R` with the same R installation that runs the demo. |
| Compiler/Make missing | Install the OS toolchain above, restart the terminal/RStudio, rerun setup. |
| CmdStan path not set | Run setup; for an existing installation, use `cmdstanr::set_cmdstan_path(...)` in R or set `CMDSTAN` for shell runs. |
| Sampling completes but no scaling plot | Read `metrics.csv`, especially `ok`, `rhat`, `min_ess`, `divergences` and `attempts`. Investigate or increase the budget. |
| Occasional rejected proposal during warmup | Inspect final diagnostics; a rejected warmup proposal alone is not a failed fit. Persistent warnings require investigation. |

When asking for help, include the command, complete error, `sessionInfo()` and
`cmdstanr::cmdstan_version()`. [VALIDATION.md](VALIDATION.md) records the tests and
pilot results obtained here; it is not a guarantee of calibration on other data.

Method references: [Tang et al., dynamic Whittle](https://doi.org/10.1080/01621459.2025.2594191),
[CmdStanR](https://mc-stan.org/cmdstanr/). This implementation uses R's built-in
`splines::splineDesign`; no code was copied from psplinePsd.
