# Alternating NUTS and Metropolis

From the repository root:

```sh
Rscript -e 'install.packages(c("nimbleHMC", "posterior"), repos="https://cloud.r-project.org")'
Rscript examples/nimble_joint.R
```

NIMBLE compiles the model and samplers with a C++ toolchain; Stan is not needed
for this example. `R/nimble_joint.R` defines the model and assigns two built-in
samplers. Each iteration updates the spline coefficients and two smoothing
precisions with NUTS, then `(A, n*f0, n^2*fdot)` with block Metropolis. Both blocks
use the current values of the other block. Adaptation ends after warmup.

The inputs are the same n=1024 LS2 realization and SNR-40, known-phase cosine
chirp as the [R blocked example](BLOCKED.md). The priors and moving-Whittle
likelihood also match. The real and imaginary Fourier observations have means
given by the proposed waveform transform and variance S/2 each. The model checks
its likelihood and normalized priors against R before sampling, including after
compilation and a change to both parameter blocks.

`R/joint_data.R` prepares the inputs; `examples/nimble_joint.R` runs four chains
and saves `case.rds`, `summary.csv`, `session.txt`, and signal/PSD comparison and
diagnostic PNGs under `results/nimble-ls2-TIMESTAMP/`. Edit `warmup` and `samples`
in that driver to change the budget. Compilation and sampling times are separate.

The saved gate requires all monitored parameters to have R-hat <1.01 and
bulk/tail ESS >100, with no post-warmup NUTS divergences. This is a convergence
screen, not interval calibration or a global search for signal modes.
NIMBLE's diagnostic interface differs from Stan's; this example does not report
Stan-style E-BFMI or tree-depth diagnostics. Console divergence totals include
warmup; the saved counts exclude it. Signal acceptance rates are also saved.

References: [NIMBLE sampler configuration](https://r-nimble.org/manual/cha-mcmc.html),
[nimbleHMC](https://cran.r-project.org/package=nimbleHMC).
