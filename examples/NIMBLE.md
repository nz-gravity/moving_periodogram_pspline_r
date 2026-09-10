# Alternating NUTS and Metropolis

From the repository root:

```sh
Rscript -e 'install.packages(c("nimbleHMC", "posterior"), repos="https://cloud.r-project.org")'
Rscript examples/nimble_joint.R
```

NIMBLE compiles the model and samplers with a C++ toolchain; Stan is not needed
for this example. `R/nimble_joint.R` defines the model and assigns two built-in
samplers. Each iteration updates standardized spline coefficients and two smoothing
precisions with NUTS, then `(A, n*f0, n^2*fdot)` with block Metropolis. Both blocks
use the current values of the other block. Adaptation ends after warmup.

The default uses n=2048 LS2 samples and an SNR-40, known-phase cosine chirp.
Use `Rscript examples/nimble_joint.R 1024 4821` to reproduce the original input
from the [R blocked example](BLOCKED.md). The priors and moving-Whittle likelihood
match; vectorized evaluation and standardized coefficients improve NUTS efficiency. The real and imaginary Fourier observations have means
given by the proposed waveform transform and variance S/2 each. The model checks
its likelihood and normalized priors against R before sampling, including after
compilation and a change to both parameter blocks.

`R/joint_data.R` prepares the inputs; `examples/nimble_joint.R` runs four chains
and saves `case.rds`, `metrics.csv`, `summary.csv`, `session.txt`, signal/PSD
comparison plots, diagnostic plots and `psd_slices.png` under `results/nimble-ls2-TIMESTAMP/`. Edit `warmup` and `samples`
in that driver to change the budget. Compilation and sampling times are separate; the diagnostic table also records
time in each sampler block. `R/joint_summary.R` computes the saved surface and
fit-quality metrics from all retained draws.

The saved gate requires all monitored parameters to have R-hat <1.01 and
bulk/tail ESS >100, with no post-warmup NUTS divergences. This is a convergence
screen, not interval calibration or a global search for signal modes.
NIMBLE's diagnostic interface differs from Stan's; this example does not report
Stan-style E-BFMI or tree-depth diagnostics. Console divergence totals include
warmup; the saved counts exclude it. Signal acceptance rates are also saved.

For a small repeated-data check, run `Rscript examples/nimble_study.R`. It fits
three seeds at n=1024 and n=2048 and saves every result, including any failed
convergence screens. Three replicates assess sensitivity to the realization;
they do not establish calibrated coverage. The displayed case is always the
first seed at n=2048, selected before fitting.

References: [NIMBLE sampler configuration](https://r-nimble.org/manual/cha-mcmc.html),
[nimbleHMC](https://cran.r-project.org/package=nimbleHMC).
