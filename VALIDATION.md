# Validation, 10 September 2026

Executed with R 4.5.3, CmdStanR 0.9.0, CmdStan 2.39.0, posterior 1.7.0.
Base-R numerical checks passed (`Rscript tests.R`). All R sources parse;
Stan compiled and executed. No Python comparison was required.

The first 200-warmup/200-draw smoke fit did not meet R-hat <1.01. The default
is now 750 warmup and 750 retained draws per chain, with four parallel chains.
The fixed-eight-time-knot endpoint pilot undercovered at n=32768 despite passing
MCMC diagnostics. The final demonstration grows the uniform time-knot count as
`round(8*(n/1024)^(1/4))`; frequency knots remain fixed at 10.

A subsequent six-size development run (`Rscript demo.R study 1`) passed all
sampling gates at every size: R-hat <1.009, minimum bulk/tail ESS >900, no retained
divergences, no maximum-tree-depth hits, and E-BFMI >0.55. Occasional non-finite
proposals were rejected during warmup; retained-chain diagnostics passed.

| n | Log-RMSE | Pointwise coverage | Mean log-width | Wall seconds |
|---:|---:|---:|---:|---:|
| 1024 | 0.2965 | 0.9281 | 1.1187 | 6.2 |
| 2048 | 0.2574 | 0.9100 | 0.8742 | 7.1 |
| 4096 | 0.2071 | 0.8969 | 0.6779 | 11.4 |
| 8192 | 0.1370 | 0.9475 | 0.5261 | 18.2 |
| 16384 | 0.1241 | 0.9019 | 0.4168 | 35.2 |
| 32768 | 0.0921 | 0.9238 | 0.3308 | 65.0 |

This is one realization per size, reusing endpoint seeds from development.
It does not establish calibration or replicate variability. No uncertainty
ribbons are drawn. A repeated study with new seeds should follow before making
scientific coverage claims. Runtime excludes compilation and file/plot writes.

Local output: `results/study-20260910-135927/` contains metrics, settings,
posterior draws/summaries and the four-panel PNG/PDF. Outputs and the compiled
Stan executable are Git-ignored; source and this validation note remain reviewable.

## Repeated smaller-n study

The follow-up uses n=512,1024,2048,4096, with eight independent realizations
per size and a fresh seed base of 300000. Local output:
`results/study-20260910-141441/`. The default driver now uses this smaller range.

| n | Mean coverage | Approximate 90% CI for mean | Mean wall seconds |
|---:|---:|:---:|---:|
| 512 | 0.9112 | [0.8752, 0.9471] | 4.67 |
| 1024 | 0.9180 | [0.9000, 0.9360] | 5.53 |
| 2048 | 0.9103 | [0.8947, 0.9259] | 9.02 |
| 4096 | 0.9131 | [0.8891, 0.9371] | 11.14 |

Ribbons use Student-t intervals across independent realizations, not correlated
grid cells. They estimate uncertainty of the mean metric, not posterior credible
bands. Eight replicates are preliminary. At n=512, individual-realization coverage
ranged from 0.8313 to 0.9631, illustrating why a one-realization curve fluctuates.
The earlier n=8192 case was not repeated in this smaller study.

All 32 fits pass: maximum R-hat 1.00942, minimum bulk/tail ESS 671, no divergences
or maximum-tree-depth hits. One n=2048 fit initially had R-hat 1.0124; the same
dataset was refit with 1500 warmup and 1500 retained iterations per chain and passed.
No dataset was dropped. `attempts` records the extension and runtime includes both
attempts, explaining the larger timing uncertainty at n=2048. Original metrics are
preserved in `metrics-before-extension.csv`. The driver now automatically allows
one such extension before withholding a scaling plot for persistent failures.

The preselected representative is n=4096, replicate=1. Its trace/ACF figure shows
both smoothing parameters, log posterior and the coefficient with largest R-hat.
Its max R-hat is 1.00685, minimum ESS 995 and minimum E-BFMI 0.712. The spectrum
figure compares analytic S and posterior median S on identical linear color limits.
The median was checked to lie within its pointwise 90% posterior interval.

A plotting array-shape error was corrected; completed samples were preserved
and the remaining fits continued. The final driver smoke run completed and produced
all three figures. R numerical tests and source parsing passed; PNGs were visually
inspected for labels, layout and shared scales. PDFs were also exported.

## Earlier independent-TF example (superseded below) and README assets

Saved figures/metric tables from the repeated study are now in `assets/`, with
source-run provenance. The standalone `examples/blocked_joint.R` uses independent
complex TF cells, a Gaussian frequency ridge and an unknown tensor-spline variance.
It uses one exact amplitude Gibbs block and three symmetric MH blocks, including
log smoothing precisions with determinant/Jacobian terms. It does not require Stan.

Independent normalized-density checks passed for the complex likelihood and prior
ratio; numerical quadrature verified the Gaussian amplitude conditional. The
initial 10000-draw run missed the gate for one coefficient (R-hat 1.0114); the
20000-draw run passes across all parameters: max R-hat 1.00684, min bulk/tail ESS
1105.22. Sampling/summaries took 8.5 seconds for four sequential chains. Plotting
uses a display subset of draws; diagnostics use all retained draws. The printed
acceptance fractions exclude warmup and all adaptation stops at the warmup boundary.

Source output: `results/blocked-20260910-162626/`. Signal and PSD figures were
visually checked; normalized axes, common row color scales and true/median labels
are explicit. The injected amplitude lies just outside the 90% interval in this
single realization, so these assets make no empirical coverage claim. This toy
uses coefficient variance units and a direct TF generator; no moving-transform
or physical chirp likelihood has been validated by this exercise.


## LS2 time-domain noise plus high-SNR chirp

The joint example now generates the same LS2 process as the noise-only study,
adds an actual cosine chirp in the time domain, and applies identical cached
moving Fourier windows to data and every signal proposal. n=1024, dt=1,
m=16, thinning=2; 8/6 interior time/frequency knots. The known phase is zero.
The true amplitude 1.93989 gives exact LS2 optimal SNR 40 using the tridiagonal
covariance. The moving-Whittle approximation gives SNR 40.89.

Tests passed for cached-transform agreement/linearity, a dense-covariance check
of the exact SNR norm, the complex-template Gaussian amplitude conditional, and
the centered/non-centered MH rescaling Jacobian. The first real-LS2 run showed
slow frequency-smoothing mixing. A second smoothing update now jointly rescales
coefficients at fixed standardized coefficients, preserving the joint target.

Final run: `results/blocked-ls2-20260910-163746/`. Four sequential chains,
2000 warmup + 20000 retained draws each, max R-hat 1.00683, minimum bulk/tail
ESS 1086.02; sampling/summaries took 70.7 seconds. Every parameter passed.

| Parameter | Injection | Posterior mean | 90% interval |
|---|---:|---:|---:|
| A | 1.93989 | 2.02885 | [1.95224, 2.10617] |
| f0 [Hz] | 0.12000 | 0.119986 | [0.119942, 0.120031] |
| fdot [Hz/s] | 0.00012000 | 0.000120035 | [0.000119921, 0.000120149] |

Amplitude is strongly separated from zero but its 90% interval misses the true
amplitude in this realization. We report that outcome without changing the seed
or claiming interval calibration. Overlapping complex moving coefficients and
local-Whittle assumptions remain approximate. The initialization search uses
only noisy data and the prior box, with local phase-aware refinement. This is
not a Bayes-factor calculation or a global posterior-mode validation.

The replacement README assets show signal Fourier power and LS2 PSD, not the
previous artificial TF ridge. Integer-rounded display centers are evaluated
at their actual times; rendering supports their slightly nonuniform spacing.
The signal/PSD and diagnostic PNGs were inspected for units, scales and layout.


## NIMBLE NUTS for PSD, block Metropolis for signal

Run: `results/nimble-ls2-20260910-170348/`, using NIMBLE 1.4.3 and nimbleHMC 0.2.5.
Four sequential chains, seeds 501–504, 1000 warmup + 3000 retained iterations each.
Both samplers stop adapting after warmup. NUTS updates all 120 spline coefficients
and two smoothing precisions; built-in `RW_block` updates amplitude and the two
cycle-scaled chirp parameters. No custom sampling algorithm is implemented.

The shared input helper exactly reproduces the earlier blocked run's data and
injection. Independent R versus NIMBLE checks of the complete normalized target
passed before and after compilation and after changing both parameter blocks;
maximum absolute difference was 1.25e-11. Existing numerical tests also passed.
All 125 monitored parameters passed the convergence screen: maximum R-hat 1.00445,
minimum bulk/tail ESS 983.42, zero post-warmup divergences. Signal acceptance rates
were 0.305, 0.263, 0.355 and 0.321. Sampling took 155.9 seconds; model setup and
compilation took 33.7 seconds. No E-BFMI/tree-depth gate is claimed for this interface.

The posterior means of A, f0, fdot and both smoothing precisions differ from the
previous blocked run by at most 0.067 of that run's corresponding posterior SD.
This is a same-data consistency check, not a calibration study. The 90% intervals
are A [1.95096, 2.10383], f0 [0.1199427, 0.1200291] Hz and fdot
[0.0001199250, 0.0001201436] Hz/s. The amplitude interval again misses the injected
1.93989; no seed or target was changed to hide this outcome. The approximate
moving-Whittle likelihood and local initialization caveats still apply.

Signal/PSD comparison and signal trace plots were visually inspected. The
implementation favors existing library samplers and a small model definition;
this run was slower than the earlier hand-written blocked sampler, whose draw
budget differs. No claim of a general speed advantage is made.


## Fit-quality and runtime audit

Audit scripts/logs: `results/nimble-audit-20260910/`; the fitted model was unchanged.
The saved NIMBLE posterior mean signal has amplitude 4.43% above the injection,
maximum phase offset 0.00415 cycles, and whitened residual signal SNR 1.84
(compared with injected SNR 40). A local maximum-likelihood fit using the exact,
known LS2 time-domain covariance gives A=2.02740, close to the joint posterior
mean (about 2.02574). Thus this realization's amplitude offset is not evidence
of a NIMBLE-specific signal-fitting failure.

The posterior-mean log PSD has grid RMSE 0.35816 and pointwise 90% interval
coverage 0.82. This is a single-realization coverage fraction, not a calibrated
repeated-simulation coverage estimate. With both smoothing precisions fixed at
the fitted posterior means, penalized MAP fits give log-RMSE 0.37023 after
subtracting the exact injected signal and 0.36594 after subtracting the estimated
signal. Their surfaces differ by log-RMS 0.01815. All MAP optimizations converged.
This check does not refit the smoothing posterior, but suggests that signal
subtraction is not the main source of the PSD discrepancy.

For context, the same fixed-smoothing MAP fit to the exact expected moving noise
power (removing realization fluctuations) has log-RMSE 0.14645. The best
least-squares representation of the log truth in the display-grid spline space
has error 0.09093; the finite-window expectation differs from the instantaneous
truth by log-RMS 0.06056. These are diagnostic comparisons, not an additive error
decomposition. The fit has 496 retained complex observations and 120 spline
coefficients. Both finite-data fluctuations and smoothing/representation bias
remain relevant; passing MCMC diagnostics did not establish accurate PSD recovery.

A separate one-chain timing probe (1000 warmup, 500 retained steps, unchanged
samplers) used NIMBLE's `run(time=TRUE)` and `getTimes()`: PSD NUTS 5.12145 s,
signal block Metropolis 0.04129 s (99.2% in PSD sampling). The last inspected NUTS
transition used 31 leapfrog steps; that is not an average over the chain. This
supports gradient/trajectory work in the PSD block as the dominant sampler cost.
The full four-chain fit previously took 155.9 s sampling plus 33.7 s setup/compile.


## Faster NIMBLE evaluation and n=2048 demonstration

`R/nimble_joint.R` now evaluates the same normalized complex-normal likelihood
as one vectorized observation node, and samples standardized Gaussian coefficients
z with c=z/sqrt(q(phi)). NUTS and block Metropolis are still NIMBLE's built-in
samplers. The independent density check includes the transformation Jacobian;
maximum absolute target error across the six completed fits was 4.39e-11.
Both blocks still stop adapting after warmup. Numerical tests passed.

At the original n=1024, seed=4821 dataset, sampling fell from 155.9 to 57.7 seconds
with the same 1000 warmup/3000 retained budget and four sequential chains. Maximum
posterior mean shift across the 125 shared quantities was 0.084 of the previous
posterior SD. Log-PSD RMSE stayed at 0.358. The minimum ESS was 857 (previously 983),
so the gain is not merely a reduction in useful samples. Including model setup
and compilation, the recorded time fell from about 190 to 97 seconds. Compilation,
input preparation, summaries and plotting are not included in sampling time.

The revised example defaults to n=2048, retains 8/6 uniform interior knots and
SNR 40, and saves pointwise interval slices plus single-realization accuracy
metrics. `R/joint_summary.R` uses all retained draws for those summaries. The
physical f0 and fdot are unchanged; the longer observation extends the track.

`results/nimble-study-20260910-233426/` contains the predeclared three seeds
4821–4823 at each size, run with the revised implementation. All six fits passed:
maximum R-hat <=1.00601, minimum bulk/tail ESS >=821, and no post-warmup divergences.
The check monitors both standardized coefficients and the physical coefficients.

| n | Mean log-PSD RMSE | Mean pointwise coverage | Mean log interval width | Mean sampling time |
|---|---:|---:|---:|---:|
| 1024 | 0.34016 | 0.83625 | 0.95250 | 58.6 s |
| 2048 | 0.23627 | 0.85750 | 0.71038 | 121.8 s |

The larger datasets reduced mean error by 30.5% in this small check. Coverage
remains below nominal 90%; neither these three realizations nor successful MCMC
diagnostics establish interval calibration. No interval inflation, seed selection,
or truth-dependent fitting was used. All cases are retained in the metrics CSV.

The representative case was fixed in advance: n=2048, seed=4821. It has log-PSD
RMSE 0.25233, coverage 0.83125 and residual signal SNR 2.116 for the waveform at
posterior-mean parameters (injected SNR 40). The updated README heatmaps, traces
and six PSD slices come from this case. The slices explicitly expose intervals
that miss the truth. Plots were visually inspected; their sampling-check label
refers to MCMC computation, not statistical accuracy.
