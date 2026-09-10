# Saved README figures

These files are deliberately outside Git-ignored `results/` so they can be shared
with the source. They are fixed snapshots, not automatically updated by new runs.
No posterior chains, compiled binaries or machine-specific runtimes are included.

## Moving-periodogram study

`scaling.png`, `spectra.png`, `convergence.png`, `metrics.csv`, `summary.csv`:
source run `results/study-20260910-141441/`, 10 September 2026.
Eight independent simulations at n=512,1024,2048,4096; seed base 300000;
four NUTS chains with 750 warmup and 750 retained draws each. One same-dataset
extension used 1500/1500 iterations and is recorded in `metrics.csv`.
The representative fit is n=4096, replicate=1, chosen in advance.
Ribbons are approximate 90% t confidence intervals for replicate means.

Regenerate using `Rscript demo.R study 8 300000`; inspect diagnostics before copying
PNGs and CSVs from the newly printed output folder to `assets/`. Software versions,
timing and Monte Carlo variation can change outputs. Tested versions: R 4.5.3,
CmdStanR 0.9.0, CmdStan 2.39.0, posterior 1.7.0. See `../VALIDATION.md`.

## Joint LS2/chirp blocked example

`blocked_joint.png`, `blocked_signal_diagnostics.png`, `blocked_summary.csv`,
`blocked_injection.csv` (parameters and SNR definitions):
source run `results/blocked-ls2-20260910-163746/`, simulation seed 4821,
chain seeds 501–504, 2000 warmup + 20000 retained iterations per chain.
These replace the earlier independent-TF toy snapshots.

The data are 1024 LS2 samples plus a phase-zero linear chirp with f0=0.12 Hz,
fdot=0.00012 Hz/s, amplitude 1.93989 (dt=1 second). The exact time-domain LS2
optimal SNR is 40; the approximate moving-Whittle SNR is 40.89. Inference uses
m=16, thinning=2, 8/6 interior time/frequency knots and a joint PSD/signal sampler.
All-parameter max R-hat=1.00683, min bulk/tail ESS=1086.02. Sampling/summaries took
70.7 seconds here. Initialization is data-driven; the truth is not supplied to it.

Regenerate with `Rscript examples/blocked_joint.R`. The PNGs are `joint.png` and
`signal_diagnostics.png` in that run, renamed here for README embedding. The signal
images display transformed waveform power; S is density per radian. PSD/interval
accuracy is not established by this one run, even though the signal is loud.

## NIMBLE alternating sampler

The `nimble_*.png` files and `nimble_summary.csv`, `nimble_diagnostics.csv`,
`nimble_injection.csv` now show the preselected n=2048, seed=4821 case from
`results/nimble-study-20260910-233426/`. The earlier n=1024 NIMBLE fit remains
documented in `../VALIDATION.md`. Regenerate the current example with
`Rscript examples/nimble_joint.R`; use `Rscript examples/nimble_study.R` for
all three seeds at each size. The median heatmap is accompanied by
`nimble_psd_slices.png`, which shows pointwise 90% intervals on fixed slices.

The signal has A=1.37621, f0=0.12 Hz, fdot=0.00012 Hz/s, known phase zero,
and dt=1 second. Exact LS2 optimal SNR is 40; moving-Whittle SNR is 41.80.
The model retains 8/6 interior time/frequency knots and 120 spline coefficients.
Four sequential chains used seeds 501–504, 1000 warmup + 3000 retained iterations.
Maximum R-hat 1.00483, minimum bulk/tail ESS 821.09, no post-warmup divergences.
Sampling took 119.2 s; model setup/compilation took 32.4 s.

This case's log-PSD RMSE is 0.2523 and pointwise coverage is 0.8313. These
figures show improved estimation, not calibrated intervals. R 4.5.3, NIMBLE 1.4.3,
nimbleHMC 0.2.5 and posterior 1.7.0 were used. Compilation and posterior summaries
can vary with software and platform.

`nimble_study_metrics.csv` retains all six fits; `nimble_study_summary.csv` gives
means over the three seeds at each n. Mean log-PSD error is 0.3402 at n=1024 and
0.2363 at n=2048; mean coverage is 0.8363 and 0.8575 respectively. This is a small
repeat check, not evidence of calibrated intervals.
