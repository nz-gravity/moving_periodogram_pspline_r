# Experimental joint signal and PSD with Stan

**Draft, not yet validated:** the first sampling run was stopped when the design
changed to alternating PSD and signal updates. This file describes the optional
all-parameter NUTS implementation, not the recommended blocked workflow.

From the repository root, after `Rscript setup.R`, run:

```sh
Rscript examples/stan_joint.R
```

The driver simulates the same LS2 noise realization and SNR-40 chirp as
`blocked_joint.R`. It builds the Fourier weights and uniform tensor spline basis,
initializes the chirp from the noisy data, and runs four joint NUTS chains.
`joint.stan` contains the complete joint model. R handles simulation and plotting;
Stan differentiates through the waveform and its moving Fourier transform.

The real time-domain signal is
`A*cos(2*pi*(f0*t + 0.5*fdot*t^2))`, with known phase zero, dt=1 second,
and n=1024. Frequency is in Hz and its derivative in Hz/s. Internally Stan samples
`n*f0` and `n^2*fdot`, which are better scaled than the physical parameters.
The likelihood is `-sum(log(S) + |D-H|^2/S)`, where D and H are the complex
moving transforms of the data and proposed signal, and `log(S)=B*c`.
This preserves signal phase; subtracting signal power from a periodogram would
not give this likelihood.

The priors match the blocked example: A ~ Normal(0,5), f0 uniform on
[0.08,0.18] Hz, fdot uniform on [0,0.0002] Hz/s, and two Gamma(2,1)
smoothing precisions controlling a Gaussian spline prior. See [BLOCKED.md](BLOCKED.md)
for the penalty and LS2/SNR definitions. NUTS samples all these parameters together;
there are no custom Gibbs or MH steps in the Stan model.

The script saves `case.rds`, `summary.csv`, `session.txt`, `stan_data.rds`,
`log_lik.rds`, `joint.png`, `signal_diagnostics.png` and `psd_diagnostics.png`
in a new `results/stan-ls2-TIMESTAMP/` directory. `case.rds` includes the NUTS
report: R-hat, bulk/tail ESS, divergences, tree-depth hits and E-BFMI.
The gate checks both signal and PSD parameters. Inspect it before interpreting
recovery. Edit `warmup` and `samples` in the driver to change the budget;
reduce `parallel_chains` if memory is limited, retaining four diagnostic chains.

Use this example when you want automatic gradient-based joint sampling. The
blocked example remains useful for learning conditional updates and running
without a C++ toolchain. Their wall times depend on chain parallelism and draw
budgets; compare effective samples and diagnostics, not just raw draw counts.

This is a local parameter-estimation demonstration initialized with a data-only
search. High SNR makes the frequency posterior narrow, so initialization matters;
four nearby chains do not establish that all remote modes were explored.
The moving-Whittle likelihood is approximate, and a single successful fit does
not establish repeated-data coverage or a detection Bayes factor.

Optional independent R-versus-Stan likelihood check after a run:

```sh
Rscript examples/check_stan_joint.R results/stan-ls2-TIMESTAMP
```

This evaluates the Stan likelihood at a saved draw with full precision and
compares it with the R complex-transform likelihood. It needs no further MCMC.
