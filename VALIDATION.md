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
