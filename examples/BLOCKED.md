# LS2 noise + a high-SNR linear chirp: Gibbs and blocked MH

Run from the repository root:

```sh
Rscript examples/blocked_joint.R
```

Requires **R and `posterior` only** (`install.packages("posterior")` in R).
No Stan or compiler is used. This uses `simulate_ls2()` and the same thinned moving
periodogram as the noise-only study; it is no longer a direct synthetic TF ridge.

## Data and SNR

Generate n=1024 samples, dt=1 second:

```
x[t] = LS2_noise[t] + A*cos(2*pi*(f0*t + 0.5*fdot*t^2))
t = 0, ..., n-1
```

The phase at t=0 is fixed to zero. `f0` is in Hz and `fdot` in Hz/s, so the
instantaneous frequency is `f0 + fdot*t`. Injection: f0=0.12 Hz and
fdot=0.00012 Hz/s. The whole track is below the 0.5 Hz Nyquist frequency.
Changing n, dt or the prior requires rechecking that condition.

The amplitude is chosen for **optimal time-domain SNR 40** using the known LS2
covariance C: `rho^2 = h' solve(C) h`. For LS2, C is tridiagonal with diagonal
`1+a[t]^2` and off-diagonal `C[t,t-1]=a[t]`, where
`a[t]=1.1*cos(1.5-cos(4*pi*t/n))`. `ls2_whiten()` evaluates that quadratic form
by a banded Cholesky solve, including the initial latent innovation.
This uses the injection model to set its strength, not to estimate the signal.

The script also prints the **approximate moving-Whittle SNR**,
`sqrt(2*sum(abs(H)^2/S_true))`, on the retained ordinates. It need not equal the
optimal SNR because the windows overlap and the likelihood is approximate.
Both definitions are saved. This is a deliberately loud parameter-estimation
example, not a calibrated detection statistic or Bayes-factor calculation.

## The likelihood retains the waveform phase

Transform the data and every proposed time-domain waveform through identical
33-sample windows (`m=16`) with block thinning 2 and the same Fourier normalization:

```
d = moving_Fourier(x)
H(A,f0,fdot) = moving_Fourier(A*cos(phase(f0,fdot)))
log L = -sum(eta + abs(d-H)^2*exp(-eta))
eta = log(S) = B*c
```

`fourier_operator()` caches the window indices and Fourier weights; it is checked
against `moving_periodogram()` and for linearity. Signal subtraction occurs in
complex coefficient space, before taking powers. Subtracting signal power from
observed power would be a different and incorrect residual likelihood.

S uses the same **density per radian** convention as the original LS2 example.
The complex-normal/dynamic-Whittle independence and local-spectrum approximations
remain approximations for these overlapping LS2 windows. Passing MCMC diagnostics
does not validate frequentist coverage of this likelihood.

## Sampling blocks

1. **Amplitude: exact Gaussian Gibbs conditional for this likelihood.**
   With complex unit template g, prior A~Normal(0,25), its precision is
   `v_inv=1/25+2*sum(abs(g)^2/S)` and mean is
   `2*sum(Re(Conj(g)*d)/S)/v_inv`. Variance is `1/v_inv`.
2. **(f0,fdot): correlated Gaussian random-walk MH.** Uniform prior
   f0 in (0.08,0.18) Hz and fdot in (0,0.0002) Hz/s. Internal coordinates
   `(n*f0,n^2*fdot)` keep proposal scales comparable; saved draws use physical units.
3. **PSD: blocked MH in preconditioned coordinates.** Uniform cubic splines with
   8 time/6 frequency interior knots give 120 coefficients. A fixed approximate
   Fisher matrix sets proposal geometry, then blocks of 12 coordinates move in
   turn. Every move uses its full likelihood and prior ratio.
4. **Smoothing: MH plus a joint coefficient-rescaling move.** Both smoothing
   precisions have Gamma(2,1) priors. First update log(phi) with coefficients fixed.
   Then propose log(phi) holding standardized coefficients `z=sqrt(q)*c` fixed,
   so `c_new=c_old*sqrt(q_old/q_new)`. This additional MH move reduces coupling.
   Its Gaussian determinant cancels the coefficient-rescaling Jacobian; the ratio
   contains the likelihood change and the Gamma priors with log-phi Jacobians.

Here `q=phi_t*lambda_t+phi_f*lambda_f+1e-6` off the joint null space; the four
null directions have precision 0.01. This proper weak null prior differs from the
main NUTS example's 1e-4. The centered prior retains its full `sum(log(q))/2`
normalizer. Independent conjugate Gamma updates would be wrong for the additive
precision. All these conditional kernels target the same joint posterior.

Only the amplitude block draws directly from its full conditional: the overall
algorithm is **Metropolis-within-Gibbs**. Proposal scales adapt only during 2,000
warmup iterations and then freeze for 20,000 retained iterations per chain.
Four chains run sequentially. Acceptance fractions exclude warmup; `psd` is the
fraction of individual PSD-block proposals accepted, averaged over blocks.

Initialization scans chirp rates across the prior box, FFTs the dechirped **noisy
data**, and refines the best candidate using the phase-aware moving transform.
It does not start at the injection. Chains start near this data-selected peak,
so this is not an independent demonstration of global multimodal exploration.

## Code and outputs

- `examples/blocked_joint.R`: LS2 data, SNR-scaled injection, data-only search,
  splines, sampler, checks and saves. Change the injection/budget here.
- `R/chirp.R`: waveform, cached Fourier operator, exact LS2 whitening and search.
- `R/blocked.R`: residual likelihood, priors, amplitude conditional and MH blocks.
- `R/blocked_plot.R`: signal-power/PSD comparisons and chain diagnostics.
- `tests.R`: transform linearity/normalization, exact covariance norm, complex
  amplitude conditional and prior/Jacobian checks against independent calculations.

Each `results/blocked-ls2-*/` folder contains data, signal/noise time series,
settings/SNRs, posterior draws and diagnostics (`case.rds`), plus `summary.csv`,
`session.txt`, `joint.png`, `signal_diagnostics.png`, `psd_diagnostics.png`.
The image displays a denser grid of actual window centers; fitting still uses
only the original scattered zigzag ordinates. Signal panels show **Fourier power**,
not a real-valued TF ridge. PSD panels share a color scale in density per radian.

Surface figures use up to 500 retained draws; displayed traces use up to 2,000
points per chain. Diagnostics use every retained draw and every parameter. Plots
can still be saved for a failing run, with `ok=FALSE`; do not interpret them as
proof of convergence. Refer to the current `VALIDATION.md` for measured results.
