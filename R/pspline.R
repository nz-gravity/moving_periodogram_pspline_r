# Cubic uniform B-splines on [0,1], with integrated squared second derivative.
# interior: nonnegative interior-knot count; returns knots, penalty, eigenvalues
# and eigenvectors. A cubic basis has interior+4 columns on [0,1].
spline_setup <- function(interior = 8L) {
  stopifnot(interior >= 0L, interior == as.integer(interior))
  knots <- c(rep(0, 4), if (interior) seq_len(interior) / (interior + 1), rep(1, 4))
  # Two-point Gauss quadrature is exact for products of linear second derivatives.
  edges <- unique(knots)
  left <- head(edges, -1); right <- tail(edges, -1)
  half <- (right - left) / 2; mid <- (right + left) / 2
  xx <- c(mid - half / sqrt(3), mid + half / sqrt(3))
  dd <- splines::splineDesign(knots, xx, ord = 4, derivs = 2)
  penalty <- crossprod(dd, dd * c(half, half))
  penalty <- penalty / sum(diag(penalty))
  ee <- eigen(penalty, symmetric = TRUE)
  values <- pmax(ee$values, 0)
  values[values <= 1e-10 * max(max(values), 1)] <- 0
  list(knots = knots, penalty = penalty, values = values, vectors = ee$vectors)
}

# Evaluate an N-by-K basis matrix for coordinates x in [0,1].
basis <- function(spline, x) splines::splineDesign(spline$knots, x, ord = 4)

# Flatten time first (R column-major); evaluate scattered observations without
# materializing the N x (Kt*Kf) tensor design.
# ordinates: dataframe from moving_periodogram; kt/kf: interior-knot counts.
# Returns time/frequency spline objects and a data list for model.stan.
prepare_model <- function(ordinates, kt = 8L, kf = 10L) {
  tt <- spline_setup(kt); ff <- spline_setup(kf)
  frequency <- sort(unique(ordinates$omega))
  lt <- rep(tt$values, length(ff$values))
  lf <- rep(ff$values, each = length(tt$values))
  list(time = tt, frequency = ff, data = list(
    N = nrow(ordinates), Kt = length(tt$values), Kf = length(ff$values),
    J = length(frequency), Bt = unname(basis(tt, ordinates$u) %*% tt$vectors),
    Bf = unname(basis(ff, frequency / pi) %*% ff$vectors),
    rung = match(ordinates$omega, frequency), y = ordinates$mi,
    lt = lt, lf = lf, is_null = as.integer(lt == 0 & lf == 0),
    null_precision = 1e-4, ridge = 1e-6, alpha_phi = 2, beta_phi = 1))
}

# prepared: prepare_model result; model: compiled cmdstanr model. Iteration
# counts are per chain. Returns a CmdStanMCMC object; does not summarize/check it.
fit_pspline <- function(prepared, model, seed = 1L, warmup = 500L,
                        samples = 500L, chains = 4L) {
  d <- prepared$data
  # Constant log-power starting surface, projected into the penalty eigenbasis.
  W <- matrix(log(mean(d$y)), d$Kt, d$Kf)
  C <- t(prepared$time$vectors) %*% W %*% prepared$frequency$vectors
  model$sample(data = d, seed = seed, chains = chains, parallel_chains = chains,
    iter_warmup = warmup, iter_sampling = samples, refresh = 0,
    adapt_delta = 0.9, max_treedepth = 12,
    init = function() list(c = as.vector(C) + rnorm(length(C), 0, 0.01),
                           phi = c(2, 2)))
}

# grid: dataframe of paired u and omega coordinates. Returns that dataframe
# plus posterior mean log S, 5%/95% log quantiles, and posterior median S.
summarize_surface <- function(fit, prepared, grid) {
  C <- posterior::as_draws_matrix(fit$draws("c"))
  C <- C[, sprintf("c[%d]", seq_len(ncol(C))), drop = FALSE]
  bt <- basis(prepared$time, grid$u) %*% prepared$time$vectors
  bf <- basis(prepared$frequency, grid$omega / pi) %*% prepared$frequency$vectors
  design <- do.call(cbind, lapply(seq_len(ncol(bf)), function(j) bt * bf[, j]))
  # Bound temporary storage when summarizing the common evaluation grid.
  result <- matrix(NA_real_, nrow(grid), 4)
  for (start in seq(1L, nrow(grid), by = 100L)) {
    ii <- start:min(start + 99L, nrow(grid))
    eta <- C %*% t(design[ii, , drop = FALSE])
    result[ii, ] <- cbind(colMeans(eta),
                          apply(eta, 2, quantile, 0.05),
                          apply(eta, 2, quantile, 0.95),
                          apply(eta, 2, median))
  }
  cbind(grid, mean_log = result[, 1], lower_log = result[, 2], upper_log = result[, 3],
        median_psd = exp(result[, 4]))
}

# Returns one row with worst-case R-hat, ESS, divergences, tree-depth hits,
# E-BFMI and an overall ok flag. Passing assesses sampling, not calibration.
fit_diagnostics <- function(fit) {
  s <- fit$summary(variables = c("c", "phi"))
  d <- fit$diagnostic_summary()
  rhat <- max(s$rhat); ess <- min(s$ess_bulk, s$ess_tail)
  divergences <- sum(d$num_divergent)
  treedepth <- sum(d$num_max_treedepth)
  ebfmi <- min(d$ebfmi)
  ok <- is.finite(rhat) && is.finite(ess) && is.finite(ebfmi) &&
    rhat < 1.01 && ess >= 100 && divergences == 0 && treedepth == 0 && ebfmi >= 0.3
  data.frame(rhat = rhat, min_ess = ess, divergences = divergences,
             treedepth = treedepth, min_ebfmi = ebfmi, ok = ok)
}
