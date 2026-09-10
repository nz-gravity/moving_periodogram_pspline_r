# Independent simulations are the sampling units for Monte Carlo uncertainty.
# Approximate t intervals describe the mean metric, not pointwise PSD uncertainty.
study_summary <- function(results) {
  do.call(rbind, lapply(sort(unique(results$n)), function(n) {
    do.call(rbind, lapply(c("rmse", "coverage", "width", "seconds"), function(key) {
      x <- results[results$n == n, key]; r <- length(x)
      se <- if (r > 1) sd(x)/sqrt(r) else NA_real_
      half <- if (r > 1) qt(0.95, r-1)*se else NA_real_
      data.frame(n = n, metric = key, replicates = r, mean = mean(x), se = se,
                 lower = max(0, mean(x)-half),
                 upper = min(if (key == "coverage") 1 else Inf, mean(x)+half))
    }))
  }))
}

plot_study <- function(results, path) {
  if (any(!results$ok)) stop("Diagnostics failed: inspect metrics.csv before plotting")
  sizes <- sort(unique(results$n))
  labels <- c("Log-spectrum RMSE [nats]", "Pointwise 90% coverage",
              "Mean log-interval width [nats]", "Wall time [s]")
  keys <- c("rmse", "coverage", "width", "seconds")
  if (endsWith(path, ".png")) png(path, width = 900, height = 1350, res = 150, bg = "white") else
    pdf(path, width = 6, height = 9, bg = "white")
  on.exit(dev.off())
  par(mfrow = c(4, 1), mar = c(1.8, 4.8, 0.5, 1), oma = c(3.5, 0, 2, 0))
  for (k in seq_along(keys)) {
    values <- lapply(sizes, function(n) results[results$n == n, keys[k]])
    center <- vapply(values, mean, numeric(1))
    half <- vapply(values, function(x) if (length(x) > 1)
      qt(0.95, length(x)-1) * sd(x) / sqrt(length(x)) else 0, numeric(1))
    lo <- pmax(0, center - half); hi <- center + half
    if (k == 2) hi <- pmin(1, hi)
    ylim <- range(lo, hi, if (k == 2) 0.9)
    if (diff(ylim) == 0) ylim <- ylim + c(-1, 1) * max(abs(ylim[1]) * 0.05, 0.01)
    plot(sizes, center, type = "n", log = if (k == 2) "x" else "xy",
         xlab = "", ylab = labels[k], xaxt = "n", ylim = ylim, las = 1)
    if (all(lengths(values) >= 2L))
      polygon(c(sizes, rev(sizes)), c(lo, rev(hi)), border = NA,
              col = adjustcolor("#2378A8", alpha.f = 0.22))
    lines(sizes, center, col = "#2378A8", lwd = 2)
    points(sizes, center, pch = 19, col = "#2378A8")
    axis(1, at = sizes, labels = if (k == 4) paste0(sizes / 1024, "k") else FALSE)
    if (k == 2) abline(h = 0.9, lty = 3, lwd = 2)
  }
  mtext("Number of observations n", side = 1, outer = TRUE, line = 0.3)
  title <- if (all(table(results$n) == 1L)) "LS2 development run: one realization per size" else
    paste0("LS2: ", min(table(results$n)), " replicates; 90% intervals for means")
  mtext(title, side = 3, outer = TRUE, line = 0.4, cex = 0.9)
}

# Saved post-warmup draws retain chain identity. Show smoothing, log posterior,
# and the coefficient with the largest R-hat; this is not cherry-picked mixing.
plot_convergence <- function(case, path) {
  cc <- case$summary[grepl("^c\\[", case$summary$variable), ]
  worst <- cc$variable[which.max(cc$rhat)]
  variables <- c("phi[1]", "phi[2]", "lp__", worst)
  draws <- posterior::as_draws_array(case$draws)[, , variables, drop = FALSE]
  if (endsWith(path, ".png")) png(path, 1500, 1500, res = 150, bg = "white") else
    pdf(path, 10, 10, bg = "white")
  on.exit(dev.off())
  par(mfrow = c(4, 2), mar = c(3, 4, 2, 1), oma = c(2.5, 0, 3, 0))
  colors <- c("#2378A8", "#D17520", "#407D50", "#9A529A")
  for (v in variables) {
    x <- matrix(draws[, , v], nrow = dim(draws)[1], ncol = dim(draws)[2])
    matplot(x, type = "l", lty = 1:4, col = colors, lwd = 0.7,
            xlab = "Post-warmup iteration", ylab = v, main = "Chain traces")
    ac <- sapply(seq_len(ncol(x)), function(j) as.numeric(acf(x[, j], lag.max = 50, plot = FALSE)$acf))
    matplot(0:50, ac, type = "l", lty = 1:4, col = colors,
            xlab = "Lag", ylab = "Autocorrelation", main = v, ylim = c(-0.2, 1))
    abline(h = 0, col = "grey70")
    if (v == variables[1]) legend("topright", paste("Chain", 1:4),
                                col = colors, lty = 1:4, bty = "n", cex = 0.8)
  }
  d <- case$diagnostics
  mtext("Representative fit: post-warmup convergence diagnostics", 3, outer = TRUE, line = 1)
  mtext(sprintf("Max R-hat %.3f | min ESS %.0f | divergences %d | depth hits %d | min E-BFMI %.2f",
                d$rhat, d$min_ess, d$divergences, d$treedepth, d$min_ebfmi),
        1, outer = TRUE, line = 0.5, cex = 0.85)
}

plot_spectra <- function(case, path) {
  s <- case$surface
  u <- sort(unique(s$u)); f <- sort(unique(s$omega)) / pi
  truth <- matrix(exp(case$truth), length(u), length(f))
  estimate <- matrix(s$median_psd, length(u), length(f))
  zlim <- range(truth, estimate)
  palette <- hcl.colors(100, "viridis")
  if (endsWith(path, ".png")) png(path, 1500, 700, res = 150, bg = "white") else
    pdf(path, 10, 14/3, bg = "white")
  on.exit(dev.off())
  layout(matrix(1:3, 1), widths = c(1, 1, 0.22))
  par(mar = c(4.5, 4.5, 2.5, 0.5))
  for (i in 1:2) {
    image(u, f, if (i == 1) truth else estimate, zlim = zlim, col = palette,
          xlab = "Rescaled time u", ylab = expression(omega/pi),
          main = if (i == 1) "True S(t, f)" else "Posterior median S(t, f)", useRaster = TRUE)
  }
  par(mar = c(4.5, 0.5, 2.5, 4))
  zz <- seq(zlim[1], zlim[2], length.out = 100)
  image(1:2, zz, matrix(rep(zz, each = 2), 2), col = palette,
        axes = FALSE, xlab = "", ylab = "", main = "S")
  axis(4, las = 1); box()
  mtext("Density / radian", 4, line = 2.6, cex = 0.8)
}
