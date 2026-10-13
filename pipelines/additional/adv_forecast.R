# ================================================================================
#  5. FORECAST — out-of-sample SARIMA backtest
# ================================================================================
#  For every virus x geography the last 12 months are held out; a seasonal ARIMA
#  (1,1,1)(1,1,0)[12] is fitted to log(1 + incidence) on the remaining months and
#  forecasts the hold-out (back-transformed with expm1, 95% prediction interval).
#  The model is the statsmodels SARIMAX specification the published results used,
#  reproduced exactly in R: the same 26-state state-space form (13 differencing +
#  13 ARMA states), approximate-diffuse initialisation (P0 = 1e6 I), log-likelihood
#  summed after a 26-observation burn-in, no stationarity/invertibility constraints,
#  conditional-sum-of-squares starting values and L-BFGS-B with forward-difference
#  gradients. (stats::arima uses a different initialisation and rejects
#  non-stationary AR fits, so it would not reproduce these numbers.)
#  Outputs (outputs/<date>/additional/):
#    forecast_backtest_metrics.csv   MAE / RMSE / MAPE / PI coverage per geography
#    forecast_vs_actual.csv          national + regional hold-out predictions
#    forecast_national.png
# ================================================================================
# Locate the project root (the folder holding R/utils.R) from this file, Rscript --file= or the working dir.
.root <- local({
  cand <- c(unlist(lapply(rev(sys.frames()), function(e) get0("ofile", envir = e, inherits = FALSE))),
            sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), file.path(getwd(), "."))
  for (p in cand) {
    d <- dirname(normalizePath(p, winslash = "/", mustWork = FALSE))
    repeat { if (file.exists(file.path(d, "R", "utils.R"))) return(d); if (dirname(d) == d) break; d <- dirname(d) }
  }
  stop("Cannot find the project root (the folder containing R/utils.R)")
})
source(file.path(.root, "R", "utils.R"))
source(path_root("pipelines", "additional", "adv_common.R"))
quiet_library("KFAS")

HOLDOUT <- 12
MIN_TRAIN <- 36
S <- 12                       # seasonal period
K_STATES <- 26                # (d + D*s) + max(p + P*s, q + 1) = 13 + 13
INIT_VAR <- 1e6
BURN <- K_STATES

# ---- state-space form of SARIMA(1,1,1)(1,1,0)[12] (statsmodels SARIMAX layout) -------------------
ss_design <- c(1, rep(0, S - 1), 1, 1, rep(0, 12))
ss_matrices <- function(phi, theta, Phi) {
  Tm <- matrix(0, K_STATES, K_STATES)
  Tm[1, 1] <- 1; Tm[1, 1 + S] <- 1; Tm[1, 14] <- 1                   # first difference
  for (r in 3:13) Tm[r, r - 1] <- 1                                    # seasonal-difference states
  Tm[2, 13] <- 1; Tm[2, 14] <- 1
  for (r in 14:25) Tm[r, r + 1] <- 1                                   # ARMA companion block
  Tm[14:26, 14] <- c(phi, rep(0, 10), Phi, -phi * Phi)                 # (1 - phi L)(1 - Phi L^12)
  Rv <- numeric(K_STATES); Rv[14] <- 1; Rv[15] <- theta
  list(T = Tm, R = Rv)
}

# Kalman filter (conventional, H = 0, P0 = 1e6 I); returns the log-likelihood after the burn-in and
# the predicted state/covariance after the last observation.
kalman <- function(y, params) {
  m <- ss_matrices(params[1], params[2], params[3])
  Tm <- m$T; RQR <- tcrossprod(m$R) * params[4]
  Z <- ss_design
  a <- numeric(K_STATES); P <- diag(INIT_VAR, K_STATES)
  ll <- 0
  for (t in seq_along(y)) {
    if (!is.na(y[t])) {
      PZ <- drop(P %*% Z); Fv <- sum(Z * PZ); v <- y[t] - sum(Z * a)
      a <- a + PZ * (v / Fv)
      P <- P - tcrossprod(PZ) / Fv
      if (t > BURN) ll <- ll - 0.5 * (log(2 * pi) + log(Fv) + v * v / Fv)
    }
    a <- drop(Tm %*% a)
    P <- Tm %*% tcrossprod(P, Tm) + RQR
    P <- (P + t(P)) / 2
  }
  list(ll = ll, a = a, P = P, T = Tm, RQR = RQR)
}

forecast_kf <- function(kf, h) {
  a <- kf$a; P <- kf$P; mean <- var <- numeric(h)
  for (i in seq_len(h)) {
    mean[i] <- sum(ss_design * a); var[i] <- drop(crossprod(ss_design, P %*% ss_design))
    a <- drop(kf$T %*% a); P <- kf$T %*% tcrossprod(P, kf$T) + kf$RQR
  }
  list(mean = mean, var = var)
}

# Fast log-likelihood for the optimiser: the same filter in KFAS's compiled code. The burn-in is
# removed as loglik(all observations) - loglik(first BURN observations), which is exact because the
# filter is sequential (verified to equal kalman()$ll).
kfas_loglik_fn <- function(y) {
  mk <- function(yy) {
    Rm <- matrix(0, K_STATES, 1); Rm[14, 1] <- 1
    SSModel(yy ~ -1 + SSMcustom(Z = matrix(ss_design, 1), T = ss_matrices(0, 0, 0)$T, R = Rm, Q = matrix(1),
                  a1 = matrix(0, K_STATES, 1), P1 = diag(INIT_VAR, K_STATES), P1inf = matrix(0, K_STATES, K_STATES)), H = matrix(0))
  }
  full <- mk(y); burn <- mk(y[seq_len(BURN)])
  upd <- function(m, p) {
    m$T[14:26, 14, 1] <- c(p[1], rep(0, 10), p[3], -p[1] * p[3]); m$R[15, 1, 1] <- p[2]; m$Q[1, 1, 1] <- p[4]
    m
  }
  function(p) as.numeric(stats::logLik(upd(full, p), check.model = FALSE)) - as.numeric(stats::logLik(upd(burn, p), check.model = FALSE))
}

# statsmodels lagmat(x, maxlag)[, maxlag] with trim="forward" (zero-padded at the start)
lag_vec <- function(x, k) c(rep(0, k), head(x, -k))
ls_fit <- function(X, Y) { b <- qr.coef(qr(X), Y); b[is.na(b)] <- 0; b }

# Conditional-sum-of-squares starting values (SARIMAX.start_params)
start_params <- function(y) {
  e <- y[(S + 1):length(y)] - y[1:(length(y) - S)]                     # seasonal difference
  e <- diff(e); e <- e[!is.na(e)]                                      # then first difference
  m <- length(e)
  # AR(1) + MA(1): long-AR(2) residuals as the MA regressor
  Y1 <- e[3:m]; X1 <- cbind(e[2:(m - 1)], e[1:(m - 2)])
  res1 <- Y1 - drop(X1 %*% ls_fit(X1, Y1))
  Y <- e[4:m]; X <- cbind(lag_vec(e, 1)[4:m], lag_vec(res1, 1)[2:length(res1)])
  b <- ls_fit(X, Y)
  res <- Y - drop(X %*% b)
  sigma2 <- max(mean(res[-1]^2), 1e-10)
  # seasonal AR(1)
  Ys <- e[(S + 1):m]; Xs <- matrix(lag_vec(e, S)[(S + 1):m])
  c(b[1], b[2], ls_fit(Xs, Ys)[1], sigma2)
}

# Maximum likelihood on -loglike/nobs in the unconstrained space (sigma2 = u4^2), as statsmodels
# parameterises it: L-BFGS-B from the CSS start with central-difference gradients and tight
# tolerances, then a BFGS polish; the better of the two is kept. (The likelihood surface is flat
# along the AR/MA ridge, so a loosely converged optimiser stops at different points; converging
# fully lands on the same maximum statsmodels reports.)
fit_sarima <- function(y) {
  nobs <- length(y)
  p0 <- start_params(y)
  u0 <- c(p0[1:3], sqrt(p0[4]))
  to_params <- function(u) c(u[1:3], u[4]^2)
  loglik <- kfas_loglik_fn(y)
  ll <- function(u) loglik(to_params(u))
  fn <- function(u) { v <- -ll(u) / nobs; if (is.finite(v)) v else 1e10 }
  gr <- function(u) {
    h <- .Machine$double.eps^(1 / 3) * pmax(abs(u), 0.1)
    -vapply(seq_along(u), function(k) { a <- u; b <- u; a[k] <- a[k] + h[k]; b[k] <- b[k] - h[k]; (ll(a) - ll(b)) / (2 * h[k]) }, 0) / nobs
  }
  o1 <- suppressWarnings(stats::optim(u0, fn, gr, method = "L-BFGS-B", control = list(maxit = 500, factr = 10, pgtol = 0, lmm = 10)))
  o2 <- tryCatch(stats::optim(o1$par, fn, gr, method = "BFGS", control = list(maxit = 500, reltol = 1e-14)), error = function(e) o1)
  to_params(if (o2$value < o1$value) o2$par else o1$par)
}

# (actual, pred, lo, hi) for the last HOLDOUT points, fitted on the rest; NULL when too short.
fit_forecast <- function(y) {
  n <- length(y)
  train <- y[1:(n - HOLDOUT)]; test <- y[(n - HOLDOUT + 1):n]
  if (sum(!is.na(train)) < MIN_TRAIN || sum(!is.na(test)) < HOLDOUT %/% 2) return(NULL)
  ytr <- log1p(train)
  params <- fit_sarima(ytr)
  fc <- forecast_kf(kalman(ytr, params), HOLDOUT)
  z <- stats::qnorm(0.975)
  list(actual = test, pred = expm1(fc$mean), lo = expm1(fc$mean - z * sqrt(fc$var)), hi = expm1(fc$mean + z * sqrt(fc$var)))
}

metrics <- function(actual, pred, lo, hi) {
  m <- is.finite(actual) & is.finite(pred)
  a <- actual[m]; p <- pred[m]
  if (!length(a)) return(NULL)
  err <- p - a; nz <- a > 0
  c(MAE = mean(abs(err)), RMSE = sqrt(mean(err^2)), MAPE_pct = if (any(nz)) mean(abs(err[nz] / a[nz])) * 100 else NA_real_,
    PI95_coverage_pct = mean(a >= lo[m] & a <= hi[m]) * 100)
}

fit_forecast_safe <- function(y) {
  if (length(y) < MIN_TRAIN + HOLDOUT) return(NULL)
  tryCatch(fit_forecast(y), error = function(e) NULL)
}

# The series are independent: fit them on a local worker pool (one R process per core, at most 8),
# falling back to a sequential loop if the pool cannot be started.
fit_all <- function(ys) {
  n_workers <- max(1L, min(8L, parallel::detectCores() - 1L, length(ys)), na.rm = TRUE)
  cl <- if (n_workers > 1) tryCatch(parallel::makeCluster(n_workers), error = function(e) NULL) else NULL
  if (is.null(cl)) return(lapply(ys, fit_forecast_safe))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(f) { source(f); invisible(NULL) }, path_root("pipelines", "additional", "adv_forecast.R"))
  cat(sprintf("  fitting %d series on %d worker processes ...\n", length(ys), n_workers))
  parallel::parLapplyLB(cl, ys, fit_forecast_safe)
}

run <- function() {
  datasets <- load_datasets()
  metric_rows <- list(); va_rows <- list()
  series <- iter_series(datasets)
  fits <- fit_all(lapply(series, function(it) it$s$Incidence))
  for (i in seq_along(series)) {
    it <- series[[i]]
    res <- fits[[i]]
    if (is.null(res)) next
    mt <- metrics(res$actual, res$pred, res$lo, res$hi)
    if (is.null(mt)) next
    metric_rows[[length(metric_rows) + 1]] <- data.table(Virus = it$virus, Level = it$level, Geo = it$geo,
      MAE = mt[["MAE"]], RMSE = mt[["RMSE"]], MAPE_pct = mt[["MAPE_pct"]], PI95_coverage_pct = mt[["PI95_coverage_pct"]],
      n_test = length(res$actual))
    if (it$level %in% c("national", "region"))
      va_rows[[length(va_rows) + 1]] <- data.table(Virus = it$virus, Level = it$level, Geo = it$geo,
        Date = tail(it$s$Date, HOLDOUT), Actual = res$actual, Pred = res$pred, PI_low = res$lo, PI_high = res$hi)
  }
  mets <- rbindlist(metric_rows); va <- rbindlist(va_rows)
  save_table(mets, "forecast_backtest_metrics.csv")
  save_table(va, "forecast_vs_actual.csv")
  for (r in seq_len(nrow(mets[Level == "national"]))) {
    x <- mets[Level == "national"][r]
    cat(sprintf("    %-12s national  MAPE=%5.1f%%  PI95 coverage=%4.0f%%\n", x$Virus, x$MAPE_pct, x$PI95_coverage_pct))
  }

  # Figure: national forecast vs actual
  panels <- lapply(seq_along(datasets), function(r) {
    virus <- names(datasets)[r]
    y <- monthly_series(datasets[[virus]], "national")$Brazil
    d <- va[Virus == virus & Level == "national"]
    p <- ggplot() + geom_line(data = y, aes(Date, Incidence, colour = "Actual"), linewidth = 1.2 / 2.13)
    if (nrow(d)) p <- p + geom_ribbon(data = d, aes(Date, ymin = PI_low, ymax = PI_high, fill = "95% PI"), alpha = 0.18) +
      geom_line(data = d, aes(Date, Pred, colour = "Forecast"), linewidth = 1.8 / 2.13)
    p <- p + scale_colour_manual(values = c(Actual = "#2c3e50", Forecast = "#e67e22"), breaks = c("Actual", "Forecast"), name = NULL) +
      scale_fill_manual(values = c("95% PI" = "#e67e22"), name = NULL) +
      labs(title = virus, x = if (r == length(datasets)) "Date" else NULL,
           y = "per 100k") + theme_adv()
    if (r == 1) p + theme(legend.position = "inside", legend.position.inside = c(0.01, 0.98), legend.justification = c(0, 1),
                          legend.text = element_text(size = 8), legend.background = element_rect(colour = "grey80", fill = "white"),
                          legend.spacing.y = unit(0, "pt"), legend.margin = margin(1, 4, 1, 4))
    else p + theme(legend.position = "none")
  })
  caption_note("forecast_national.png", "Former title", "<virus> — national 12-month hold-out forecast")
  save_fig(wrap_plots(panels, ncol = 1), "forecast_national.png", 10, 2.6 * length(datasets))
  write_caption_notes("adv_forecast.R")
  cat("  Forecasting done.\n")
}

if (sys.nframe() == 0L) run()
