# ================================================================================
#  GLM COMMON HELPERS — Poisson / quasi-Poisson / Negative Binomial
# ================================================================================
#  Shared regression machinery for the three count-data GLM scripts.
#
#  Model: Cases ~ climate predictors + C(Month)  with  offset(log Population)
#  Month fixed effects absorb the seasonal cycle, so climate IRRs describe
#  departures from the usual seasonal pattern rather than seasonality itself.
#
#  Three model specifications per model family:
#    spec="per_state" : one GLM per (Disease, State) on its monthly time series
#                       Predictors: Temp_Avg, Rain_mm, C(Month)
#    spec="lagged"    : per (Disease, State) with lag-0/1/2 Temp and Rain, plus C(Month)
#    spec="pooled_fe" : one GLM per Disease, pooled across all states with
#                       State as a categorical fixed effect.
#                       Predictors: Temp_Avg, Rain_mm, C(Month), C(Location_Name)
#
#  Standard errors (column `cov_type`):
#    - Poisson / Negative Binomial, per-state and lagged : HAC (Newey-West, 3 lags,
#      Bartlett kernel, no small-sample adjustment) for within-state autocorrelation.
#    - Poisson / Negative Binomial, pooled_fe : cluster-robust by state (HC1 x G/(G-1)).
#    - Quasi-Poisson : model-based SEs inflated by sqrt(Pearson chi2 / df), with
#      t-distribution p-values and CIs (R's quasipoisson).
#
#  Per-state and lagged results carry Benjamini-Hochberg q-values
#  (`q_value`, `Significant_FDR`) within each (Virus, Term) family of 27 states.
#
#  Negative Binomial: NB2 GLM (Var = mu + alpha mu^2) with alpha estimated by profile
#  likelihood over log(alpha) in [1e-4, 100] (the estimator of MASS::glm.nb);
#  AIC/BIC count alpha as a parameter.
#  Term names follow the patsy convention of the original outputs
#  (Intercept, C(Month)[T.2], C(Location_Name)[T.Alagoas], Temp_Avg, ...).
#  Sourced after R/utils.R and pipelines/data_loader.R.
# ================================================================================


ALPHA <- 0.05
HAC_MAXLAGS <- 3
IRLS_MAXITER <- 1000
ALPHA_BOUNDS <- c(1e-4, 100)

BASE_PREDICTORS   <- c("Temp_Avg", "Rain_mm")
LAGGED_PREDICTORS <- c("Temp_Avg", "Temp_Lag1", "Temp_Lag2", "Rain_mm", "Rain_Lag1", "Rain_Lag2")
FORMULA_BASE   <- "Cases ~ Temp_Avg + Rain_mm + C(Month)"
FORMULA_LAGGED <- "Cases ~ Temp_Avg + Temp_Lag1 + Temp_Lag2 + Rain_mm + Rain_Lag1 + Rain_Lag2 + C(Month)"
FORMULA_FE     <- "Cases ~ Temp_Avg + Rain_mm + C(Month) + C(Location_Name)"

SPECS <- list(
  list(spec = "per_state", predictors = BASE_PREDICTORS,   formula = FORMULA_BASE,   pooled = FALSE, min_n = 24),
  list(spec = "lagged",    predictors = LAGGED_PREDICTORS, formula = FORMULA_LAGGED, pooled = FALSE, min_n = 30),
  list(spec = "pooled_fe", predictors = BASE_PREDICTORS,   formula = FORMULA_FE,     pooled = TRUE,  min_n = 50))

# Drop NA on predictors and dependent variables; require Population > 0.
# Row order (state, year, month) is preserved, which HAC SEs rely on.
prepare_glm_frame <- function(df, predictors) {
  cols <- c("Cases", "Population", "Location_Name", "Month", predictors)
  sub <- as.data.frame(df)[, cols]
  sub <- sub[stats::complete.cases(sub), , drop = FALSE]
  sub <- sub[sub$Population > 0, , drop = FALSE]
  sub$log_pop <- log(as.numeric(sub$Population))
  rownames(sub) <- NULL
  sub
}

# Design matrix with patsy-style column names: Intercept, C(Month)[T.m], C(Location_Name)[T.x], predictors
design_matrix <- function(sub, predictors, pooled) {
  X <- cbind(Intercept = rep(1, nrow(sub)))
  months <- sort(unique(sub$Month))
  for (mo in months[-1]) X <- cbind(X, as.numeric(sub$Month == mo))
  colnames(X)[-1] <- sprintf("C(Month)[T.%d]", months[-1])
  if (pooled) {
    locs <- sort_c(unique(sub$Location_Name))
    L <- vapply(locs[-1], function(l) as.numeric(sub$Location_Name == l), numeric(nrow(sub)))
    if (!is.matrix(L)) L <- matrix(L, nrow = nrow(sub))
    colnames(L) <- sprintf("C(Location_Name)[T.%s]", locs[-1])
    X <- cbind(X, L)
  }
  cbind(X, as.matrix(sub[, predictors, drop = FALSE]))
}

# ---- GLM engine -----------------------------------------------------------------------------------
# Log-link Poisson (alpha = 0) / NB2 (alpha > 0) GLM fitted exactly as statsmodels' GLM.fit does it, so
# results agree to rounding: IRLS from mu0 = (y + mean(y)) / 2, stopping when the deviance changes by
# <= 1e-8 (absolute), with a Newton fallback; robust covariances are sandwiches of the analytic score
# and the observed Hessian at the estimates.
glm_deviance <- function(y, mu, alpha) {
  ymu <- pmax(y / mu, .Machine$double.eps)
  if (alpha == 0) return(sum(2 * (y * log(ymu) - (y - mu))))
  ya <- y + 1 / alpha; ma <- mu + 1 / alpha
  sum(2 * (y * log(ymu) - ya * log(ya / ma)))
}
glm_loglik <- function(y, mu, alpha) {
  if (alpha == 0) return(sum(y * log(mu) - mu - lgamma(y + 1)))
  sum(y * log(alpha * mu) - (y + 1 / alpha) * log(1 + alpha * mu) + lgamma(y + 1 / alpha) - lgamma(1 / alpha) - lgamma(y + 1))
}
wls_solve <- function(X, z, w) {
  sw <- sqrt(w)
  qr.coef(qr(X * sw, tol = 1e-12), z * sw)
}

irls_fit <- function(X, y, offset, alpha = 0, maxiter = IRLS_MAXITER) {
  mu <- (y + mean(y)) / 2
  lin_pred <- log(mu)
  dev_old <- glm_deviance(y, mu, alpha)
  if (is.nan(dev_old)) stop("The first guess on the deviance function returned a nan.")
  converged <- FALSE; params <- rep(NA_real_, ncol(X)); w <- NULL
  for (it in seq_len(maxiter)) {
    w <- mu / (1 + alpha * mu)
    z <- lin_pred + (y - mu) / mu - offset
    params <- wls_solve(X, z, w)
    if (any(!is.finite(params))) stop("IRLS diverged")
    lin_pred <- drop(X %*% params) + offset
    mu <- exp(lin_pred)
    dev <- glm_deviance(y, mu, alpha)
    if (isTRUE(abs(dev - dev_old) <= 1e-8)) { converged <- TRUE; break }
    dev_old <- dev
  }
  list(params = setNames(params, colnames(X)), mu = mu, converged = converged, weights = w, method = "IRLS")
}

score_factor <- function(y, mu, alpha) (y - mu) / (1 + alpha * mu)
oim_factor <- function(y, mu, alpha) mu / (1 + alpha * mu) * (1 + alpha * score_factor(y, mu, alpha))

# statsmodels LikelihoodModel Newton: x <- x - H^-1 g on (-loglike/n), ridge 1e-10, |dx| <= 1e-8
newton_fit <- function(X, y, offset, alpha, start, maxiter = 200, tol = 1e-8, ridge = 1e-10) {
  n <- length(y); old <- rep(Inf, length(start)); new <- start; it <- 0
  while (it < maxiter && any(abs(new - old) > tol)) {
    mu <- exp(drop(X %*% new) + offset)
    g <- -crossprod(X, score_factor(y, mu, alpha)) / n
    H <- crossprod(X, X * oim_factor(y, mu, alpha)) / n
    diag(H) <- diag(H) + ridge
    old <- new
    new <- old - drop(solve(H, g))
    it <- it + 1
  }
  mu <- exp(drop(X %*% new) + offset)
  list(params = setNames(new, colnames(X)), mu = mu, converged = it < maxiter, weights = mu / (1 + alpha * mu), method = "newton")
}

glm_fit <- function(X, sub, alpha = 0) {
  f <- irls_fit(X, sub$Cases, sub$log_pop, alpha)
  f$alpha <- alpha
  f$llf <- glm_loglik(sub$Cases, f$mu, alpha)
  f
}
fit_poisson <- function(X, sub) glm_fit(X, sub, 0)

pearson_dispersion <- function(f, y, k) {
  df_resid <- length(y) - k
  if (df_resid > 0) sum((y - f$mu)^2 / (f$mu + f$alpha * f$mu^2)) / df_resid else NaN
}

# NB2 GLM with alpha by bounded profile likelihood (Brent, xatol 1e-6 on log alpha; R's optimize and
# scipy's minimize_scalar(method="bounded") run the same algorithm). IRLS fits that fail are redone
# by Newton from the Poisson solution.
fit_negbin <- function(X, sub) {
  start <- fit_poisson(X, sub)$params
  fit <- function(alpha) {
    f <- tryCatch(glm_fit(X, sub, alpha), error = function(e) NULL)
    if (!is.null(f) && f$converged && is.finite(f$llf)) return(f)
    f <- newton_fit(X, sub$Cases, sub$log_pop, alpha, start)
    f$alpha <- alpha
    f$llf <- glm_loglik(sub$Cases, f$mu, alpha)
    f
  }
  neg_llf <- function(log_alpha) {
    f <- tryCatch(fit(exp(log_alpha)), error = function(e) NULL)
    if (is.null(f) || !is.finite(f$llf)) return(Inf)
    -f$llf
  }
  lo <- log(ALPHA_BOUNDS[1]); hi <- log(ALPHA_BOUNDS[2])
  opt <- suppressWarnings(stats::optimize(neg_llf, c(lo, hi), tol = 1e-6))
  alpha <- exp(opt$minimum)
  f <- fit(alpha)
  at_bound <- min(opt$minimum - lo, hi - opt$minimum) < 1e-3
  f$profile_converged <- is.finite(opt$objective) && is.finite(f$llf) && !at_bound
  if (!(f$profile_converged && f$converged))
    cat(sprintf("    WARNING: NegBin fit did not converge (alpha=%.4g); estimates unreliable\n", alpha))
  list(m = f, alpha = alpha)
}

# Robust covariance from score_obs and the observed Hessian (statsmodels cov_hac / cov_cluster):
# HAC Newey-West (Bartlett, HAC_MAXLAGS, no correction); cluster with G/(G-1)*(N-1)/(N-K).
robust_vcov <- function(f, X, sub, pooled) {
  y <- sub$Cases
  xu <- X * score_factor(y, f$mu, f$alpha)
  Hinv <- solve(-crossprod(X, X * oim_factor(y, f$mu, f$alpha)))
  if (pooled) {
    g <- sub$Location_Name
    sums <- rowsum(xu, g, reorder = FALSE)
    S <- crossprod(sums)
    G <- nrow(sums); n <- nrow(xu); k <- ncol(xu)
    V <- Hinv %*% S %*% t(Hinv) * (G / (G - 1) * ((n - 1) / (n - k)))
  } else {
    S <- crossprod(xu)
    for (lag in seq_len(HAC_MAXLAGS)) {
      wl <- 1 - lag / (HAC_MAXLAGS + 1)
      s <- crossprod(xu[(lag + 1):nrow(xu), , drop = FALSE], xu[1:(nrow(xu) - lag), , drop = FALSE])
      S <- S + wl * (s + t(s))
    }
    V <- Hinv %*% S %*% t(Hinv)
  }
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

# Model-based covariance of the final IRLS step: (X' W X)^-1 with the weights of the last iteration
model_vcov <- function(f, X) {
  V <- solve(crossprod(X, X * f$weights))
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

# Per-coefficient rows. `scale` (quasi-Poisson) inflates model-based SEs and switches
# inference to the t distribution.
summarize_fit <- function(coefs, V, scale = NULL, df_resid = NA) {
  se <- sqrt(diag(V))
  if (is.null(scale)) {
    z <- coefs / se; p <- 2 * stats::pnorm(-abs(z)); crit <- stats::qnorm(0.975)
  } else {
    se <- se * sqrt(scale); z <- coefs / se
    p <- 2 * stats::pt(abs(z), df_resid, lower.tail = FALSE); crit <- stats::qt(0.975, df_resid)
  }
  data.table(Term = names(coefs), Coef = unname(coefs), SE = unname(se), IRR = exp(unname(coefs)),
             IRR_CI_L = exp(unname(coefs - crit * se)), IRR_CI_U = exp(unname(coefs + crit * se)),
             z = unname(z), p_value = unname(p), Significant = unname(p) < ALPHA)
}

information_criteria <- function(llf, k, n) list(AIC = -2 * llf + 2 * k, BIC = -2 * llf + k * log(n))

# One model; returns list(rows, aux, llf, converged)
fit_family <- function(model_family, sub, predictors, pooled) {
  X <- design_matrix(sub, predictors, pooled)
  n <- nrow(sub)
  k <- qr(X)$rank
  if (model_family %in% c("poisson", "quasi-poisson")) {
    f <- fit_poisson(X, sub)
    disp <- pearson_dispersion(f, sub$Cases, k)
    if (model_family == "poisson") {
      rows <- summarize_fit(f$params, robust_vcov(f, X, sub, pooled))
      aux <- c(list(dispersion = disp, alpha = NA_real_, overdispersion = disp, cov_type = if (pooled) "cluster" else "HAC"),
               information_criteria(f$llf, k, n))
    } else {
      rows <- summarize_fit(f$params, model_vcov(f, X), scale = disp, df_resid = n - k)
      aux <- c(list(dispersion = disp, alpha = NA_real_, overdispersion = disp, cov_type = "quasi (Pearson scale, t)"),
               information_criteria(f$llf, k, n))
    }
    return(list(rows = rows, aux = aux, llf = f$llf, converged = isTRUE(f$converged)))
  }
  if (model_family == "negbin") {
    nb <- fit_negbin(X, sub)
    f <- nb$m
    aux <- c(list(dispersion = pearson_dispersion(f, sub$Cases, k), alpha = nb$alpha, overdispersion = nb$alpha,
                  cov_type = if (pooled) "cluster" else "HAC"),
             information_criteria(f$llf, k + 1, n))
    return(list(rows = summarize_fit(f$params, robust_vcov(f, X, sub, pooled)), aux = aux, llf = f$llf,
                converged = isTRUE(f$converged) && isTRUE(f$profile_converged)))
  }
  stop(model_family)
}

# Run all three model specifications for the given family
# (model_family in "poisson", "quasi-poisson", "negbin").
run_models <- function(model_family) {
  datasets <- load_pipeline_datasets()
  results <- list()
  for (sp in SPECS) {
    cat(sprintf("\n[%s | %s] %s\n", model_family, sp$spec, sp$formula))
    rows <- list(); diag_rows <- list()
    for (virus in names(datasets)) {
      df <- datasets[[virus]]
      units <- if (sp$pooled) list(ALL = df) else split(df, df$Location_Name)[sort_c(unique(df$Location_Name))]
      for (uf in names(units)) {
        sub <- prepare_glm_frame(units[[uf]], sp$predictors)
        if (nrow(sub) < sp$min_n || sum(sub$Cases) == 0) next
        res <- tryCatch(fit_family(model_family, sub, sp$predictors, sp$pooled),
                        error = function(e) { cat(sprintf("    fit failed for %s %s: %s\n", virus, uf, conditionMessage(e))); NULL })
        if (is.null(res)) next
        common <- c(list(Virus = virus, Location_Name = uf, Model = model_family, Spec = sp$spec, n = nrow(sub)),
                    res$aux[c("dispersion", "alpha", "overdispersion", "cov_type", "AIC")])
        rows[[length(rows) + 1]] <- cbind(res$rows, as.data.table(common))
        diag_rows[[length(diag_rows) + 1]] <- as.data.table(c(common, list(BIC = res$aux$BIC, logLik = res$llf, converged = res$converged)))
      }
    }
    out <- rbindlist(rows)
    if (!sp$pooled && nrow(out)) out <- as.data.table(add_fdr(out, by = c("Virus", "Term")))   # 27 states per (virus, term)
    results[[sp$spec]] <- out
    results[[paste0(sp$spec, "_diag")]] <- rbindlist(diag_rows)
  }
  results
}

# Write three CSVs + one Excel workbook with all spec results.
export_results <- function(model_family, results) {
  base <- file.path(OUTPUT_DIR, paste0("glm_", gsub("-", "_", model_family)))
  for (spec in c("per_state", "lagged", "pooled_fe")) {
    df <- results[[spec]]
    if (!nrow(df)) next
    write_csv_py(df, sprintf("%s_%s.csv", base, spec))
    cat(sprintf("    Saved -> %s_%s.csv (%d rows)\n", basename(base), spec, nrow(df)))
  }
  sheets <- list()
  for (spec in c("per_state", "lagged", "pooled_fe")) {
    if (nrow(results[[spec]])) sheets[[spec]] <- results[[spec]]
    dg <- results[[paste0(spec, "_diag")]]
    if (!is.null(dg) && nrow(dg)) sheets[[paste0(spec, "_diagnostics")]] <- dg
  }
  write_xlsx_sheets(sheets, sprintf("%s_results.xlsx", base))
  cat(sprintf("    Saved -> %s_results.xlsx\n", basename(base)))
}

# Console summary: number of significant Temp / Rain coefficients by spec and virus.
summary_climate_only <- function(results, model_family) {
  cat(sprintf("\n  [%s] summary of climate IRRs (Temp_Avg, Rain_mm):\n", model_family))
  for (spec in c("per_state", "lagged", "pooled_fe")) {
    df <- results[[spec]]
    if (!nrow(df)) next
    clim <- df[Term %in% c("Temp_Avg", "Rain_mm")]
    for (virus in sort_c(unique(clim$Virus))) {
      g <- clim[Virus == virus]
      n_total <- uniqueN(g$Location_Name)
      fdr <- ""
      if ("Significant_FDR" %in% names(g))
        fdr <- sprintf("  | FDR: Temp %d/%d  Rain %d/%d", sum(g$Term == "Temp_Avg" & g$Significant_FDR), n_total,
                       sum(g$Term == "Rain_mm" & g$Significant_FDR), n_total)
      cat(sprintf("    [%-10s] %-12s Temp sig=%d/%d  Rain sig=%d/%d%s\n", spec, virus,
                  sum(g$Term == "Temp_Avg" & g$Significant), n_total, sum(g$Term == "Rain_mm" & g$Significant), n_total, fdr))
    }
  }
}

glm_main <- function(model_family, banner) {
  cat(strrep("=", 62), "\n", sprintf("  %s", banner), "\n", strrep("=", 62), "\n", sep = "")
  results <- run_models(model_family)
  export_results(model_family, results)
  summary_climate_only(results, model_family)
  cat("\nDone.\n\n")
}
