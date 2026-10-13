# ================================================================================
#  1. DLNM — non-linear distributed-lag climate model
# ================================================================================
#  Gap addressed: the project's GLMs enter climate (and its lags) as LINEAR terms,
#  but temperature-arbovirus relationships are non-linear (an Aedes transmission
#  optimum near 26-30 C). This builds a distributed-lag non-linear crossbasis
#  (natural cubic spline on the exposure x a linear lag basis over lags 0-2 - the
#  only lags present in the data) and fits one pooled Negative-Binomial GLM per
#  virus with month + state fixed effects and a log-population offset.
#  The crossbasis is dlnm::crossbasis(argvar = ns(knots), arglag = "lin") written
#  out with splines::ns: knots at the 33.3/66.7 % (temperature) or 50 % (rain)
#  quantiles of the distinct pooled lag 0-2 exposure values, boundary knots at its range.
#  Outputs (outputs/<date>/additional/):
#    dlnm_exposure_response.csv   cumulative IRR curves vs Temp / Rain, with 95% CI
#    dlnm_summary.csv             per-virus temperature optimum + overdispersion
#    dlnm_exposure_response.png   figure (3 viruses x {Temperature, Rainfall})
#  Limitation (written into dlnm_summary.csv): only lags 0-2 are available, so the
#  lag dimension of the crossbasis is 3 with a linear lag basis.
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
source(path_root("pipelines", "glm_common.R"))

LAGS <- 0:2
LAG_BASIS <- cbind(1, LAGS)                       # [1, lag]
PREDICTORS <- list(Temp = list(cols = c("Temp_Avg", "Temp_Lag1", "Temp_Lag2"), df = 4, unit = "°C"),
                   Rain = list(cols = c("Rain_mm", "Rain_Lag1", "Rain_Lag2"), df = 3, unit = "mm"))
GRID_N <- 60

# Natural cubic spline exposure basis with df-2 inner knots at evenly spaced quantiles of the
# DISTINCT pooled exposure values (patsy cr(x, df) knot placement); returns x -> basis matrix.
var_basis <- function(values, df) {
  n_inner <- df - 2
  probs <- seq(0, 1, length.out = n_inner + 2)[-c(1, n_inner + 2)]
  knots <- stats::quantile(unique(values), probs, type = 7, names = FALSE)
  bounds <- range(values)
  function(x) splines::ns(x, knots = knots, Boundary.knots = bounds, intercept = FALSE)
}

# n x (p*q) crossbasis: var-basis(exposure at each lag) combined with the linear lag basis
# and summed over lags.
crossbasis <- function(sub, cols, vb, prefix) {
  B <- lapply(cols, function(cc) unclass(vb(sub[[cc]])))
  p <- ncol(B[[1]]); q <- ncol(LAG_BASIS)
  cb <- matrix(0, nrow(sub), p * q)
  nm <- character(p * q)
  for (i in seq_len(p)) for (j in seq_len(q)) {
    k <- (i - 1) * q + j
    cb[, k] <- Reduce(`+`, lapply(seq_along(LAGS), function(l) B[[l]][, i] * LAG_BASIS[l, j]))
    nm[k] <- sprintf("%s_cb%d_%d", prefix, i - 1, j - 1)
  }
  colnames(cb) <- nm
  cb
}

# Cumulative design rows (grid minus reference) for a constant exposure across lags.
crosspred_rows <- function(vb, grid, ref) {
  colsum <- colSums(LAG_BASIS)
  row <- function(x) { bx <- unclass(vb(x))[1, ]; as.vector(t(outer(bx, colsum))) }
  r0 <- row(ref)
  t(vapply(grid, function(x) row(x) - r0, numeric(length(r0))))
}

run <- function() {
  datasets <- load_datasets()
  curves <- list(); summary <- list()
  for (virus in names(datasets)) {
    df <- datasets[[virus]]
    cols_needed <- c("Cases", "Population", "Location_Name", "Month", PREDICTORS$Temp$cols, PREDICTORS$Rain$cols)
    sub <- as.data.frame(df)[, cols_needed]
    sub <- sub[stats::complete.cases(sub) & sub$Population > 0, ]
    if (nrow(sub) < 100) { cat(sprintf("  [%s] too few rows, skipped\n", virus)); next }
    sub$log_pop <- log(as.numeric(sub$Population))
    vbs <- list(); cbs <- list()
    for (pred in names(PREDICTORS)) {
      cfg <- PREDICTORS[[pred]]
      vbs[[pred]] <- var_basis(unlist(sub[cfg$cols], use.names = FALSE), cfg$df)
      cbs[[pred]] <- crossbasis(sub, cfg$cols, vbs[[pred]], tolower(pred))
    }
    X <- cbind(design_matrix(sub, character(), pooled = TRUE), cbs$Temp, cbs$Rain)
    cat(sprintf("  [%s] fitting pooled NegBin DLNM (n=%d) ...\n", virus, nrow(sub)))
    nb <- tryCatch(fit_negbin(X, sub), error = function(e) { cat("    FAILED:", conditionMessage(e), "\n"); NULL })
    if (is.null(nb)) next
    f <- nb$m
    # Covariance = inverse observed information at the estimates (what the published curves used:
    # statsmodels fitted this rank-deficient patsy design by Newton, whose covariance is -H^-1).
    beta <- f$params
    V <- solve(crossprod(X, X * oim_factor(sub$Cases, f$mu, nb$alpha)))
    dimnames(V) <- list(colnames(X), colnames(X))
    for (pred in names(PREDICTORS)) {
      cfg <- PREDICTORS[[pred]]
      terms <- colnames(cbs[[pred]])
      x0 <- sub[[cfg$cols[1]]]
      grid <- seq(stats::quantile(x0, 0.01, type = 7, names = FALSE), stats::quantile(x0, 0.99, type = 7, names = FALSE),
                  length.out = GRID_N)
      ref <- stats::median(x0)
      D <- crosspred_rows(vbs[[pred]], grid, ref)
      logrr <- drop(D %*% beta[terms])
      se <- sqrt(pmax(rowSums((D %*% V[terms, terms]) * D), 0))
      irr <- exp(logrr)
      curves[[length(curves) + 1]] <- data.table(Virus = virus, Predictor = pred, Value = grid, IRR = irr,
        IRR_CI_L = exp(logrr - 1.96 * se), IRR_CI_U = exp(logrr + 1.96 * se), Reference = ref, Unit = cfg$unit)
      if (pred == "Temp")
        summary[[length(summary) + 1]] <- data.table(Virus = virus, Temp_optimum_C = round(grid[which.max(irr)], 2),
          IRR_at_optimum = round(max(irr), 3), overdispersion_alpha = round(nb$alpha, 4),
          dispersion = round(pearson_dispersion(f, sub$Cases, qr(X)$rank), 3), n = nrow(sub),
          note = "lags 0-2 only; linear lag basis")
    }
  }
  curves <- rbindlist(curves); summary <- rbindlist(summary)
  save_table(curves, "dlnm_exposure_response.csv")
  save_table(summary, "dlnm_summary.csv")

  # Figure
  panels <- list()
  for (virus in names(datasets)) for (pred in c("Temp", "Rain")) {
    d <- curves[Virus == virus & Predictor == pred]
    lbl <- if (pred == "Temp") "Temperature" else "Rainfall"
    panels[[length(panels) + 1]] <- if (!nrow(d)) plot_spacer() else
      ggplot(d, aes(Value, IRR)) +
        geom_ribbon(aes(ymin = IRR_CI_L, ymax = IRR_CI_U), fill = "#c0392b", alpha = 0.15) +
        geom_hline(yintercept = 1, colour = "grey50", linetype = "dashed", linewidth = 0.8 / 2.13) +
        geom_line(colour = "#c0392b", linewidth = 2 / 2.13) +
        labs(title = virus, x = sprintf("%s (%s)", lbl, d$Unit[1]), y = "Cumulative IRR") +
        theme_adv()
  }
  caption_note("dlnm_exposure_response.png", "Former title", "DLNM cumulative exposure–response (lags 0–2, state+month FE)")
  save_fig(wrap_plots(panels, ncol = 2), "dlnm_exposure_response.png", 11, 3.2 * length(datasets))
  write_caption_notes("adv_dlnm.R")
  cat("  DLNM done.\n")
}

if (sys.nframe() == 0L) run()
