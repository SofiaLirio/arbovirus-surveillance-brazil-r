# ================================================================================
#  3. SPATIAL REGRESSION — spatial-lag / spatial-error models
# ================================================================================
#  The main pipeline shows strong spatial autocorrelation (Moran's I) but never
#  models it. Per virus, state-level log(1 + mean annual incidence) is regressed on
#  the state's mean temperature and rainfall with:
#    OLS (+ LM-lag / LM-error diagnostics and residual Moran's I),
#    ML spatial-lag (y = rho W y + X b + e) and ML spatial-error (u = lambda W u + e),
#  using the KNN-5 row-standardised weights of the main pipeline.
#  AIC follows the spreg conventions of the original outputs:
#    OLS 2k + n(log(2 pi RSS/n) + 1) with k = 3; ML_Lag 2*4 - 2 logLik; ML_Error 2*3 - 2 logLik.
#  Outputs (outputs/<date>/additional/):
#    spatial_regression_coefficients.csv, spatial_regression_diagnostics.csv,
#    spatial_regression.png
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
quiet_library("spdep", "spatialreg")

# Per-state mean annual incidence + mean climate, per virus.
state_table <- function(datasets) {
  annual <- compute_annual(datasets)
  lapply(setNames(names(datasets), names(datasets)), function(virus) {
    clim <- datasets[[virus]][, .(Temp_Avg = mean(Temp_Avg, na.rm = TRUE), Rain_mm = mean(Rain_mm, na.rm = TRUE)), by = Location_Name]
    inc <- annual[[virus]][, .(Annual_Inc_100k = mean(Annual_Inc_100k, na.rm = TRUE)), by = Location_Name]
    merge(clim, inc, by = "Location_Name")
  })
}

wald_p <- function(est, se) 2 * stats::pnorm(-abs(est / se))

run <- function() {
  datasets <- load_datasets()
  gdf <- get_brazil_states()
  ord <- gdf$Location_Name
  tables <- state_table(datasets)
  coef_rows <- list(); diag_rows <- list()
  set.seed(SEED)
  for (virus in names(tables)) {
    t <- tables[[virus]][match(ord, Location_Name)]
    t <- t[stats::complete.cases(t)]
    if (nrow(t) < 20) { cat(sprintf("  [%s] too few states (%d), skipped\n", virus, nrow(t))); next }
    g2 <- gdf[gdf$Location_Name %in% t$Location_Name, ]
    t <- t[match(g2$Location_Name, Location_Name)]
    w <- knn_listw(cbind(g2$centroid_x, g2$centroid_y), k = min(5, nrow(g2) - 1))
    lw <- w$spdep
    d <- data.frame(y = log1p(t$Annual_Inc_100k), Temp_Avg = t$Temp_Avg, Rain_mm = t$Rain_mm)
    n <- nrow(d)
    cat(sprintf("  [%s] spatial regression on %d states ...\n", virus, n))
    add_coefs <- function(model, est, p) {
      terms <- c("CONSTANT", "Temp_Avg", "Rain_mm", setdiff(names(est), c("(Intercept)", "Temp_Avg", "Rain_mm")))
      coef_rows[[length(coef_rows) + 1]] <<- data.table(Virus = virus, Model = model, Term = terms,
                                                        Coef = unname(est), p_value = unname(p))
    }
    # OLS + diagnostics
    ols <- stats::lm(y ~ Temp_Avg + Rain_mm, data = d)
    so <- summary(ols)$coefficients
    add_coefs("OLS", setNames(so[, 1], rownames(so)), so[, 4])
    rss <- sum(stats::residuals(ols)^2); k <- 3
    mo <- moran_global(stats::residuals(ols), w, permutations = 999)
    lm_tests <- spdep::lm.RStests(ols, lw, test = c("RSlag", "RSerr"))
    diag_rows[[length(diag_rows) + 1]] <- data.table(Virus = virus, Model = "OLS", n = n, r2 = summary(ols)$r.squared,
      AIC = 2 * k + n * (log(2 * pi * rss / n) + 1), logLik = NA_real_, resid_Moran_I = mo$I, resid_Moran_p = mo$p_sim,
      LM_lag_p = lm_tests$RSlag$p.value[[1]], LM_error_p = lm_tests$RSerr$p.value[[1]], spatial_param = NA_real_)
    # ML spatial lag / spatial error
    for (label in c("ML_Lag", "ML_Error")) {
      m <- tryCatch(if (label == "ML_Lag") spatialreg::lagsarlm(y ~ Temp_Avg + Rain_mm, d, lw, method = "eigen", quiet = TRUE)
                    else spatialreg::errorsarlm(y ~ Temp_Avg + Rain_mm, d, lw, method = "eigen", quiet = TRUE),
                    error = function(e) { cat(sprintf("    %s failed: %s\n", label, conditionMessage(e))); NULL })
      if (is.null(m)) next
      b <- m$coefficients; se <- m$rest.se
      if (label == "ML_Lag") {
        est <- c(b, rho = unname(m$rho)); p <- c(wald_p(b, se), wald_p(m$rho, m$rho.se))
        predy <- drop(cbind(1, d$Temp_Avg, d$Rain_mm) %*% b) + m$rho * spatial_lag(w, d$y)
        k_aic <- 4; sp <- unname(m$rho)
      } else {
        est <- c(b, lambda = unname(m$lambda)); p <- c(wald_p(b, se), wald_p(m$lambda, m$lambda.se))
        predy <- drop(cbind(1, d$Temp_Avg, d$Rain_mm) %*% b)
        k_aic <- 3; sp <- unname(m$lambda)
      }
      add_coefs(label, est, p)
      ll <- as.numeric(m$LL)
      diag_rows[[length(diag_rows) + 1]] <- data.table(Virus = virus, Model = label, n = n,
        r2 = stats::cor(d$y, predy)^2, AIC = 2 * k_aic - 2 * ll, logLik = ll, resid_Moran_I = NA_real_,
        resid_Moran_p = NA_real_, LM_lag_p = NA_real_, LM_error_p = NA_real_, spatial_param = sp)
    }
  }
  coefs <- rbindlist(coef_rows); diags <- rbindlist(diag_rows)
  save_table(coefs, "spatial_regression_coefficients.csv")
  save_table(diags, "spatial_regression_diagnostics.csv")

  # Figure: AIC per model/virus + residual Moran's I from OLS
  viruses <- names(tables); models <- c("OLS", "ML_Lag", "ML_Error")
  a <- diags[, .(AIC = mean(AIC)), by = .(Virus, Model)]
  a[, x := match(Virus, viruses) + (match(Model, models) - 2) * 0.26]
  p1 <- ggplot(a, aes(x, AIC, fill = factor(Model, models))) + geom_col(width = 0.26) +
    scale_fill_manual(values = setNames(TAB10[1:3], models), name = NULL) +
    scale_x_continuous(breaks = seq_along(viruses), labels = viruses) +
    labs(x = NULL, y = "AIC (lower = better)") + theme_adv() +
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1),
          legend.text = element_text(size = 8), legend.background = element_rect(colour = "grey80", fill = "white"))
  ov <- diags[Model == "OLS"]
  p2 <- ggplot(ov, aes(factor(Virus, Virus), resid_Moran_I, fill = resid_Moran_p < 0.05)) + geom_col(width = 0.8) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.8 / 2.13) +
    scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#95a5a6"), guide = "none") +
    labs(x = NULL, y = "Moran's I of OLS residuals") + theme_adv()
  caption_note("spatial_regression.png", "Former title", "Model fit by virus", "left")
  caption_note("spatial_regression.png", "Former title", "Residual spatial autocorrelation", "right")
  caption_note("spatial_regression.png", "Key (former title)", "Red bar = residual Moran's I with p < 0.05", "right")
  save_fig(p1 + p2, "spatial_regression.png", 12, 4.2)
  write_caption_notes("adv_spatial_regression.R")
  cat("  Spatial regression done.\n")
}

if (sys.nframe() == 0L) run()
