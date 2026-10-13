# ================================================================================
#  PLOT GLM RESULTS — Poisson / quasi-Poisson / Negative Binomial
# ================================================================================
#  Reads the CSVs produced by glm_common.R (run_models / export_results) and
#  generates publication-style figures of the climate IRRs.
#  Outputs (PNG, 300 dpi) -> figures/:
#    glm_fig1_pooled_climate_irr.png   : pooled-FE IRR bars (Temp, Rain) by virus and model family
#    glm_fig2_per_state_heatmap.png    : per-state IRR heatmap (negbin), Temp/Rain
#    glm_fig3_lag_profile.png          : median IRR by lag (0/1/2) for Temp & Rain
#    glm_fig4_aic_comparison.png       : model-fit AIC across family x spec x virus
#    glm_fig5_forest_significant.png   : forest plot of significant per-state IRRs
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
source(path_root("pipelines", "data_loader.R"))
quiet_library("ggplot2", "patchwork")

VIRUS_ORDER  <- c("Dengue", "Chikungunya", "Zika")
VIRUS_COLOR  <- c(Dengue = "#d7263d", Chikungunya = "#1b998b", Zika = "#2e86ab")
FAMILY_ORDER <- c("poisson", "quasi_poisson", "negbin")
FAMILY_LABEL <- c(poisson = "Poisson", quasi_poisson = "Quasi-Poisson", negbin = "Negative Binomial")
FAMILY_COLOR <- setNames(mpl_cmap("viridis")(0.2 + 0.3 * 0:2), FAMILY_LABEL)
LAG_TERMS_TEMP <- c("Temp_Avg", "Temp_Lag1", "Temp_Lag2")
LAG_TERMS_RAIN <- c("Rain_mm", "Rain_Lag1", "Rain_Lag2")
LAG_LABELS <- c("Lag 0", "Lag 1", "Lag 2")
BAR_W <- 0.26

theme_glm <- function() {
  theme_mpl(9, grid = FALSE) + theme(plot.title = element_text(size = 11), axis.title = element_text(size = 10),
                                     legend.text = element_text(size = 9))
}
read_glm <- function(family, spec) fread(file.path(OUTPUT_DIR, sprintf("glm_%s_%s.csv", family, spec)))
save_fig <- function(p, name, w, h) {
  out <- file.path(FIGURES_DIR, name)
  save_png(p, out, w, h, dpi = 300)
  cat(sprintf("  saved -> %s\n", out))
}

fig1_pooled_climate_irr <- function() {
  data <- rbindlist(lapply(FAMILY_ORDER, function(fam) read_glm(fam, "pooled_fe")[Term %in% c("Temp_Avg", "Rain_mm")][, Family := fam]))
  data[, `:=`(Fam = factor(FAMILY_LABEL[Family], FAMILY_LABEL), Virus = factor(Virus, VIRUS_ORDER))]
  ylabs <- c(Temp_Avg = "IRR per +1 °C (mean temperature)", Rain_mm = "IRR per +1 mm (rainfall)")
  panels <- lapply(names(ylabs), function(term) {
    sub <- data[Term == term]
    sub[, pad := 0.01 * max(IRR_CI_U, na.rm = TRUE), by = Family]
    sub[, x := as.integer(Virus) + (match(Family, FAMILY_ORDER) - 2) * BAR_W]   # matplotlib-style offsets
    p <- ggplot(sub, aes(x, IRR, fill = Fam)) +
      geom_col(width = BAR_W, colour = "black", linewidth = 0.6 / 2.13) +
      geom_errorbar(aes(ymin = IRR_CI_L, ymax = IRR_CI_U), width = 0.08, linewidth = 0.8 / 2.13) +
      geom_text(data = sub[Significant == TRUE], aes(y = IRR_CI_U + pad, label = "*"),
                vjust = 0, size = 12 / .pt, fontface = "bold") +
      geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.8 / 2.13, linetype = "dashed") +
      scale_x_continuous(breaks = seq_along(VIRUS_ORDER), labels = VIRUS_ORDER) +
      scale_fill_manual(values = FAMILY_COLOR, name = NULL) +
      labs(x = NULL, y = ylabs[[term]]) + theme_glm()
    if (term == "Temp_Avg") p + theme(legend.position = "inside", legend.position.inside = c(0.5, 1), legend.justification = c(0.5, 1),
                                      legend.background = element_blank(), legend.key.height = unit(10, "pt"))
    else p + theme(legend.position = "none")
  })
  caption_note("glm_fig1_pooled_climate_irr.png", "Former title", "Pooled GLM with state + month fixed effects — climate IRRs")
  caption_note("glm_fig1_pooled_climate_irr.png", "Key (former title)", "95% CI; cluster-robust SEs for Poisson / NegBin; * p < 0.05")
  save_fig(wrap_plots(panels, nrow = 1), "glm_fig1_pooled_climate_irr.png", 12, 5)
}

fig2_per_state_heatmap <- function() {
  df <- read_glm("negbin", "per_state")
  xlabs <- c(Temp_Avg = "Mean temperature (IRR per +1 °C)", Rain_mm = "Rainfall (IRR per +1 mm)")
  panels <- lapply(names(xlabs), function(term) {
    sub <- df[Term == term, .(Location_Name, Virus, IRR, Significant)]
    sub[, lg := log2(IRR)]
    vmax <- max(abs(sub$lg), na.rm = TRUE)
    sub[, `:=`(Location_Name = factor(Location_Name, rev(sort_c(unique(Location_Name)))), Virus = factor(Virus, VIRUS_ORDER),
               lab = sprintf("%.2f%s", IRR, ifelse(Significant, "*", "")))]
    ggplot(sub, aes(Virus, Location_Name, fill = lg)) + geom_tile() +
      geom_text(aes(label = lab, colour = abs(lg) > vmax * 0.55), size = 7 / .pt, family = PLOT_FONT) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
      scale_fill_gradientn(colours = mpl_colours("RdBu_r"), limits = c(-vmax, vmax), name = "log₂(IRR)",
                           guide = guide_colourbar(barheight = unit(6.5, "in"))) +
      scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
      labs(x = xlabs[[term]], y = NULL) + theme_glm() + theme(axis.line = element_blank())
  })
  caption_note("glm_fig2_per_state_heatmap.png", "Former title", "Per-state Negative-Binomial GLM — climate IRRs by state and virus")
  caption_note("glm_fig2_per_state_heatmap.png", "Key (former title)", "Month FE, HAC SEs; * significant at p < 0.05; red = positive, blue = negative")
  save_fig(wrap_plots(panels, nrow = 1), "glm_fig2_per_state_heatmap.png", 14, 9)
}

fig3_lag_profile <- function() {
  df <- read_glm("negbin", "lagged")
  panels <- list(list(LAG_TERMS_TEMP, "IRR per +1 °C (mean temperature)"),
                 list(LAG_TERMS_RAIN, "IRR per +1 mm (rainfall)"))
  plots <- lapply(panels, function(pn) {
    d <- rbindlist(lapply(VIRUS_ORDER, function(v) rbindlist(lapply(seq_along(pn[[1]]), function(i) {
      x <- df[Virus == v & Term == pn[[1]][i]]$IRR; x <- x[is.finite(x)]
      data.table(Virus = v, Lag = i - 1, med = if (length(x)) median(x) else NA_real_,
                 lo = if (length(x)) quantile(x, 0.25, names = FALSE) else NA_real_,
                 hi = if (length(x)) quantile(x, 0.75, names = FALSE) else NA_real_)
    }))))
    d[, Virus := factor(Virus, VIRUS_ORDER)]
    ggplot(d, aes(Lag, med, colour = Virus, fill = Virus)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
      geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.8 / 2.13, linetype = "dashed") +
      geom_line(linewidth = 2 / 2.13) + geom_point(size = 7 / .pt * 1.3) +
      scale_colour_manual(values = VIRUS_COLOR, name = NULL) + scale_fill_manual(values = VIRUS_COLOR, name = NULL) +
      scale_x_continuous(breaks = 0:2, labels = LAG_LABELS) +
      labs(x = "Lag (months)", y = pn[[2]]) + theme_glm() +
      theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1))
  })
  caption_note("glm_fig3_lag_profile.png", "Former title", "Per-state NegBin GLM with lagged climate predictors and month FE — median IRR across states")
  caption_note("glm_fig3_lag_profile.png", "Key (former title)", "Line = median IRR across states; shaded = IQR")
  save_fig(wrap_plots(plots, nrow = 1), "glm_fig3_lag_profile.png", 12, 5)
}

fig4_aic_comparison <- function() {
  specs <- c("per_state", "lagged", "pooled_fe")
  agg <- rbindlist(lapply(FAMILY_ORDER, function(fam) rbindlist(lapply(specs, function(spec) {
    u <- unique(read_glm(fam, spec)[, .(Virus, Location_Name, AIC)])[is.finite(AIC)]
    u[, .(AIC_median = median(AIC), AIC_n = .N), by = Virus][, `:=`(Family = fam, Spec = spec)]
  }))))
  agg[, `:=`(Fam = factor(FAMILY_LABEL[Family], FAMILY_LABEL),
             SpecLab = factor(c(per_state = "per-state", lagged = "lagged", pooled_fe = "pooled FE")[Spec],
                              c("per-state", "lagged", "pooled FE")))]
  agg[, x := as.integer(SpecLab) + (match(Family, FAMILY_ORDER) - 2) * BAR_W]
  panels <- lapply(VIRUS_ORDER, function(v) {
    ggplot(agg[Virus == v], aes(x, AIC_median, fill = Fam)) +
      geom_col(width = BAR_W, colour = "black", linewidth = 0.6 / 2.13) +
      scale_x_continuous(breaks = 1:3, labels = levels(agg$SpecLab)) +
      scale_y_log10() + scale_fill_manual(values = FAMILY_COLOR, name = NULL) +
      labs(title = v, x = NULL, y = "Median AIC  (log scale; lower = better fit)") + theme_glm() +
      theme(plot.title = element_text(colour = VIRUS_COLOR[[v]]), legend.text = element_text(size = 8),
            legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1))
  })
  caption_note("glm_fig4_aic_comparison.png", "Former title", "Model fit comparison — median AIC across states (log scale)")
  caption_note("glm_fig4_aic_comparison.png", "Key (former title)",
               "Negative Binomial fits dengue/chikv/zika monthly counts vastly better than Poisson due to overdispersion")
  save_fig(wrap_plots(panels, nrow = 1), "glm_fig4_aic_comparison.png", 14, 4.5)
}

fig5_forest_significant <- function() {
  df <- read_glm("negbin", "per_state")[Term %in% c("Temp_Avg", "Rain_mm") & Significant == TRUE]
  if (!nrow(df)) { cat("  no significant per-state rows; skipping forest plot\n"); return(invisible()) }
  df[, `:=`(log_irr = log2(IRR), log_lo = log2(IRR_CI_L), log_hi = log2(IRR_CI_U))]
  xlabs <- c(Temp_Avg = "log₂(IRR), mean temperature   (right = positive,  left = negative)",
             Rain_mm = "log₂(IRR), rainfall   (right = positive,  left = negative)")
  panels <- lapply(names(xlabs), function(term) {
    sub <- df[Term == term]
    if (!nrow(sub)) return(ggplot() + annotate("text", x = 0, y = 0, size = 9 / .pt, family = PLOT_FONT,
                                               label = sprintf("No significant %s association", if (term == "Temp_Avg") "temperature" else "rainfall")) +
                             theme_void())
    sub[, label := sprintf("%s  (%s)", Location_Name, Virus)]
    sub <- sub[order_c(Virus, log_irr)]
    sub[, y := seq_len(.N)]
    ggplot(sub, aes(log_irr, y)) +
      geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.8 / 2.13, linetype = "dashed") +
      geom_errorbar(aes(xmin = log_lo, xmax = log_hi), orientation = "y", colour = "grey50", width = 0.3, linewidth = 0.8 / 2.13) +
      geom_point(aes(fill = factor(Virus, VIRUS_ORDER)), shape = 21, colour = "black", size = mpl_s_to_size(30), stroke = 0.4 / 2.13) +
      scale_fill_manual(values = VIRUS_COLOR, name = NULL, drop = FALSE) +
      scale_y_continuous(breaks = sub$y, labels = sub$label, expand = expansion(add = 0.6)) +
      labs(x = xlabs[[term]], y = NULL) + theme_glm() +
      theme(axis.text.y = element_text(size = 7))
  })
  p <- wrap_plots(panels, nrow = 1, guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 9))
  p <- p & guides(fill = guide_legend(override.aes = list(size = 3)))
  caption_note("glm_fig5_forest_significant.png", "Former title",
               "Forest plot of significant per-state climate effects (Negative-Binomial GLM, p < 0.05)")
  save_fig(p, "glm_fig5_forest_significant.png", 13, 11)
}

main <- function() {
  cat("Plotting GLM results ...\n")
  fig1_pooled_climate_irr()
  fig2_per_state_heatmap()
  fig3_lag_profile()
  fig4_aic_comparison()
  fig5_forest_significant()
  write_caption_notes("plot_glm_results.R")
  cat("Done.\n")
}

if (sys.nframe() == 0L) main()
