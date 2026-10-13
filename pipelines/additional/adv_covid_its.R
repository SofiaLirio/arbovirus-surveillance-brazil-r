# ================================================================================
#  4. COVID-19 INTERRUPTED TIME SERIES
# ================================================================================
#  Notifications fell during the pandemic (care-seeking and surveillance were
#  disrupted), which can bias long-term trends. For every virus x geography the
#  monthly log(1 + incidence) is modelled by segmented regression:
#    y = b0 + b1 t + b2 pandemic + b3 t_since_pandemic + harmonics(K=2) + e
#  with the pandemic window Mar 2020 - Dec 2021 and Newey-West (HAC, 6 lags,
#  Bartlett, no small-sample correction) standard errors with normal p-values.
#  State-level p-values get Benjamini-Hochberg q-values within each virus.
#  A trend-sensitivity table compares the linear slope on t with and without the
#  pandemic months.
#  Outputs (outputs/<date>/additional/):
#    its_pandemic_effects.csv, its_trend_sensitivity.csv, its_national.png
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
quiet_library("sandwich")

PANDEMIC_START <- as.Date("2020-03-01")
PANDEMIC_END <- as.Date("2021-12-01")
HAC_LAGS <- 6

design <- function(s) {
  d <- copy(s)
  d[, y := log1p(as.numeric(Incidence))]
  d[, pandemic := as.numeric(Date >= PANDEMIC_START & Date <= PANDEMIC_END)]
  months_in <- (as.integer(format(d$Date, "%Y")) - 2020L) * 12 + (as.integer(format(d$Date, "%m")) - 3L)
  d[, t_since_pandemic := ifelse(pandemic > 0, months_in, 0)]
  cbind(d, harmonic_terms(d$Month, K = 2))
}

run <- function() {
  datasets <- load_datasets()
  eff_rows <- list(); sens_rows <- list()
  for (it in iter_series(datasets)) {
    d <- design(it$s)
    d <- d[!is.na(y)]
    if (nrow(d) < MIN_OBS || sum(d$pandemic) < 6) next
    m <- tryCatch(stats::lm(y ~ t + pandemic + t_since_pandemic + sin1 + cos1 + sin2 + cos2, data = d), error = function(e) NULL)
    if (is.null(m)) next
    V <- sandwich::NeweyWest(m, lag = HAC_LAGS, prewhite = FALSE, adjust = FALSE)
    b <- stats::coef(m)
    p <- 2 * stats::pnorm(-abs(b / sqrt(diag(V))))
    coef <- b[["pandemic"]]
    eff_rows[[length(eff_rows) + 1]] <- data.table(Virus = it$virus, Level = it$level, Geo = it$geo,
      pandemic_log_effect = coef, pandemic_pct_change = (exp(coef) - 1) * 100, p_value = p[["pandemic"]],
      slope_change_per_month = b[["t_since_pandemic"]], n = nrow(d))
    # trend sensitivity: slope on t with vs without pandemic months
    slope <- function(dd) tryCatch(unname(stats::coef(stats::lm(y ~ t, data = dd))[["t"]]), error = function(e) NA_real_)
    s_all <- slope(d); s_excl <- slope(d[pandemic == 0])
    sens_rows[[length(sens_rows) + 1]] <- data.table(Virus = it$virus, Level = it$level, Geo = it$geo,
      slope_with_pandemic = s_all, slope_without_pandemic = s_excl,
      abs_diff = if (is.finite(s_all) && is.finite(s_excl)) abs(s_all - s_excl) else NA_real_)
  }
  eff <- rbindlist(eff_rows); sens <- rbindlist(sens_rows)
  st <- eff[Level == "state"]
  if (nrow(st)) eff <- rbindlist(list(eff[Level != "state"], as.data.table(add_fdr(st, by = "Virus"))), fill = TRUE)
  save_table(eff, "its_pandemic_effects.csv")
  save_table(sens, "its_trend_sensitivity.csv")

  # Figure: national series per virus with pandemic window
  panels <- lapply(seq_along(datasets), function(r) {
    virus <- names(datasets)[r]
    nat <- monthly_series(datasets[[virus]], "national")$Brazil
    row <- eff[Virus == virus & Level == "national"]
    if (nrow(row)) caption_note("its_national.png", "Pandemic level effect (p)",
                                sprintf("%+.0f%% (%s)", row$pandemic_pct_change[1],
                                        if (row$p_value[1] < 0.001) "p < 0.001" else sprintf("p = %.3f", row$p_value[1])), virus)
    p <- ggplot(nat, aes(Date, Incidence)) +
      annotate("rect", xmin = PANDEMIC_START, xmax = PANDEMIC_END, ymin = -Inf, ymax = Inf, fill = "#e74c3c", alpha = 0.12) +
      geom_rect(data = data.frame(x = PANDEMIC_START), aes(xmin = x, xmax = x, ymin = 0, ymax = 0, fill = "COVID-19 window"),
                inherit.aes = FALSE, alpha = 0.12) +
      geom_line(colour = "#2c3e50", linewidth = 1.4 / 2.13) +
      scale_fill_manual(values = c("COVID-19 window" = "#e74c3c"), name = NULL) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.25))) +
      labs(title = virus, x = if (r == length(datasets)) "Date" else NULL,
           y = "per 100k") + theme_adv()
    if (r == 1) p + theme(legend.position = "inside", legend.position.inside = c(0.99, 0.98), legend.justification = c(1, 1),
                          legend.text = element_text(size = 8), legend.background = element_rect(colour = "grey80", fill = "white"))
    else p + theme(legend.position = "none")
  })
  caption_note("its_national.png", "Former title", "<virus> — national monthly incidence, with pandemic effect")
  save_fig(wrap_plots(panels, ncol = 1), "its_national.png", 10, 2.6 * length(datasets))
  write_caption_notes("adv_covid_its.R")
  cat("  COVID ITS done.\n")
}

if (sys.nframe() == 0L) run()
