# ================================================================================
#  6. EPIDEMIC THRESHOLDS — Serfling-type baseline
# ================================================================================
#  For every virus x geography: OLS of log(1 + incidence) on a linear trend and
#  annual harmonics (K=2), fitted on all reported months. Expected = expm1(fit);
#  epidemic threshold = expm1(fit + 1.96 * sd(residuals)) (sd with divisor n);
#  a month is "epidemic" when the observed incidence exceeds the threshold.
#  Outputs (outputs/<date>/additional/):
#    epidemic_threshold_monthly.csv   national + regional monthly series with flags
#    epidemic_threshold_summary.csv   per geography: months flagged and which
#    epidemic_threshold_national.png
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

Z <- 1.96
MIN_OBS_T <- 36

# Series with Expected, Threshold, Epidemic columns (NULL when too short).
serfling <- function(s) {
  d <- cbind(copy(s), harmonic_terms(s$Month, K = 2))
  d[, y := log1p(as.numeric(Incidence))]
  ok <- !is.na(d$y)
  if (sum(ok) < MIN_OBS_T) return(NULL)
  m <- stats::lm(y ~ t + sin1 + cos1 + sin2 + cos2, data = d[ok])
  fitted <- stats::predict(m, newdata = d)
  sd0 <- sqrt(mean(stats::residuals(m)^2))           # numpy std (ddof = 0); residuals have mean 0
  d[, `:=`(Expected = expm1(fitted), Threshold = expm1(fitted + Z * sd0))]
  d[, Epidemic := !is.na(Incidence) & Incidence > Threshold]
  d
}

run <- function() {
  datasets <- load_datasets()
  monthly_rows <- list(); summary_rows <- list()
  for (it in iter_series(datasets)) {
    d <- serfling(it$s)
    if (is.null(d)) next
    n_obs <- sum(!is.na(d$Incidence)); n_ep <- sum(d$Epidemic)
    summary_rows[[length(summary_rows) + 1]] <- data.table(Virus = it$virus, Level = it$level, Geo = it$geo,
      n_months = n_obs, n_epidemic_months = n_ep, pct_epidemic = py_round(100 * n_ep / max(n_obs, 1), 1),
      epidemic_months = paste(format(d$Date[d$Epidemic], "%Y-%m"), collapse = ";"))
    if (it$level %in% c("national", "region"))
      monthly_rows[[length(monthly_rows) + 1]] <- data.table(Virus = it$virus, Level = it$level, Geo = it$geo, Date = d$Date,
        Observed = d$Incidence, Expected = d$Expected, Threshold = d$Threshold, Epidemic = d$Epidemic)
  }
  monthly <- rbindlist(monthly_rows); summary <- rbindlist(summary_rows)
  save_table(monthly, "epidemic_threshold_monthly.csv")
  save_table(summary, "epidemic_threshold_summary.csv")

  # Sanity: did dengue flag the 2024 surge nationally?
  dn <- summary[Virus == "Dengue" & Level == "national"]
  if (nrow(dn)) {
    flagged <- grep("^2024", strsplit(dn$epidemic_months[1], ";", fixed = TRUE)[[1]], value = TRUE)
    cat(sprintf("    Dengue national epidemic months in 2024: [%s]\n", paste(sprintf("'%s'", flagged), collapse = ", ")))
  }

  # Figure: national series + threshold + flags
  panels <- lapply(seq_along(datasets), function(r) {
    virus <- names(datasets)[r]
    d <- monthly[Virus == virus & Level == "national"]
    if (!nrow(d)) return(plot_spacer())
    long <- rbind(data.frame(Date = d$Date, v = d$Observed, s = "Observed"), data.frame(Date = d$Date, v = d$Threshold, s = "Epidemic threshold"))
    p <- ggplot() +
      geom_line(data = long, aes(Date, v, colour = s, linetype = s, linewidth = s)) +
      geom_point(data = d[Epidemic == TRUE], aes(Date, Observed, shape = "Epidemic month"), colour = "#c0392b", size = mpl_s_to_size(14)) +
      scale_colour_manual(values = c(Observed = "#2c3e50", "Epidemic threshold" = "#c0392b"), breaks = c("Observed", "Epidemic threshold"), name = NULL) +
      scale_linetype_manual(values = c(Observed = "solid", "Epidemic threshold" = "dashed"), breaks = c("Observed", "Epidemic threshold"), name = NULL) +
      scale_linewidth_manual(values = c(Observed = 1.3 / 2.13, "Epidemic threshold" = 1 / 2.13), breaks = c("Observed", "Epidemic threshold"), name = NULL) +
      scale_shape_manual(values = c("Epidemic month" = 16), name = NULL) +
      labs(title = virus, x = if (r == length(datasets)) "Date" else NULL,
           y = "per 100k") + theme_adv()
    if (r == 1) p + theme(legend.position = "inside", legend.position.inside = c(0.01, 0.98), legend.justification = c(0, 1),
                          legend.text = element_text(size = 8), legend.background = element_rect(colour = "grey80", fill = "white"),
                          legend.spacing.y = unit(0, "pt"), legend.margin = margin(1, 4, 1, 4))
    else p + theme(legend.position = "none")
  })
  caption_note("epidemic_threshold_national.png", "Former title", "<virus> — national incidence vs epidemic threshold")
  save_fig(wrap_plots(panels, ncol = 1), "epidemic_threshold_national.png", 10, 2.6 * length(datasets))
  write_caption_notes("adv_epidemic_threshold.R")
  cat("  Epidemic threshold done.\n")
}

if (sys.nframe() == 0L) run()
