# ================================================================================
#  ADDITIONAL-ANALYSES ORCHESTRATOR
# ================================================================================
#  Runs the six additional analysis modules that extend the main pipeline:
#    1. adv_dlnm.R                non-linear distributed-lag climate model (DLNM)
#    2. adv_cocirculation.R       cross-arbovirus co-circulation / synchrony
#    3. adv_spatial_regression.R  spatial-lag / spatial-error regression
#    4. adv_covid_its.R           COVID-19 interrupted time series
#    5. adv_forecast.R            out-of-sample SARIMA backtest
#    6. adv_epidemic_threshold.R  epidemic thresholds (Serfling-type baseline)
#  All outputs land in outputs/<RUN_DATE>/additional/ (see adv_common.R OUT_ROOT).
#  Each stage runs in its own R process; a failure is logged and the run
#  continues, exiting 1 at the end if any stage failed. Mirrors pipelines/run_all.R.
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

SCRIPTS <- c("adv_dlnm.R", "adv_cocirculation.R", "adv_spatial_regression.R", "adv_covid_its.R",
             "adv_forecast.R", "adv_epidemic_threshold.R")

ok <- setNames(logical(length(SCRIPTS)), SCRIPTS)
for (s in SCRIPTS) {
  cat(sprintf("\n%s\n>>> Running %s\n%s\n", strrep("=", 70), s, strrep("=", 70)))
  rc <- run_script(path_root("pipelines", "additional", s))
  ok[[s]] <- rc == 0
  if (rc != 0) cat(sprintf("    !!! %s FAILED (exit %d)\n", s, rc))
}
cat(sprintf("\n%s\n  ADDITIONAL ANALYSES DONE - %d/%d succeeded\n", strrep("=", 70), sum(ok), length(SCRIPTS)))
for (s in SCRIPTS) cat(sprintf("    %s  %s\n", if (ok[[s]]) "OK  " else "FAIL", s))
cat(strrep("=", 70), "\n")
quit(status = if (all(ok)) 0 else 1)
