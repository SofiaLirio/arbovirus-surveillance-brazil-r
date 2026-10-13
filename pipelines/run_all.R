# Run the full normalized pipeline end-to-end (each stage in its own R process):
#   1. arbovirus_pipeline_en.R   — descriptive + trend + Spearman + Moran
#   2. arbovirus_geospatial.R    — choropleth, LISA, trend maps
#   3. graphpad_pipeline.R       — clinical / demographic plots
#   4. bootstrap_ci.R            — bootstrap 95% CIs alongside normal-approx
#   5. glm_poisson.R             — Poisson regression
#   6. glm_quasi_poisson.R       — quasi-Poisson regression
#   7. glm_negbin.R              — Negative Binomial regression
#   8. plot_glm_results.R        — GLM figures
#   9. merge_figures.R           — composite publication figures
# A failed stage is reported and the run continues; the exit status is 1 if any stage failed.
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

SCRIPTS <- c("arbovirus_pipeline_en.R", "arbovirus_geospatial.R", "graphpad_pipeline.R", "bootstrap_ci.R",
             "glm_poisson.R", "glm_quasi_poisson.R", "glm_negbin.R", "plot_glm_results.R", "merge_figures.R")

failed <- character()
for (s in SCRIPTS) {
  cat(sprintf("\n%s\n>>> Running %s\n%s\n", strrep("=", 70), s, strrep("=", 70)))
  rc <- run_script(path_root("pipelines", s))
  if (rc != 0) {
    cat(sprintf("!!! %s failed (exit %d)\n", s, rc))
    failed <- c(failed, s)
  }
}
cat(sprintf("\n%s\n  ALL DONE — %d/%d succeeded\n", strrep("=", 70), length(SCRIPTS) - length(failed), length(SCRIPTS)))
if (length(failed)) cat("  Failed:", paste(failed, collapse = ", "), "\n")
cat(strrep("=", 70), "\n")
quit(status = if (length(failed)) 1 else 0)
