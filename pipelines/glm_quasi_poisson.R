# ================================================================================
#  QUASI-POISSON REGRESSION
# ================================================================================
#  Same point estimates as Poisson, but with SEs inflated by sqrt(dispersion)
#  where dispersion = Pearson chi-squared / residual df.  This handles
#  over-dispersion without specifying a different distribution
#  (R's family = quasipoisson()).
#  Output:
#    outputs/glm_quasi_poisson_per_state.csv
#    outputs/glm_quasi_poisson_lagged.csv
#    outputs/glm_quasi_poisson_pooled_fe.csv
#    outputs/glm_quasi_poisson_results.xlsx
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
source(path_root("pipelines", "glm_common.R"))

if (sys.nframe() == 0L) glm_main("quasi-poisson", "QUASI-POISSON GLM — Poisson with Pearson dispersion correction")
