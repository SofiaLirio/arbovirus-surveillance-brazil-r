# ================================================================================
#  NEGATIVE BINOMIAL REGRESSION (NB2 parameterization)
# ================================================================================
#  Cases ~ Climate + offset(log Population), family = NegativeBinomial.
#  Handles over-dispersion via a quadratic mean-variance relationship:
#    Var(Y) = mu + alpha * mu^2
#  where alpha > 0 is the dispersion parameter estimated by maximum likelihood.
#  Output:
#    outputs/glm_negbin_per_state.csv
#    outputs/glm_negbin_lagged.csv
#    outputs/glm_negbin_pooled_fe.csv
#    outputs/glm_negbin_results.xlsx
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

if (sys.nframe() == 0L) glm_main("negbin", "NEGATIVE BINOMIAL GLM — Cases ~ Climate + offset(log Pop)")
