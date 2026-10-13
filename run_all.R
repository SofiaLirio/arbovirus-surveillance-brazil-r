# Full rebuild from raw data: packages -> raw->db ETL -> all analyses.
# Usage (from any folder):  Rscript run_all.R [--skip-etl] [--skip-install]
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
skip_install <- "--skip-install" %in% commandArgs(TRUE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!skip_install) {
  rc <- system2(rscript, shQuote(file.path(.root, "R", "install_packages.R")))
  if (rc != 0) stop("package installation failed")
}
source(file.path(.root, "R", "utils.R"))

steps <- c(if (!cli_flag("skip-etl")) "etl/build_db.R", "pipelines/run_all.R",
           "pipelines/additional/run_additional_analyses.R", "pipelines/alt_results.R", "pipelines/build_source_data.R")
for (s in steps) {
  cat(sprintf("\n######## %s ########\n", s))
  rc <- run_script(path_root(s))
  if (rc != 0) { cat(sprintf("%s failed (exit %d)\n", s, rc)); quit(status = 1) }
}
cat("\nALL STEPS FINISHED. Results: outputs/, figures/, source_data/, clinical_tables/\n")
