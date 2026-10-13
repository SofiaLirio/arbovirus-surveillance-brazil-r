# Build the analysis inputs from raw data:
#   reference/ -> db/ (dim_location, brazil_states.geojson, graphpad_data/)
#   raw/ibge   -> db/normalized_fact_population.csv
#   raw/sinan  -> db/normalized_fact_cases.csv + clinical_tables/
#   raw/inmet  -> db/normalized_fact_climate.csv
# then runs check_db.R. Run from anywhere: Rscript etl/build_db.R
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
source(path_root("etl", "_paths.R"))

file.copy(file.path(REF, "dim_location.csv"), file.path(DB, "normalized_dim_location.csv"), overwrite = TRUE, copy.date = TRUE)
file.copy(file.path(REF, "brazil_states.geojson"), file.path(DB, "brazil_states.geojson"), overwrite = TRUE, copy.date = TRUE)
dir.create(file.path(DB, "graphpad_data"), showWarnings = FALSE)
file.copy(list.files(file.path(REF, "graphpad_data"), full.names = TRUE), file.path(DB, "graphpad_data"),
          overwrite = TRUE, copy.date = TRUE)
cat("reference files copied into db/\n")

for (script in c("build_population.R", "build_cases_clinical.R", "build_climate.R", "check_db.R")) {
  cat(sprintf("\n%s\n>>> %s\n%s\n", strrep("=", 70), script, strrep("=", 70)))
  rc <- run_script(path_root("etl", script))
  if (rc != 0) {
    cat(sprintf("!!! %s failed (exit %d)\n", script, rc))
    quit(status = 1)
  }
}
cat("\nDB BUILD DONE\n")
