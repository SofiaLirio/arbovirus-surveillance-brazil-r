# Sanity checks on the built db/ + clinical_tables/, and a sha256 comparison with
# reference/expected_hashes.csv (the files the published 2016-2026 results were built from).
# The CSVs must match exactly; xlsx files embed a write timestamp, so a mismatch there is reported
# but not treated as a failure.
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
quiet_library("data.table")

fc   <- fread(file.path(DB, "normalized_fact_cases.csv"))
pop  <- fread(file.path(DB, "normalized_fact_population.csv"))
clim <- fread(file.path(DB, "normalized_fact_climate.csv"))
cl   <- fread(file.path(CLIN, "clinical_demographic_long.csv"))

cat("=== FILE SHAPES ===\n")
for (nm in c("fact_cases", "fact_population", "fact_climate", "clinical_long")) {
  d <- switch(nm, fact_cases = fc, fact_population = pop, fact_climate = clim, clinical_long = cl)
  cat(sprintf("  %-16s rows=%6d years=%d-%d\n", nm, nrow(d), min(d$Year), max(d$Year)))
}

sex <- cl[Dimension == "sex", .(sex_total = sum(Count)), by = .(Disease, Year)]
cas <- fc[, .(cases = sum(Cases)), by = .(Disease, Year)]
rec <- merge(cas, sex, by = c("Disease", "Year"), all = TRUE)
cat("\nclinical 'sex' totals reconcile with fact_cases:",
    if (isTRUE(all(abs(rec$cases - rec$sex_total) <= 1e-8 + 1e-5 * abs(rec$sex_total)))) "True" else "False", "\n")
cat("states with 2026 climate:", uniqueN(clim[Year == 2026]$Location_ID), "/ 27\n")

exp_path <- file.path(REF, "expected_hashes.csv")
if (!file.exists(exp_path)) {
  cat("\nno reference/expected_hashes.csv: reproduction check skipped\n")
  quit(status = 0)
}
expected <- fread(exp_path)
cat("\n=== reproduction check vs reference/expected_hashes.csv ===\n")
bad <- 0
for (i in seq_len(nrow(expected))) {
  rel <- expected$file[i]; want <- expected$sha256[i]
  p <- do.call(path_root, as.list(strsplit(rel, "/", fixed = TRUE)[[1]]))
  got <- if (file.exists(p)) digest::digest(file = p, algo = "sha256") else "missing"
  ok <- identical(got, want)
  status <- if (ok) "OK" else if (grepl("\\.xlsx$", rel) && got != "missing") "differs (xlsx timestamp, ignored)" else "MISMATCH"
  if (!ok && status == "MISMATCH") bad <- bad + 1
  cat(sprintf("  %-34s %s\n", status, rel))
}
quit(status = if (bad) 1 else 0)
