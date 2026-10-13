# ================================================================================
#  BOOTSTRAP 95% CONFIDENCE INTERVALS
# ================================================================================
#  Adds non-parametric bootstrap CIs (1000 resamples, percentile method) alongside
#  the normal-approximation CIs.  More robust to the right-skew typical
#  of count and incidence data.
#  Output files (in ./outputs/):
#    bootstrap_ci_95_cases.csv
#    bootstrap_ci_95_incidence.csv
#    ci_comparison.xlsx     (joins normal-approx and bootstrap CIs in one table)
#  Resampling uses R's RNG (seed 42), so the bootstrap bounds differ from the
#  numpy-based Python version by Monte-Carlo error only.
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

N_BOOT <- 1000
RNG_SEED <- 42

# Percentile bootstrap CI for the mean (the generator is re-seeded on every call).
bootstrap_ci <- function(values, n_boot = N_BOOT, ci = 95, seed = RNG_SEED) {
  set.seed(seed)
  values <- values[is.finite(values)]
  n <- length(values)
  if (n < 2) return(c(NA_real_, NA_real_, n))
  boot_means <- vapply(seq_len(n_boot), function(i) mean(values[sample.int(n, n, replace = TRUE)]), numeric(1))
  q <- stats::quantile(boot_means, c((100 - ci) / 2, 100 - (100 - ci) / 2) / 100, type = 7, names = FALSE)
  c(q, n)
}

ci_row <- function(code, virus, x) {
  x <- x[!is.na(x)]
  m <- mean(x); se <- stats::sd(x) / sqrt(length(x))
  b <- bootstrap_ci(x)
  data.table(Disease = code, Disease_Name = virus, Mean = m, Normal_CI_Lower = m - 1.96 * se, Normal_CI_Upper = m + 1.96 * se,
             Bootstrap_CI_Lower = b[1], Bootstrap_CI_Upper = b[2], n = length(x), method = "bootstrap-percentile-1000")
}

main <- function() {
  cat(strrep("=", 62), "\n  BOOTSTRAP 95% CIs (1000 resamples, percentile method)\n", strrep("=", 62), "\n", sep = "")
  datasets <- load_pipeline_datasets()
  virus_to_code <- setNames(names(DISEASE_KEY), DISEASE_KEY)
  rows_cases <- list(); rows_inc <- list()
  for (virus in names(datasets)) {
    df <- datasets[[virus]]
    rc <- ci_row(virus_to_code[[virus]], virus, df$Cases)
    rows_cases[[virus]] <- rc
    rows_inc[[virus]] <- ci_row(virus_to_code[[virus]], virus, df$Monthly_Inc_per100k)
    cat(sprintf("  [%-12s] cases mean=%.2f normal=(%.2f,%.2f) bootstrap=(%.2f,%.2f)\n", virus, rc$Mean,
                rc$Normal_CI_Lower, rc$Normal_CI_Upper, rc$Bootstrap_CI_Lower, rc$Bootstrap_CI_Upper))
  }
  df_cases <- rbindlist(rows_cases); df_inc <- rbindlist(rows_inc)
  write_csv_py(df_cases, file.path(OUTPUT_DIR, "bootstrap_ci_95_cases.csv"))
  write_csv_py(df_inc, file.path(OUTPUT_DIR, "bootstrap_ci_95_incidence.csv"))
  cat("\n  Saved -> bootstrap_ci_95_cases.csv\n  Saved -> bootstrap_ci_95_incidence.csv\n")
  write_xlsx_sheets(list(Cases = df_cases, Incidence_per100k = df_inc), file.path(OUTPUT_DIR, "ci_comparison.xlsx"))
  cat("  Saved -> ci_comparison.xlsx\n\nDone.\n\n")
}

if (sys.nframe() == 0L) main()
