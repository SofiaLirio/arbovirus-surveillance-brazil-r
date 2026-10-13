# Install the CRAN packages this project needs (binary builds on Windows/macOS, so no Rtools).
# Usage (from any folder):  Rscript R/install_packages.R
# Writes R_SESSION_INFO.txt at the project root with the R and package versions used.

PACKAGES <- c(
  "data.table", "readxl", "writexl", "xml2", "stringi", "jsonlite", "digest", "curl",
  "sf", "spdep", "spatialreg", "classInt",
  "sandwich", "lmtest", "MASS", "dlnm", "KFAS",
  "ggplot2", "patchwork", "scales", "ragg", "png"
)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
ROOT <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()

# The system library (C:/Program Files/R/...) is usually read-only: install into the per-user
# library R already looks for (R_LIBS_USER), creating it on first use.
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && file.access(.libPaths()[1], 2) != 0) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}

repos <- getOption("repos")
if (is.null(repos) || !length(repos) || any(repos %in% c("@CRAN@", ""))) options(repos = c(CRAN = "https://cloud.r-project.org"))
pkg_type <- if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") "binary" else getOption("pkgType")

missing <- PACKAGES[!vapply(PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("Installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, lib = .libPaths()[1], type = pkg_type, dependencies = c("Depends", "Imports", "LinkingTo"))
}
still <- PACKAGES[!vapply(PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
if (length(still)) {
  cat("FAILED to install:", paste(still, collapse = ", "), "\n")
  quit(status = 1)
}

info <- c(R.version.string, "",
          sprintf("%-12s %s", PACKAGES, vapply(PACKAGES, function(p) as.character(packageVersion(p)), "")))
writeLines(info, file.path(ROOT, "R_SESSION_INFO.txt"))
cat("All", length(PACKAGES), "packages available.\n")
