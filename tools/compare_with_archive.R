# Compare freshly downloaded raw sources against the archived copies the paper was
# built from (data/updated_data.rar, data/metereologicos.rar) by size and CRC32,
# without extracting the archives (uses WinRAR's `rar lt`). Also integrity-tests every
# downloaded zip by extracting it to a temporary folder.
#
# Usage:
#   Rscript tools/compare_with_archive.R --data DIR [--raw DIR]
# Writes <raw>/archive_comparison.csv and prints a summary.
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
quiet_library("data.table")

RAR <- local({ p <- Sys.which("rar"); if (nzchar(p)) p else "C:/Program Files/WinRAR/Rar.exe" })
ARCHIVES <- list("updated_data.rar" = c("sinan", "ibge"), "metereologicos.rar" = "inmet")
# SINAN files still receiving notifications; differences are expected
PRELIMINARY <- c("25", "26")

# {name: list(size, crc32)} from `rar lt` (works with localized RAR output)
rar_listing <- function(path) {
  out <- iconv(system2(RAR, c("lt", shQuote(path)), stdout = TRUE, stderr = FALSE), "latin1", "UTF-8")
  entries <- list(); cur <- list()
  for (line in out) {
    if (grepl("^\\s*(Name|Nome):", line)) { cur <- list(name = trimws(sub("^\\s*(Name|Nome):\\s*", "", line))); next }
    if (grepl("^\\s*(Size|Tamanho):\\s*\\d+\\s*$", line) && !is.null(cur$name)) {
      cur$size <- as.numeric(sub("^\\s*(Size|Tamanho):\\s*(\\d+)\\s*$", "\\2", line)); next
    }
    if (grepl("^\\s*CRC32:\\s*[0-9A-Fa-f]{8}", line) && !is.null(cur$name)) {
      crc <- toupper(sub("^\\s*CRC32:\\s*([0-9A-Fa-f]{8}).*$", "\\1", line))
      entries[[basename(gsub("\\\\", "/", cur$name))]] <- list(size = cur$size, crc = crc)
      cur <- list()
    }
  }
  entries
}

crc32_file <- function(path) {
  x <- digest::digest(file = path, algo = "crc32", serialize = FALSE)
  toupper(paste0(strrep("0", max(0, 8 - nchar(x))), x))
}

zip_test <- function(path) {
  tmp <- tempfile("ziptest"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  ok <- tryCatch({ res <- utils::unzip(path, exdir = tmp); length(res) > 0 }, warning = function(w) FALSE, error = function(e) FALSE)
  if (ok) "ok" else "CORRUPT"
}

main <- function() {
  raw <- normalizePath(cli_value("raw", path_root("raw")), winslash = "/", mustWork = FALSE)
  data <- cli_value("data")
  if (is.null(data)) stop("--data DIR is required: the folder holding updated_data.rar and metereologicos.rar (not bundled)")
  rows <- list()
  for (archive in names(ARCHIVES)) {
    apath <- file.path(data, archive)
    if (!file.exists(apath)) stop(sprintf("archive not found: %s", apath))
    listing <- rar_listing(apath)
    files <- unlist(lapply(ARCHIVES[[archive]], function(g) list.files(file.path(raw, g), full.names = TRUE)))
    files <- files[!dir.exists(files)]
    fresh <- setNames(files, basename(files))
    for (name in sort_c(union(names(listing), names(fresh)))) {
      a <- listing[[name]]; p <- if (name %in% names(fresh)) fresh[[name]] else NULL
      f_size <- if (!is.null(p)) file.size(p) else NA
      f_crc <- if (!is.null(p)) crc32_file(p) else NA_character_
      zt <- if (!is.null(p) && grepl("\\.zip$", p)) zip_test(p) else ""
      status <- if (is.null(a)) "new (not in archive)" else if (is.null(p)) "MISSING from download" else if (identical(a$crc, f_crc)) "identical" else {
        yy <- regmatches(name, regexec("BR(\\d\\d)\\.csv\\.zip$", name))[[1]]
        prelim <- (length(yy) && yy[2] %in% PRELIMINARY) || name == "2026.zip"
        if (prelim) "differs (preliminary/updated release)" else "differs (source re-published)"
      }
      rows[[length(rows) + 1]] <- data.table(archive = archive, file = name,
        archive_bytes = if (is.null(a)) NA_character_ else format(a$size, scientific = FALSE),
        download_bytes = if (is.na(f_size)) NA_character_ else format(f_size, scientific = FALSE),
        archive_crc32 = if (is.null(a)) NA_character_ else a$crc, download_crc32 = f_crc, zip_test = zt, status = status)
    }
  }
  tab <- rbindlist(rows)
  out <- file.path(raw, "archive_comparison.csv")
  write_csv_py(tab, out)
  cat(sprintf("Wrote %s\n", out))
  counts <- table(tab$status)
  for (k in sort_c(names(counts))) cat(sprintf("  %3d  %s\n", counts[[k]], k))
  for (i in which(tab$status != "identical" | tab$zip_test == "CORRUPT"))
    cat(sprintf("  - %s: %s (%s -> %s bytes) %s\n", tab$file[i], tab$status[i], tab$archive_bytes[i], tab$download_bytes[i], tab$zip_test[i]))
}

if (sys.nframe() == 0L) main()
