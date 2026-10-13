# Re-download every raw input of the arbovirus paper from its official source
# and record provenance (URL, size, sha256, Last-Modified, access time).
#
# Sources (verified 2026-10-08):
#   SINAN  - OpenDataSUS notification microdata (CSV zips), DENV/CHIKV/ZIKV 2016-2026
#   INMET  - Banco de Dados Historicos, automatic stations, yearly zips 2016-2026
#   IBGE   - population estimates (estimativa_dou / POP / serie TCU) + 2024 projections
#   GEO    - IBGE malhas API, 27 UF boundaries (GeoJSON)
#
# Usage:
#   Rscript tools/fetch_raw_sources.R [--out DIR] [--only sinan,inmet,ibge,geo]
#
# Downloads are streamed to <file>.part and resumed with HTTP Range; files
# already complete are skipped (their hash is still recorded in the manifest).
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
quiet_library("curl", "data.table")

UA <- "Mozilla/5.0 (arbovirus-surveillance-brazil data fetch)"
YEARS <- 2016:2026
SINAN_S3 <- "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SINAN"
SINAN_DISEASES <- c(DENG = "Dengue", CHIK = "Chikungunya", ZIKA = "Zikavirus")
INMET <- "https://portal.inmet.gov.br/uploads/dadoshistoricos"
IBGE_EST <- "https://ftp.ibge.gov.br/Estimativas_de_Populacao"
IBGE_PROJ <- "https://ftp.ibge.gov.br/Projecao_da_Populacao/Projecao_da_Populacao_2024"
GEO_API <- "https://servicodados.ibge.gov.br/api/v3/malhas/estados/%s?formato=application/vnd.geo+json&qualidade=intermediaria"
UF_CODES <- c("11", "12", "13", "14", "15", "16", "17", "21", "22", "23", "24", "25", "26", "27", "28", "29",
              "31", "32", "33", "35", "41", "42", "43", "50", "51", "52", "53")

# Exactly the IBGE files the results were built from, by FTP folder year.
# 2022/2023 have no estimate release (Censo 2022 year); the 2024 projections cover them.
IBGE_FILES <- list(
  "2016" = c("estimativa_dou_2016_20160913.xlsx", "estimativa_TCU_2016_20170614.xls", "serie_2001_2016_TCU.xls"),
  "2017" = c("estimativa_dou_2017.xls", "POP2017_20220905.xls", "serie_2001_2017_TCU.xls"),
  "2018" = c("estimativa_dou_2018_20181019.xls", "POP2018_20220905.xls", "serie_2001_2018_TCU.xls"),
  "2019" = c("estimativa_dou_2019.xls", "POP2019_20220905.xls", "serie_2001_2019_TCU.xls"),
  "2020" = c("estimativa_dou_2020.xls", "POP2020_20220711.xls", "POP2020_20220905.xls", "serie_2001_2020_TCU.xls"),
  "2021" = c("estimativa_dou_2021.xls", "POP2021_20240624.xls"),
  "2024" = c("estimativa_dou_2024.xls", "POP2024_20241101.xls", "POP2024_20241230.xls"),
  "2025" = c("estimativa_dou_2025.xls", "POP2025_20251031.xls", "POP2025_20260113.xls", "POP2025_20260828.xls"),
  "2026" = c("estimativa_dou_2026.xlsx"))
IBGE_PROJ_FILES <- c("projecoes_2024_tab1_idade_simples.xlsx", "projecoes_2024_tab4_indicadores.xlsx")

build_jobs <- function() {
  jobs <- list()
  add <- function(group, url, rel) jobs[[length(jobs) + 1]] <<- list(group = group, url = url, rel = rel)
  for (code in names(SINAN_DISEASES)) for (y in YEARS) {
    name <- sprintf("%sBR%02d.csv.zip", code, y %% 100)
    add("sinan", sprintf("%s/%s/csv/%s", SINAN_S3, SINAN_DISEASES[[code]], name), file.path("sinan", name))
  }
  for (y in YEARS) add("inmet", sprintf("%s/%d.zip", INMET, y), sprintf("inmet/%d.zip", y))
  for (y in names(IBGE_FILES)) for (name in IBGE_FILES[[y]])
    add("ibge", sprintf("%s/Estimativas_%s/%s", IBGE_EST, y, name), file.path("ibge", name))
  for (name in IBGE_PROJ_FILES) add("ibge", sprintf("%s/%s", IBGE_PROJ, name), file.path("ibge", name))
  for (code in UF_CODES) add("geo", sprintf(GEO_API, code), sprintf("geo/uf_%s.geojson", code))
  jobs
}

new_handle_ua <- function(...) curl::new_handle(useragent = UA, followlocation = TRUE, ...)

# Size and Last-Modified via a 1-byte ranged GET (INMET rejects HEAD).
remote_meta <- function(url) {
  r <- curl::curl_fetch_memory(url, handle = new_handle_ua(range = "0-0", timeout = 60))
  if (!r$status_code %in% c(200, 206)) stop(sprintf("HTTP %d", r$status_code))
  h <- curl::parse_headers_list(r$headers)
  size <- NA_real_
  cr <- h[["content-range"]]
  if (!is.null(cr) && grepl("/", cr) && !grepl("/\\*$", cr)) size <- as.numeric(sub("^.*/", "", cr))
  else if (r$status_code == 200 && !is.null(h[["content-length"]])) size <- as.numeric(h[["content-length"]])
  list(size = size, last_modified = if (is.null(h[["last-modified"]])) "" else h[["last-modified"]], ranged = r$status_code == 206)
}

download <- function(url, dest, retries = 5) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  meta <- remote_meta(url)
  if (file.exists(dest) && (is.na(meta$size) || file.size(dest) == meta$size)) return(list(status = "skipped", last_modified = meta$last_modified))
  part <- paste0(dest, ".part"); chunk <- paste0(dest, ".chunk")
  done <- FALSE
  for (attempt in seq_len(retries)) {
    have <- if (file.exists(part)) file.size(part) else 0
    if (have && !meta$ranged) { unlink(part); have <- 0 }
    h <- if (have) new_handle_ua(range = sprintf("%.0f-", have), timeout = 0) else new_handle_ua(timeout = 0)
    ok <- tryCatch({
      r <- curl::curl_fetch_disk(url, chunk, handle = h)
      if (r$status_code == 416) { unlink(chunk); TRUE }                      # already complete
      else {
        if (r$status_code >= 400) stop(sprintf("HTTP %d", r$status_code))
        if (r$status_code == 206 && file.exists(part)) { file.append(part, chunk); unlink(chunk) }
        else file.rename(chunk, part)
        is.na(meta$size) || file.size(part) == meta$size
      }
    }, error = function(e) { cat(sprintf("    retry %d/%d: %s\n", attempt, retries, conditionMessage(e))); Sys.sleep(min(60, 5 * attempt)); FALSE })
    if (ok) { done <- TRUE; break }
  }
  if (!done) stop("download incomplete after retries")
  if (!is.na(meta$size) && file.size(part) != meta$size) stop(sprintf("size mismatch %.0f != %.0f", file.size(part), meta$size))
  file.rename(part, dest)
  list(status = "downloaded", last_modified = meta$last_modified)
}

main <- function() {
  out <- normalizePath(cli_value("out", path_root("raw")), winslash = "/", mustWork = FALSE)
  groups <- strsplit(cli_value("only", "sinan,inmet,ibge,geo"), ",", fixed = TRUE)[[1]]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  manifest <- file.path(out, "sources_manifest.csv")
  rows <- list()
  if (file.exists(manifest)) {
    old <- fread(manifest, colClasses = "character", na.strings = NULL)
    for (i in seq_len(nrow(old))) rows[[old$file[i]]] <- as.list(old[i])
  }
  failures <- character()
  jobs <- Filter(function(j) j$group %in% groups, build_jobs())
  for (i in seq_along(jobs)) {
    j <- jobs[[i]]; dest <- file.path(out, j$rel)
    cat(sprintf("[%d/%d] %s\n", i, length(jobs), j$rel)); flush.console()
    res <- tryCatch(download(j$url, dest), error = function(e) { cat(sprintf("    FAILED: %s\n", conditionMessage(e))); NULL })
    if (is.null(res)) { failures <- c(failures, j$rel); next }
    prev <- rows[[j$rel]]
    rows[[j$rel]] <- list(source = j$group, file = j$rel, url = j$url, bytes = format(file.size(dest), scientific = FALSE),
      sha256 = digest::digest(file = dest, algo = "sha256"),
      last_modified = if (nzchar(res$last_modified)) res$last_modified else if (!is.null(prev)) prev$last_modified else "",
      accessed_utc = if (res$status == "skipped" && !is.null(prev)) prev$accessed_utc else format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
    cat(sprintf("    %s (%s bytes)\n", res$status, format(file.size(dest), big.mark = ",", scientific = FALSE)))
    # write after every file so an interrupted run keeps its provenance
    tab <- rbindlist(rows[sort_c(names(rows))])
    write_csv_py(tab, manifest)
  }
  cat(sprintf("\nManifest: %s (%d files)\n", manifest, length(rows)))
  if (length(failures)) { cat("FAILED:\n", paste0("  ", failures, "\n"), sep = ""); quit(status = 1) }
}

if (sys.nframe() == 0L) main()
