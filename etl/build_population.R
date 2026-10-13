# Build db/normalized_fact_population.csv (2016-2026) entirely from raw IBGE files:
#   2016-2023  IBGE Projecoes da Populacao, Revisao 2024 (raw/ibge/projecoes_2024_tab4_indicadores.xlsx,
#              POP_T). Census-2022-consistent; IBGE published no estimates for 2022-2023.
#   2024-2026  IBGE official 1-July estimates published in the DOU (raw/ibge/estimativa_dou_<year>.xls(x)).
# The pre-census DOU estimates for 2016-2021 are not used: they overstate several states (Amapa +11%)
# and would create a ~2.4 M drop between 2021 and 2022.
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
quiet_library("data.table", "readxl")

PROJ_YEARS <- 2016:2023
DOU_YEARS  <- 2024:2026

# '3.273.227(1)' / '3.212.180 (*)' / 1787279 -> integer; NA if not a population value.
to_int <- function(v) {
  if (is.numeric(v)) return(trunc(v))
  s <- if (inherits(v, "POSIXt")) format(v, "%Y-%m-%d %H:%M:%S") else as.character(v)
  s <- gsub("\\(.*?\\)", "", s, perl = TRUE)          # drop footnote marks before keeping digits
  s <- gsub("\\D", "", s, perl = TRUE)
  if (nzchar(s)) as.numeric(s) else NA_real_
}

dim_loc <- fread(DIM_LOCATION, encoding = "UTF-8")
name2id <- setNames(dim_loc$Location_ID, strip_accents(dim_loc$Location_Name))

parse_dou <- function(path) {
  sheet <- suppressMessages(read_excel(path, sheet = 1, col_names = FALSE, col_types = "list"))
  rows <- list()
  for (i in seq_len(nrow(sheet))) {
    first <- sheet[[1]][[i]]
    key <- strip_accents(if (is.null(first) || (length(first) == 1 && is.na(first))) "nan" else first)
    if (!key %in% names(name2id)) next
    cells <- lapply(sheet[i, -1], `[[`, 1)
    cells <- cells[!vapply(cells, function(v) is.null(v) || (length(v) == 1 && is.na(v)), logical(1))]
    vals <- vapply(cells, to_int, numeric(1))
    vals <- vals[!is.na(vals) & vals != 0]
    id <- as.character(name2id[[key]])
    if (length(vals) && is.null(rows[[id]])) rows[[id]] <- vals[[1]]   # first block is the UF table
  }
  rows
}

parts <- list()
cat(sprintf("%4s  %-40s %6s %14s\n", "year", "source", "states", "Brazil total"))
fmt_total <- function(x) formatC(x, format = "d", big.mark = ",")

# ---- 2016-2023: projections ----
proj_file <- "projecoes_2024_tab4_indicadores.xlsx"
proj <- as.data.table(suppressMessages(read_excel(file.path(IBGE_DIR, proj_file), skip = 6)))
proj[, Location_ID := unname(name2id[strip_accents(LOCAL)])]
proj[, Year := suppressWarnings(as.numeric(ANO))]
proj <- proj[!is.na(Location_ID) & !is.na(Year)]
for (y in PROJ_YEARS) {
  df <- proj[Year == y]
  df <- data.table(Location_ID = as.integer(df$Location_ID), Year = as.integer(y),
                   Population = as.integer(py_round(as.numeric(df$POP_T))))
  if (nrow(df) != 27) stop(sprintf("%s: %d states for %d, expected 27", proj_file, nrow(df), y))
  cat(sprintf("%4d  %-40s %6d %14s\n", y, proj_file, nrow(df), fmt_total(sum(as.numeric(df$Population)))))
  parts[[length(parts) + 1]] <- df
}

# ---- 2024-2026: DOU estimates ----
for (y in DOU_YEARS) {
  files <- sort_c(list.files(IBGE_DIR, pattern = paste0("^estimativa_dou_", y, ".*\\.xls.*$")))
  if (!length(files)) stop(sprintf("missing raw/ibge/estimativa_dou_%d.xls(x)", y))
  rows <- parse_dou(file.path(IBGE_DIR, files[1]))
  if (length(rows) != 27) stop(sprintf("%s: parsed %d states, expected 27", files[1], length(rows)))
  df <- data.table(Location_ID = as.integer(names(rows)), Year = as.integer(y),
                   Population = as.integer(unlist(rows)))
  cat(sprintf("%4d  %-40s %6d %14s\n", y, files[1], nrow(df), fmt_total(sum(as.numeric(df$Population)))))
  parts[[length(parts) + 1]] <- df
}

out <- rbindlist(parts)
setorder(out, Location_ID, Year)
write_csv_py(out, file.path(DB, "normalized_fact_population.csv"))
cat(sprintf("\nnormalized_fact_population.csv: (%d, %d) years %d - %d\n", nrow(out), ncol(out), min(out$Year), max(out$Year)))

# ---- comparison with the series used before ----
old <- fread(file.path(ORIG_DB, "normalized_fact_population.csv"))
m <- merge(out, old, by = c("Location_ID", "Year"), suffixes = c("_new", "_old"))
cat(sprintf("vs previous series (2016-2025): %d/%d state-years identical\n",
            sum(m$Population_new == m$Population_old), nrow(m)))
