# Aggregate raw/inmet/<year>.zip hourly station CSVs -> db/normalized_fact_climate.csv (state x month, 2016-2026).
# Temp_Avg = mean across UF stations of each station's monthly mean hourly air temp.
# Rain_mm  = mean across UF stations of each station's monthly total precip.
# Coverage filter: temp needs >=240 valid hourly obs/month; rain needs >=360 hourly obs/month.
# Keeps only the latest snapshot per station when a year zip holds several.
# Sums and means use the same compensated (Kahan) summation, in the same row order, as the
# pandas groupby the published db/ was built with, so rounding to 2 decimals lands identically.
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

dim_loc <- fread(DIM_LOCATION, encoding = "UTF-8")
name2id <- setNames(dim_loc$Location_ID, dim_loc$Location_Name)
NAME_TO_ABBR <- c("Rondonia" = "RO", "Acre" = "AC", "Amazonas" = "AM", "Roraima" = "RR", "Para" = "PA", "Amapa" = "AP",
  "Tocantins" = "TO", "Maranhao" = "MA", "Piaui" = "PI", "Ceara" = "CE", "Rio Grande do Norte" = "RN", "Paraiba" = "PB",
  "Pernambuco" = "PE", "Alagoas" = "AL", "Sergipe" = "SE", "Bahia" = "BA", "Minas Gerais" = "MG", "Espirito Santo" = "ES",
  "Rio de Janeiro" = "RJ", "Sao Paulo" = "SP", "Parana" = "PR", "Santa Catarina" = "SC", "Rio Grande do Sul" = "RS",
  "Mato Grosso do Sul" = "MS", "Mato Grosso" = "MT", "Goias" = "GO", "Distrito Federal" = "DF")
abbr2id <- setNames(unname(name2id[names(NAME_TO_ABBR)]), unname(NAME_TO_ABBR))
MIN_T <- 240; MIN_R <- 360

# Decimal-comma numbers ("21,5", ",4", "-9999") parsed as integer mantissa / 10^k, which is how
# pandas' C parser converts them (correctly rounded). Blank -> NA.
parse_decimal_comma <- function(tok) {
  u <- unique(tok)
  m <- regmatches(u, regexec("^(-?)([0-9]*),?([0-9]*)$", u))
  val <- vapply(m, function(p) {
    if (!length(p) || !nzchar(paste0(p[3], p[4]))) return(NA_real_)
    v <- as.numeric(paste0(p[3], p[4])) / 10^nchar(p[4])
    if (p[2] == "-") -v else v
  }, numeric(1))
  val[match(tok, u)]
}

# Year / month from "YYYY-MM-DD" or "YYYY/MM/DD" (pandas to_datetime(format="mixed")); NA if invalid.
parse_ym <- function(tok) {
  u <- unique(tok)
  ok <- !is.na(u) & grepl("^[0-9]{4}[-/][0-9]{2}[-/][0-9]{2}$", u)
  d <- as.Date(rep(NA_character_, length(u)))
  d[ok] <- as.Date(gsub("/", "-", u[ok], fixed = TRUE), format = "%Y-%m-%d")
  y <- as.integer(format(d, "%Y")); mo <- as.integer(format(d, "%m"))
  i <- match(tok, u)
  list(Y = y[i], M = mo[i])
}

# A year zip can hold several cumulative snapshots of the same station
# (2026.zip: ..._A_31-08-2026.CSV and ..._A_30-09-2026.CSV). Keep the one ending latest
# (the first one wins a tie), listed in order of each station's first appearance.
latest_per_station <- function(names) {
  names <- names[grepl("\\.CSV$", toupper(names))]
  parts <- strsplit(basename(names), "_", fixed = TRUE)
  ok <- lengths(parts) >= 4
  names <- names[ok]; parts <- parts[ok]
  code <- vapply(parts, `[[`, "", 4)
  m <- regmatches(names, regexec("_A_([0-9]{2})-([0-9]{2})-([0-9]{4})", names))
  end <- vapply(m, function(p) if (length(p)) as.numeric(paste0(p[4], p[3], p[2])) else -Inf, numeric(1))
  keep <- character()
  for (cd in unique(code)) {
    i <- which(code == cd)
    keep <- c(keep, names[i[which.max(end[i])]])     # which.max returns the first maximum
  }
  keep
}

rows <- list()
t0 <- Sys.time()
tmp <- file.path(tempdir(), "inmet")
for (yzip in sort_c(list.files(INMET_DIR, pattern = "\\.zip$", full.names = TRUE))) {
  yr <- as.integer(substr(basename(yzip), 1, 4))
  all_csv <- unzip(yzip, list = TRUE)$Name
  all_csv <- all_csv[grepl("\\.CSV$", toupper(all_csv))]
  keep <- latest_per_station(all_csv)
  if (length(keep) < length(all_csv))
    cat(sprintf("  %d: %d older duplicate station snapshots skipped\n", yr, length(all_csv) - length(keep)))
  parts <- strsplit(basename(keep), "_", fixed = TRUE)
  abbr <- vapply(parts, function(p) if (length(p) >= 3) p[3] else NA_character_, "")
  use <- !is.na(abbr) & abbr %in% names(abbr2id)
  keep <- keep[use]; abbr <- abbr[use]
  unlink(tmp, recursive = TRUE); dir.create(tmp)
  unzip(yzip, files = keep, exdir = tmp)
  stations <- vector("list", length(keep))
  for (s in seq_along(keep)) {
    d <- tryCatch(fread(file.path(tmp, keep[s]), skip = 8, sep = ";", header = TRUE, select = c(1L, 3L, 8L),
                        colClasses = "character", encoding = "Latin-1", na.strings = NULL, strip.white = FALSE,
                        showProgress = FALSE),
                  error = function(e) NULL)
    if (is.null(d)) next
    setnames(d, c("Data", "Rain", "Temp"))
    d[, s := s]
    stations[[s]] <- d
  }
  unlink(tmp, recursive = TRUE)
  st <- rbindlist(stations)
  ym <- parse_ym(st$Data)
  st[, `:=`(Y = ym$Y, M = ym$M, Rain = parse_decimal_comma(Rain), Temp = parse_decimal_comma(Temp))]
  # INMET codes missing as -9999 (and blanks); mask sentinels / impossible values
  st[!is.na(Rain) & Rain < 0, Rain := NA_real_]
  st[!is.na(Temp) & (Temp < -50 | Temp > 60), Temp := NA_real_]
  st <- st[!is.na(M)]
  st[, gid := frankv(list(s, Y, M), ties.method = "dense")]
  G <- max(st$gid)
  tm <- kahan_by(st$Temp, st$gid, G, "mean")
  rs <- kahan_by(st$Rain, st$gid, G, "sum")
  grp <- unique(st[, .(gid, s, Y, M)])[order(gid)]
  rows[[length(rows) + 1]] <- data.table(Abbr = abbr[grp$s], Year = grp$Y, Month = grp$M,
                                         Temp_mean = tm$value[grp$gid], Temp_n = tm$n[grp$gid],
                                         Rain_sum = rs$value[grp$gid], Rain_n = rs$n[grp$gid])
  cat(sprintf("  %d: %d stations parsed  (elapsed %.0fs)\n", yr, uniqueN(grp$s),
              as.numeric(Sys.time() - t0, units = "secs")))
  flush.console()
}

st <- rbindlist(rows)
st <- st[Year >= 2016 & Year <= 2026]
station_mean <- function(d, col, out) {
  d <- copy(d)
  d[, gid := frankv(list(Abbr, Year, Month), ties.method = "dense")]
  k <- kahan_by(d[[col]], d$gid, max(d$gid), "mean")
  g <- unique(d[, .(gid, Abbr, Year, Month)])[order(gid)]
  g[, (out) := k$value[gid]]
  g[, gid := NULL]
  g
}
temp <- station_mean(st[Temp_n >= MIN_T], "Temp_mean", "Temp_Avg")
rain <- station_mean(st[Rain_n >= MIN_R], "Rain_sum", "Rain_mm")
clim <- merge(temp, rain, by = c("Abbr", "Year", "Month"), all = TRUE)
clim[, Location_ID := unname(abbr2id[Abbr])]
clim <- clim[!is.na(Location_ID)]
clim[, Location_ID := as.integer(Location_ID)]
clim <- clim[order(Location_ID, Year, Month)]
clim[, `:=`(Temp_Avg = py_round(Temp_Avg, 2), Rain_mm = py_round(Rain_mm, 2))]
# lags within state (positional: the previous row of the same state)
clim[, `:=`(Temp_Lag1 = py_round(shift(Temp_Avg, 1), 2), Temp_Lag2 = py_round(shift(Temp_Avg, 2), 2),
            Rain_Lag1 = py_round(shift(Rain_mm, 1), 2), Rain_Lag2 = py_round(shift(Rain_mm, 2), 2)), by = Location_ID]
clim <- clim[, .(Location_ID, Year = as.integer(Year), Month = as.integer(Month), Temp_Avg, Rain_mm,
                 Temp_Lag1, Temp_Lag2, Rain_Lag1, Rain_Lag2)]
write_csv_py(clim, file.path(DB, "normalized_fact_climate.csv"))
cat(sprintf("\nnormalized_fact_climate.csv: (%d, %d) years %d-%d\n", nrow(clim), ncol(clim), min(clim$Year), max(clim$Year)))
cat("2026 months present:", paste(sort(unique(clim[Year == 2026]$Month)), collapse = ", "), "\n")

# ---- calibration vs committed ----
comm <- fread(file.path(ORIG_DB, "normalized_fact_climate.csv"))
m <- merge(comm, clim, by = c("Location_ID", "Year", "Month"), suffixes = c("_old", "_new"))
cat(sprintf("\n=== calibration vs committed (overlap cells=%d) ===\n", nrow(m)))
for (cc in c("Temp_Avg", "Rain_mm")) {
  o <- m[[paste0(cc, "_old")]]; n <- m[[paste0(cc, "_new")]]
  ok <- !is.na(o) & !is.na(n)
  cat(sprintf("  %s: corr=%.4f  MAE=%.2f  mean_old=%.2f mean_new=%.2f ratio_new/old=%.3f\n", cc, cor(o[ok], n[ok]),
              mean(abs(o[ok] - n[ok])), mean(o[ok]), mean(n[ok]), mean(n[ok]) / max(mean(o[ok]), 1e-9)))
}
cat(sprintf("elapsed %.0fs\n", as.numeric(Sys.time() - t0, units = "secs")))
