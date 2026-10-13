# One pass over raw/sinan/*.csv.zip -> db/normalized_fact_cases.csv + clinical_tables/.
# Cases are grouped by SG_UF_NOT (notification state) and month of DT_NOTIFIC; Year = file year.
# Each zip is unpacked to a temp file, only the needed columns are read (as text), and the temp file
# is deleted again, so at most one uncompressed SINAN file (<= ~1.9 GB) is on disk at a time.
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
quiet_library("data.table", "writexl")

IBGE_CODE_TO_NAME <- c("11" = "Rondonia", "12" = "Acre", "13" = "Amazonas", "14" = "Roraima", "15" = "Para", "16" = "Amapa",
  "17" = "Tocantins", "21" = "Maranhao", "22" = "Piaui", "23" = "Ceara", "24" = "Rio Grande do Norte", "25" = "Paraiba",
  "26" = "Pernambuco", "27" = "Alagoas", "28" = "Sergipe", "29" = "Bahia", "31" = "Minas Gerais", "32" = "Espirito Santo",
  "33" = "Rio de Janeiro", "35" = "Sao Paulo", "41" = "Parana", "42" = "Santa Catarina", "43" = "Rio Grande do Sul",
  "50" = "Mato Grosso do Sul", "51" = "Mato Grosso", "52" = "Goias", "53" = "Distrito Federal")
dim_loc <- fread(DIM_LOCATION, encoding = "UTF-8")
name2id <- setNames(dim_loc$Location_ID, dim_loc$Location_Name)
id2name <- setNames(dim_loc$Location_Name, dim_loc$Location_ID)
code2id <- setNames(unname(name2id[IBGE_CODE_TO_NAME]), names(IBGE_CODE_TO_NAME))
DISEASE <- c(DENG = "DENV", CHIK = "CHIKV", ZIKA = "ZIKV")

# pandas read_csv default NA markers (the Python version read with dtype=str and these defaults)
PD_NA <- c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN",
           "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null")

# ---- SINAN code decodings (documented in README) ----
MAP_SEX  <- c(M = "Masculino", F = "Feminino", I = "Ignorado")
MAP_CRIT <- c("1" = "Laboratorial", "2" = "Clinico-epidemiologico", "3" = "Em investigacao")
MAP_EVO  <- c("1" = "Cura", "2" = "Obito pelo agravo", "3" = "Obito por outras causas", "4" = "Obito em investigacao", "9" = "Ignorado")
MAP_HOSP <- c("1" = "Sim", "2" = "Nao", "9" = "Ignorado")
MAP_SORO <- c("1" = "DENV1", "2" = "DENV2", "3" = "DENV3", "4" = "DENV4")
MAP_GEST <- c("1" = "1o trimestre", "2" = "2o trimestre", "3" = "3o trimestre", "4" = "Idade gestacional ignorada",
              "5" = "Nao", "6" = "Nao se aplica", "9" = "Ignorado", "0" = "Nao")
MAP_CLASSI <- c("1" = "Dengue classico (schema antigo)", "2" = "Dengue c/ complicacoes (antigo)", "3" = "FHD (antigo)",
  "4" = "SCD (antigo)", "5" = "Descartado", "8" = "Inconclusivo", "10" = "Dengue", "11" = "Dengue c/ sinais de alarme",
  "12" = "Dengue grave", "13" = "Chikungunya")
MAP_LAB <- c("1" = "Reagente/Positivo", "2" = "Nao reagente/Negativo", "3" = "Inconclusivo", "4" = "Nao realizado", "9" = "Ignorado")
AGE_BINS   <- c(0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 999)
AGE_LABELS <- c("<1", "1-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

# NU_IDADE_N "4yyy" = years, "1/2/3xxx" = hours/days/months (< 1 year). Anything else -> NA.
# Applied to the distinct raw values only; returns the age band label (NA = unknown).
age_band <- function(u) {
  x <- py_strip(ifelse(is.na(u), "", u))
  ok <- grepl("^[0-9]{1,4}$", x)
  x <- ifelse(ok, paste0(strrep("0", pmax(0, 4 - nchar(x))), x), NA_character_)   # str.zfill(4)
  unit <- substr(x, 1, 1)
  val <- suppressWarnings(as.numeric(substr(x, 2, 4)))
  yrs <- ifelse(!is.na(unit) & unit == "4", val, ifelse(!is.na(unit) & unit %in% c("1", "2", "3"), 0, NA))
  k <- findInterval(yrs, AGE_BINS)            # [b_i, b_i+1)
  k[!is.na(k) & (k < 1 | k > length(AGE_LABELS))] <- NA
  AGE_LABELS[k]
}

# Python str.isdigit() over the Latin-1 range: 0-9 plus superscript 1/2/3
is_py_digit <- function(v) grepl("^[0-9²³¹]+$", v, perl = TRUE)

# fillna('NA') -> strip -> '' as 'NA' -> label / 'Ignorado/NA' / 'codigo_<digits>'
mapped <- function(u, mapping) {
  s <- ifelse(is.na(u), "NA", u)
  s <- py_strip(s)
  s[s == ""] <- "NA"
  out <- unname(mapping[s])
  miss <- is.na(out)
  out[miss] <- ifelse(s[miss] %in% c("NA", "nan"), "Ignorado/NA",
                      ifelse(is_py_digit(s[miss]), paste0("codigo_", s[miss]), "Ignorado/NA"))
  out
}

# Month of an ISO date string (pandas to_datetime(errors="coerce")); NA when not a valid date.
iso_month <- function(u) {
  m <- rep(NA_integer_, length(u))
  ok <- !is.na(u) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", u)
  d <- as.Date(u[ok], format = "%Y-%m-%d")
  mm <- as.integer(format(d, "%m"))
  m[ok] <- mm
  m
}

DISCARD <- list(DENV = "5", CHIKV = character(), ZIKV = character())
COMMON <- list(list("sex", "CS_SEXO", MAP_SEX), list("criterion", "CRITERIO", MAP_CRIT),
               list("evolution", "EVOLUCAO", MAP_EVO), list("classification", "CLASSI_FIN", MAP_CLASSI))
EXTRA <- list(DENV = list(list("hospitalization", "HOSPITALIZ", MAP_HOSP), list("serotype", "SOROTIPO", MAP_SORO),
                          list("ns1_result", "RESUL_NS1", MAP_LAB), list("pcr_result", "RESUL_PCR_", MAP_LAB),
                          list("serology_igm_result", "RESUL_SORO", MAP_LAB)),
              ZIKV = list(list("pregnancy", "CS_GESTANT", MAP_GEST)),
              CHIKV = list())

cases_parts <- list(); clin_rows <- list(); dropped <- 0; total_rows <- 0
t0 <- Sys.time()
tmpdir <- file.path(tempdir(), "sinan"); dir.create(tmpdir, showWarnings = FALSE)
zips <- sort_c(list.files(SINAN_DIR, pattern = "BR.*\\.csv\\.zip$", full.names = TRUE))
non_iso <- 0
for (zf in zips) {
  base <- basename(zf); pref <- substr(base, 1, 4); yr <- 2000L + as.integer(substr(base, 7, 8)); dis <- DISEASE[[pref]]
  entry <- unzip(zf, list = TRUE)$Name[1]
  unzip(zf, files = entry, exdir = tmpdir, junkpaths = TRUE)
  csv <- file.path(tmpdir, basename(entry))
  header <- names(fread(csv, nrows = 0, encoding = "Latin-1"))
  dims <- Filter(function(d) d[[2]] %in% header, c(COMMON, EXTRA[[dis]]))
  agecol <- if ("NU_IDADE_N" %in% header) "NU_IDADE_N" else NULL
  need <- unique(c("DT_NOTIFIC", "SG_UF_NOT", vapply(dims, `[[`, "", 2), agecol))
  need <- need[need %in% header]
  ch <- fread(csv, select = need, colClasses = "character", encoding = "Latin-1", na.strings = PD_NA,
              strip.white = FALSE, showProgress = FALSE)
  unlink(csv)
  nrows <- nrow(ch); n_disc <- 0
  drop_codes <- DISCARD[[dis]]
  # apply the per-disease case definition to the whole analysis population
  if (length(drop_codes) && "CLASSI_FIN" %in% names(ch)) {
    u <- unique(ch$CLASSI_FIN)
    drop_u <- py_strip(ifelse(is.na(u), "nan", u)) %in% drop_codes
    drop <- drop_u[match(ch$CLASSI_FIN, u)]
    n_disc <- sum(drop); ch <- ch[!drop]
  }
  # cases (valid state + month) and federative_unit (valid state)
  g <- ch[, .N, by = .(SG_UF_NOT, DT_NOTIFIC)]
  g[, lid := unname(code2id[SG_UF_NOT])]
  ud <- unique(g$DT_NOTIFIC); um <- iso_month(ud)
  non_iso <- non_iso + sum(g$N[!is.na(g$DT_NOTIFIC) & is.na(um[match(g$DT_NOTIFIC, ud)])])
  g[, Month := um[match(DT_NOTIFIC, ud)]]
  c_cases <- g[!is.na(lid) & !is.na(Month), .(Cases = sum(N)), by = .(lid, Month)]
  c_fed <- g[!is.na(lid), .(Count = sum(N)), by = lid]
  if (nrow(c_cases))
    cases_parts[[length(cases_parts) + 1]] <- data.table(Location_ID = as.integer(c_cases$lid), Disease = dis, Year = yr,
                                                         Month = as.integer(c_cases$Month), Cases = as.integer(c_cases$Cases))
  # clinical dims
  for (d in dims) {
    vc <- ch[, .N, by = c(d[[2]])]
    vc[, Category := mapped(vc[[1]], d[[3]])]
    vc <- vc[, .(Count = sum(N)), by = Category]
    clin_rows[[length(clin_rows) + 1]] <- data.table(Disease = dis, Year = yr, Dimension = d[[1]], vc)
  }
  # age: every band is reported (zero counts included), plus 'Ignorado/NA' when any age is unknown
  if (!is.null(agecol)) {
    va <- ch[, .N, by = c(agecol)]
    va[, band := age_band(va[[1]])]
    counts <- setNames(integer(length(AGE_LABELS)), AGE_LABELS)
    ag <- va[!is.na(band), .(Count = sum(N)), by = band]
    counts[ag$band] <- ag$Count
    age <- data.table(Category = AGE_LABELS, Count = as.integer(counts))
    n_na <- sum(va$N[is.na(va$band)])
    if (n_na > 0) age <- rbind(age, data.table(Category = "Ignorado/NA", Count = as.integer(n_na)))
    clin_rows[[length(clin_rows) + 1]] <- data.table(Disease = dis, Year = yr, Dimension = "age_group", age)
  }
  if (nrow(c_fed))
    clin_rows[[length(clin_rows) + 1]] <- data.table(Disease = dis, Year = yr, Dimension = "federative_unit",
                                                     Category = unname(id2name[as.character(c_fed$lid)]), Count = as.integer(c_fed$Count))
  total_rows <- total_rows + nrows
  valid <- sum(c_cases$Cases); dropped <- dropped + (nrows - n_disc - valid)
  cat(sprintf("  %-18s %-5s %d notif=%8d discarded=%8d case_pop=%8d\n", base, dis, yr, nrows, n_disc, valid))
  flush.console()
  rm(ch, g); invisible(gc(verbose = FALSE))
}
if (non_iso) cat(sprintf("  note: %d notifications with a non-ISO DT_NOTIFIC were left without a month\n", non_iso))

# ---- fact_cases.csv ----
fc <- rbindlist(cases_parts)
pop <- fread(file.path(DB, "normalized_fact_population.csv"))
fc <- merge(fc, pop, by = c("Location_ID", "Year"), all.x = TRUE, sort = FALSE)
fc[, Monthly_Inc_per100k := py_round(Cases / Population * 1e5, 2)]
fc <- fc[order_c(Disease, Location_ID, Year, Month), .(Location_ID, Disease, Year, Month, Cases, Monthly_Inc_per100k)]
fc[, Cases := as.double(Cases)]
write_csv_py(fc, file.path(DB, "normalized_fact_cases.csv"))
cat(sprintf("\nnormalized_fact_cases.csv: (%d, %d) | rows_read=%s dropped_no_state_or_month=%s\n", nrow(fc), ncol(fc),
            format(total_rows, big.mark = ","), format(dropped, big.mark = ",")))

# ---- clinical long + per-disease xlsx ----
cl <- rbindlist(clin_rows)
cl <- cl[, .(Count = as.integer(sum(Count))), by = .(Disease, Year, Dimension, Category)]
cl <- cl[order_c(Disease, Dimension, Year, Category)]
write_csv_py(cl, file.path(CLIN, "clinical_demographic_long.csv"))
cat("clinical_demographic_long.csv: (", nrow(cl), ", ", ncol(cl), ") dims: ",
    paste(sort_c(unique(cl$Dimension)), collapse = ", "), "\n", sep = "")
for (dis in c("DENV", "CHIKV", "ZIKV")) {
  sub <- cl[Disease == dis]
  sheets <- list()
  for (dname in sort_c(unique(sub$Dimension))) {
    g <- sub[Dimension == dname]
    piv <- dcast(g, Category ~ Year, value.var = "Count", fun.aggregate = sum, fill = 0L)
    piv <- piv[order_c(Category)]
    sheets[[dname]] <- piv
  }
  write_xlsx_sheets(sheets, file.path(CLIN, sprintf("clinical_%s.xlsx", dis)))
  cat(sprintf("  clinical_%s.xlsx written (%d sheets)\n", dis, length(sheets)))
}

# ---- verification vs committed ----
comm <- fread(file.path(ORIG_DB, "normalized_fact_cases.csv"))
cmp <- merge(fc[, .(new = sum(Cases)), by = .(Disease, Year)],
             comm[, .(committed = sum(Cases)), by = .(Disease, Year)], by = c("Disease", "Year"), all.x = TRUE)
cmp[, `diff_%` := py_round((new - committed) / committed * 100, 2)]
cat("\n=== fact_cases overlap-year totals vs committed (per-disease definition) ===\n")
print_df(cmp)
ov <- cmp[!is.na(committed)]
cat(sprintf("\noverlap years: %d | within 1%%: %d | within 5%%: %d | >5%%: %d\n", nrow(ov),
            sum(abs(ov$`diff_%`) <= 1), sum(abs(ov$`diff_%`) <= 5), sum(abs(ov$`diff_%`) > 5)))
cat("years >5% (data-vintage differences):\n")
print_df(ov[abs(`diff_%`) > 5])
cat(sprintf("\nelapsed %.0fs\n", as.numeric(Sys.time() - t0, units = "secs")))
