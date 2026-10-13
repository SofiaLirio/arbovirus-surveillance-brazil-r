# ================================================================================
#  ADV_COMMON — shared helpers for the additional-analyses pipeline
# ================================================================================
#  Reuses pipelines/data_loader.R for data access. Every analysis module sources
#  this to resolve the dated output root, build national/regional/state monthly
#  series generically, add seasonal harmonic terms, and save tables/figures.
#  Output root: outputs/<RUN_DATE>/additional/   (CSV + PNG together)
#  Sourced after R/utils.R.
# ================================================================================
source(path_root("pipelines", "data_loader.R"))
quiet_library("ggplot2", "patchwork")

RUN_DATE <- "2026-10-08"
OUT_ROOT <- file.path(OUTPUT_DIR, RUN_DATE, "additional")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

SEED <- 42
set.seed(SEED)

CLIM_COLS <- c("Temp_Avg", "Rain_mm", "Temp_Lag1", "Temp_Lag2", "Rain_Lag1", "Rain_Lag2")
LEVELS <- c("national", "region", "state")
# Minimum reported months for a per-geography time-series model to be attempted.
MIN_OBS <- 36

# ── saving ────────────────────────────────────────────────────────────────────
save_fig <- function(p, name, width, height, dpi = 150) {
  save_png(p, file.path(OUT_ROOT, name), width, height, dpi = dpi)
  cat(sprintf("  saved -> %s\n", name))
}
save_table <- function(df, name) {
  write_csv_py(df, file.path(OUT_ROOT, name))
  cat(sprintf("  saved -> %s (%d rows)\n", name, nrow(df)))
}

# ── series construction ───────────────────────────────────────────────────────
# Collapse a (geo) subset of the long per-state table into one row per Date.
# Cases/Population summed across states (NA only if all missing); climate is
# population-weighted; a complete monthly calendar is added so time-series models see gaps.
agg_to_monthly <- function(sub) {
  sub <- as.data.table(sub)
  wmean <- function(v, w) { ok <- !is.na(v) & !is.na(w) & w > 0; if (any(ok)) sum(v[ok] * w[ok]) / sum(w[ok]) else NA_real_ }
  res <- sub[, c(list(Year = Year[1], Month = Month[1], Cases = sum_min1(Cases), Population = sum_min1(as.numeric(Population))),
                 lapply(setNames(CLIM_COLS, CLIM_COLS), function(cc) wmean(get(cc), as.numeric(Population)))),
             by = Date][order(Date)]
  res[, Incidence := ifelse(!is.na(Population) & Population > 0 & !is.na(Cases), Cases / Population * 1e5, NA_real_)]
  if (!nrow(res)) return(res)
  full <- data.table(Date = seq(min(res$Date), max(res$Date), by = "month"))
  res <- merge(full, res, by = "Date", all.x = TRUE)[order(Date)]
  res[, `:=`(Year = as.integer(format(Date, "%Y")), Month = as.integer(format(Date, "%m")), t = seq_len(.N) - 1L)]
  setcolorder(res, c("Date", "Year", "Month", "Cases", "Population", CLIM_COLS, "Incidence", "t"))
  res
}

# Named list {geo: monthly table} for one virus at the given level.
monthly_series <- function(df, level) {
  groups <- switch(level,
    national = list(Brazil = df),
    region = split(df, df$Region)[sort_c(unique(df$Region))],
    state = split(df, df$Location_Name)[sort_c(unique(df$Location_Name))],
    stop(level))
  out <- lapply(groups, agg_to_monthly)
  out[vapply(out, nrow, 0L) > 0]
}

# List of list(virus, level, geo, s) across all viruses/levels.
iter_series <- function(datasets, levels = LEVELS) {
  out <- list()
  for (virus in names(datasets)) for (level in levels) {
    ms <- monthly_series(datasets[[virus]], level)
    for (geo in names(ms)) out[[length(out) + 1]] <- list(virus = virus, level = level, geo = geo, s = ms[[geo]])
  }
  out
}

# ── seasonality ───────────────────────────────────────────────────────────────
# sin/cos Fourier terms of annual period (12 months)
harmonic_terms <- function(month, K = 2) {
  month <- as.numeric(month)
  cols <- list()
  for (k in seq_len(K)) {
    cols[[paste0("sin", k)]] <- sin(2 * pi * k * month / 12)
    cols[[paste0("cos", k)]] <- cos(2 * pi * k * month / 12)
  }
  as.data.frame(cols)
}

load_datasets <- function() {
  cat("Loading datasets via data_loader ...\n")
  load_pipeline_datasets()
}

theme_adv <- function() theme_mpl(9) + theme(plot.title = element_text(size = 10, face = "plain"))
