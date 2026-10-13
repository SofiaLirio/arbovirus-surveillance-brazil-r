# ================================================================================
#  2. CO-CIRCULATION — cross-arbovirus synchrony
# ================================================================================
#  Dengue, Chikungunya and Zika share the Aedes vector, so their incidence series
#  co-move. For every geography (national, 5 regions, 27 states) and virus pair:
#    - cross-correlation of the deseasonalised series (12-month difference,
#      z-scored) for lags -6..+6 months (+k = first virus leads by k months),
#    - the lag with the largest |CCF| and a naive t-test p-value (FDR across states),
#    - co-occurrence of "active" months (incidence above the series' own median).
#  Outputs (outputs/<date>/additional/):
#    cocirculation_ccf.csv, cocirculation_bestlag.csv, cocirculation_cooccurrence.csv,
#    cocirculation_heatmap.png (national CCF curves + co-activity bars)
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
source(path_root("pipelines", "additional", "adv_common.R"))

PAIRS <- list(c("Dengue", "Chikungunya"), c("Dengue", "Zika"), c("Chikungunya", "Zika"))
MAX_LAG <- 6          # months
MIN_OVERLAP <- 36     # overlapping months required

# 12-month difference to strip the shared seasonal cycle, then z-score (sd with n-1).
deseason <- function(s) {
  x <- c(rep(NA_real_, 12), diff(as.numeric(s), lag = 12))[seq_along(s)]
  x <- x - mean(x, na.rm = TRUE)
  sd <- stats::sd(x, na.rm = TRUE)
  if (!is.na(sd) && sd > 0) x / sd else x
}

# Cross-correlation of a vs b for lag in [-max_lag, max_lag].
# Positive lag k = a leads b by k months (corr of a[t] with b[t+k]).
ccf_ab <- function(a, b, max_lag = MAX_LAG) {
  n <- length(a)
  out <- setNames(rep(NA_real_, 2 * max_lag + 1), -max_lag:max_lag)
  for (k in -max_lag:max_lag) {
    if (k >= 0) { x <- a[seq_len(n - k)]; y <- b[(k + 1):n] } else { x <- a[(-k + 1):n]; y <- b[seq_len(n + k)] }
    m <- is.finite(x) & is.finite(y)
    if (sum(m) >= MIN_OVERLAP) { r <- suppressWarnings(stats::cor(x[m], y[m])); out[as.character(k)] <- if (is.finite(r)) r else NA }
  }
  out
}

# {geo: table with Date + one incidence column per virus}, full outer join on Date
series_by_level <- function(datasets, level) {
  per_virus <- lapply(datasets, function(df) monthly_series(df, level))
  geos <- sort_c(unique(unlist(lapply(per_virus, names))))
  out <- list()
  for (geo in geos) {
    cols <- list()
    for (v in names(datasets)) {
      s <- per_virus[[v]][[geo]]
      if (!is.null(s)) cols[[v]] <- setnames(s[, .(Date, Incidence)], "Incidence", v)
    }
    if (length(cols) >= 2) out[[geo]] <- Reduce(function(x, y) merge(x, y, by = "Date", all = TRUE), cols)[order(Date)]
  }
  out
}

run <- function() {
  datasets <- load_datasets()
  ccf_rows <- list(); best_rows <- list(); cooc_rows <- list()
  for (level in LEVELS) {
    sbl <- series_by_level(datasets, level)
    for (geo in names(sbl)) {
      wide <- sbl[[geo]]
      for (pr in PAIRS) {
        v1 <- pr[1]; v2 <- pr[2]
        if (!(v1 %in% names(wide) && v2 %in% names(wide))) next
        a <- deseason(wide[[v1]]); b <- deseason(wide[[v2]])
        cc <- ccf_ab(a, b)
        if (all(is.na(cc))) next
        pair <- paste0(v1, "-", v2)
        ccf_rows[[length(ccf_rows) + 1]] <- data.table(Level = level, Geo = geo, Pair = pair, Lag_months = as.integer(names(cc)), CCF = unname(cc))
        fin <- cc[!is.na(cc)]
        best_lag <- as.integer(names(fin)[which.max(abs(fin))])    # first maximum, scanning from -6
        r <- fin[[as.character(best_lag)]]
        n_eff <- sum(is.finite(a))
        dfree <- max(n_eff - 2, 1)
        p <- if (abs(r) < 1) 2 * stats::pt(abs(r) * sqrt(dfree / max(1 - r^2, 1e-9)), dfree, lower.tail = FALSE) else 0
        best_rows[[length(best_rows) + 1]] <- data.table(Level = level, Geo = geo, Pair = pair, Best_lag_months = best_lag,
                                                         Peak_CCF = r, p_value = p, n = n_eff)
        # co-occurrence on raw incidence
        raw <- wide[!is.na(get(v1)) & !is.na(get(v2)), c(v1, v2), with = FALSE]
        if (nrow(raw) >= MIN_OVERLAP) {
          act1 <- raw[[v1]] > stats::median(raw[[v1]]); act2 <- raw[[v2]] > stats::median(raw[[v2]])
          cooc_rows[[length(cooc_rows) + 1]] <- data.table(Level = level, Geo = geo, Pair = pair,
            Frac_both_active = mean(act1 & act2), Jaccard_active = sum(act1 & act2) / max(sum(act1 | act2), 1),
            Corr_raw = stats::cor(raw[[v1]], raw[[v2]]), n_months = nrow(raw))
        }
      }
    }
  }
  ccf <- rbindlist(ccf_rows); best <- rbindlist(best_rows); cooc <- rbindlist(cooc_rows)
  # FDR across states within (Level=state, Pair)
  st <- best[Level == "state"]
  if (nrow(st)) {
    st <- as.data.table(add_fdr(st, by = "Pair"))
    best <- rbindlist(list(best[Level != "state"], st), fill = TRUE)
  }
  save_table(ccf, "cocirculation_ccf.csv")
  save_table(best, "cocirculation_bestlag.csv")
  save_table(cooc, "cocirculation_cooccurrence.csv")

  # Figure: national CCF curves + co-occurrence bars
  nat <- ccf[Level == "national"]
  pairs_sorted <- sort_c(unique(nat$Pair))
  p1 <- ggplot(nat, aes(Lag_months, CCF, colour = factor(Pair, pairs_sorted))) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.8 / 2.13) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.8 / 2.13, linetype = "dashed") +
    geom_line(linewidth = 1.5 / 2.13) + geom_point(size = 3 / .pt * 1.4) +
    scale_colour_manual(values = setNames(TAB10[seq_along(pairs_sorted)], pairs_sorted), name = NULL) +
    labs(x = "Lag (months, + = first virus leads)", y = "Cross-correlation (deseasonalised)") + theme_adv() +
    theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0),
          legend.text = element_text(size = 8), legend.background = element_rect(colour = "grey80", fill = "white"))
  nc <- cooc[Level == "national"]
  p2 <- if (nrow(nc)) ggplot(nc, aes(factor(Pair, Pair), Frac_both_active)) + geom_col(fill = "#2c7fb8", width = 0.8) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Fraction of months both active") +
    theme_adv() + theme(axis.text.x = element_text(angle = 20, hjust = 1)) else plot_spacer()
  caption_note("cocirculation_heatmap.png", "Former title", "National cross-correlation by virus pair", "left")
  caption_note("cocirculation_heatmap.png", "Former title", "National co-activity (both > own median)", "right")
  save_fig(p1 + p2, "cocirculation_heatmap.png", 12, 4.2)
  write_caption_notes("adv_cocirculation.R")
  cat("  Co-circulation done.\n")
}

if (sys.nframe() == 0L) run()
