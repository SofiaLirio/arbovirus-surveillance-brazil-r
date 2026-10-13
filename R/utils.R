# Shared helpers for every script in the project.
# Each script finds the project root and sources this file first (see the header of any script).
# Several helpers reproduce pandas/numpy behaviour exactly, so the rebuilt db/ CSVs are
# byte-identical to the ones the published results came from (reference/expected_hashes.csv).

if (!exists(".root")) stop("source R/utils.R through the script header, which sets .root")
ROOT <- normalizePath(.root, winslash = "/")

# Packages live in the per-user library on Windows (R/install_packages.R creates it).
local({
  ul <- Sys.getenv("R_LIBS_USER")
  if (nzchar(ul) && dir.exists(ul) && !(normalizePath(ul, winslash = "/") %in% normalizePath(.libPaths(), winslash = "/")))
    .libPaths(c(ul, .libPaths()))
})

options(stringsAsFactors = FALSE, warn = 1, datatable.print.class = FALSE, width = 120)
Sys.setenv(TZ = "UTC")

path_root <- function(...) file.path(ROOT, ...)

quiet_library <- function(...) {
  for (p in c(...)) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Run another project script in a fresh R process (mirrors subprocess.run in the Python version).
# Returns the exit status.
RSCRIPT <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run_script <- function(path, args = character()) {
  system2(RSCRIPT, c(shQuote(path), args), stdout = "", stderr = "")
}

# ---------------------------------------------------------------------------------------------
# numpy / pandas numerics
# ---------------------------------------------------------------------------------------------

# numpy round(x, d) = rint(x * 10^d) / 10^d, half-to-even on the scaled double.
# R's round(y) with 0 digits is also half-to-even (IEC 60559), so this matches bit for bit.
py_round <- function(x, digits = 0) {
  if (digits == 0) return(round(x))
  f <- 10^digits
  round(x * f) / f
}

# pandas groupby sum/mean use Kahan-compensated summation in row order and skip NaN.
# `g` is an integer group id (1..G); values are accumulated in the order they appear.
# Vectorised across groups: step p updates the p-th non-missing value of every group at once.
kahan_by <- function(x, g, G = max(g), fun = c("sum", "mean")) {
  fun <- match.arg(fun)
  force(G)   # avalia max(g) com o g ORIGINAL, antes do filtro de NA abaixo (no-op quando G e passado)
  ok <- !is.na(x)
  x <- x[ok]; g <- g[ok]
  s <- numeric(G); comp <- numeric(G); n <- integer(G)
  if (length(x)) {
    pos <- data.table::rowid(g)
    o <- order(pos, method = "radix")
    brk <- split(o, pos[o])
    for (idx in brk) {
      gi <- g[idx]; v <- x[idx]
      y <- v - comp[gi]
      t <- s[gi] + y
      cc <- (t - s[gi]) - y
      cc[is.nan(cc)] <- 0
      comp[gi] <- cc
      s[gi] <- t
    }
    n <- tabulate(g, nbins = G)
  }
  if (fun == "sum") return(list(value = s, n = n))
  m <- s / n
  m[n == 0] <- NaN
  list(value = m, n = n)
}

# pandas Series.sum(min_count=1): NA when there is no non-missing value.
sum_min1 <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

# Benjamini-Hochberg q-values within groups (statsmodels multipletests "fdr_bh"); NA p skipped.
add_fdr <- function(df, p_col = "p_value", by = NULL, alpha = 0.05) {
  df <- as.data.frame(df)
  q <- rep(NA_real_, nrow(df))
  grp <- if (length(by)) interaction(df[by], drop = TRUE, lex.order = TRUE) else factor(rep(1L, nrow(df)))
  for (lv in levels(grp)) {
    i <- which(grp == lv & !is.na(df[[p_col]]))
    if (length(i)) q[i] <- stats::p.adjust(df[[p_col]][i], method = "BH")
  }
  df$q_value <- q
  df$Significant_FDR <- !is.na(q) & q < alpha
  df
}

# ---------------------------------------------------------------------------------------------
# Python string semantics
# ---------------------------------------------------------------------------------------------

# str.strip(): Python whitespace, including NBSP and the ASCII separators.
py_strip <- function(x) gsub("^[\\s\\x{a0}\\x{1c}-\\x{1f}\\x{85}]+|[\\s\\x{a0}\\x{1c}-\\x{1f}\\x{85}]+$", "", x, perl = TRUE)

# unicodedata NFKD + drop combining marks + strip()
strip_accents <- function(x) py_strip(stringi::stri_replace_all_regex(stringi::stri_trans_nfkd(as.character(x)), "\\p{Mn}", ""))

# Code-point (C locale) ordering, as pandas sorts strings.
sort_c <- function(x) sort(x, method = "radix")
order_c <- function(...) order(..., method = "radix")

# ---------------------------------------------------------------------------------------------
# CSV output identical to pandas DataFrame.to_csv(index=False) as run on Windows
# ---------------------------------------------------------------------------------------------

# Python repr() of a float: shortest of 15/16/17 significant digits that round-trips,
# with ".0" appended to integral values. NA -> "".
py_float_str <- function(x) {
  out <- rep("", length(x))
  fin <- which(is.finite(x))
  s <- sprintf("%.15g", x[fin])
  bad <- as.numeric(s) != x[fin]
  if (any(bad)) s[bad] <- sprintf("%.16g", x[fin][bad])
  bad <- as.numeric(s) != x[fin]
  if (any(bad)) s[bad] <- sprintf("%.17g", x[fin][bad])
  plain <- !grepl("[.eE]", s)
  s[plain] <- paste0(s[plain], ".0")
  out[fin] <- s
  inf <- which(is.infinite(x))
  out[inf] <- ifelse(x[inf] > 0, "inf", "-inf")
  out
}

csv_quote <- function(s) {
  need <- !is.na(s) & grepl("[,\"\r\n]", s)
  s[need] <- paste0("\"", gsub("\"", "\"\"", s[need], fixed = TRUE), "\"")
  s[is.na(s)] <- ""
  s
}

py_col_str <- function(v) {
  if (inherits(v, "Date")) return(ifelse(is.na(v), "", format(v, "%Y-%m-%d")))
  if (inherits(v, "POSIXt")) return(ifelse(is.na(v), "", format(v, "%Y-%m-%d %H:%M:%S")))
  if (is.factor(v)) v <- as.character(v)
  if (is.logical(v)) return(ifelse(is.na(v), "", ifelse(v, "True", "False")))
  if (is.integer(v)) return(ifelse(is.na(v), "", as.character(v)))
  if (is.double(v)) return(py_float_str(v))
  csv_quote(enc2utf8(as.character(v)))
}

# Columns that pandas holds as int64 must be integer here (as.integer) to print without ".0".
write_csv_py <- function(df, path) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  header <- paste(csv_quote(enc2utf8(names(df))), collapse = ",")
  body <- if (nrow(df)) do.call(paste, c(lapply(df, py_col_str), sep = ",")) else character()
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(c(header, body), con, sep = "\r\n", useBytes = TRUE)
  invisible(path)
}

read_csv_dt <- function(path, ...) data.table::fread(path, encoding = "UTF-8", na.strings = "", ...)

# ---------------------------------------------------------------------------------------------
# Excel output (one sheet per list element; sheet names cut to Excel's 31 characters)
# ---------------------------------------------------------------------------------------------
write_xlsx_sheets <- function(sheets, path) {
  names(sheets) <- substr(names(sheets), 1, 31)
  sheets <- lapply(sheets, function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE, check.names = FALSE)
    for (j in seq_along(d)) if (is.factor(d[[j]])) d[[j]] <- as.character(d[[j]])
    d
  })
  writexl::write_xlsx(sheets, path)
  invisible(path)
}

# ---------------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------------
PLOT_FONT <- local({
  fam <- "sans"
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    have <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
    if ("DejaVu Sans" %in% have) fam <- "DejaVu Sans"
  }
  fam
})

# Save a ggplot / patchwork / grid object as PNG with the matplotlib figsize (inches) and dpi.
save_png <- function(plot, path, width, height, dpi = 150, bg = "white") {
  ragg::agg_png(path, width = width, height = height, units = "in", res = dpi, background = bg)
  on.exit(grDevices::dev.off())
  if (inherits(plot, c("gg", "patchwork"))) print(plot) else if (inherits(plot, "grob")) grid::grid.draw(plot) else print(plot)
  invisible(path)
}

# ggplot theme close to the matplotlib style the figures were designed with
# (white background, left/bottom spines only, light grid, bold centred titles).
theme_mpl <- function(base_size = 10, grid = TRUE) {
  ggplot2::theme_classic(base_size = base_size, base_family = PLOT_FONT) +
    ggplot2::theme(
      panel.grid.major = if (grid) ggplot2::element_line(colour = "#E7E7E7", linewidth = 0.3) else ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black", size = ggplot2::rel(1)),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      legend.key = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# Composite figure panel of each stage PNG (see pipelines/merge_figures.R); others are supplementary.
FIGURE_OF <- c(
  py_01_annual_incidence = "Figure 1A", py_03_boxplot_regional = "Figure 1B",
  py_04_mannkendall_tau = "Figure 2A", py_05_sens_slope = "Figure 2B",
  py_06_spearman_summary = "Figure 3A", py_08_lag_profile_region = "Figure 3B",
  py_10_moran_scatterplot = "Figure 4",
  Dengue_04_age_group_violin = "Figure 5A", Chikungunya_04_age_group_violin = "Figure 5A", Zika_04_age_group_violin = "Figure 5A",
  Dengue_03_sex = "Figure 5B", Chikungunya_03_sex = "Figure 5B", Zika_03_sex = "Figure 5B",
  geo_07_combined_map = "Figure 6",
  geo_02_temporal_dengue = "Figure S1", geo_03_temporal_chikv = "Figure S1", geo_04_temporal_zikv = "Figure S1",
  geo_05_lisa_clusters = "Figure S2A", geo_06_trend_map = "Figure S2B",
  Dengue_02_ELISA = "Figure 7A", Dengue_02_IgM = "Figure 7B", Dengue_02_hospitalization = "Figure 7C", Dengue_03_serotype = "Figure 7D",
  Chikungunya_03_classification = "Figure 8A", Chikungunya_03_criteria = "Figure 8A", Chikungunya_03_evolution = "Figure 8A",
  Zika_03_classification = "Figure 8B", Zika_03_criteria = "Figure 8B", Zika_03_evolution = "Figure 8B", Zika_03_pregnancy = "Figure 8B",
  glm_fig1_pooled_climate_irr = "GLM Figure 1", glm_fig2_per_state_heatmap = "GLM Figure 2", glm_fig3_lag_profile = "GLM Figure 3",
  glm_fig4_aic_comparison = "GLM Figure 4", glm_fig5_forest_significant = "GLM Figure 5")
figure_of <- function(output_file) {
  key <- sub("\\.png$", "", basename(output_file))
  ifelse(key %in% names(FIGURE_OF), unname(FIGURE_OF[key]), "Supplementary")
}

# Figures carry no in-plot titles; the text they used to hold (former titles, statistics, colour keys)
# is kept for the manuscript captions in outputs/figure_caption_values.csv. caption_note() collects
# rows; write_caption_notes(script) replaces that script's rows in the file.
.caption_rows <- new.env()
.caption_rows$rows <- list()
caption_note <- function(output_file, item, value, panel = "") {
  .caption_rows$rows[[length(.caption_rows$rows) + 1]] <- data.frame(
    Output_file = output_file, Figure = figure_of(output_file), Panel = panel, Item = item, Value = value)
  invisible()
}
write_caption_notes <- function(script) {
  path <- path_root("outputs", "figure_caption_values.csv")
  new <- do.call(rbind, .caption_rows$rows)
  new <- cbind(Script = rep(script, NROW(new)), new)
  old <- if (file.exists(path)) data.table::fread(path, encoding = "UTF-8", colClasses = "character", na.strings = NULL) else NULL
  if (!is.null(old)) old <- as.data.frame(old)[old$Script != script, names(new), drop = FALSE]
  out <- rbind(old, new)
  out <- out[order(out$Script, method = "radix"), ]   # stable: keeps each script's own row order
  write_csv_py(out, path)
  .caption_rows$rows <- list()
  cat(sprintf("  Caption values -> figure_caption_values.csv (%d rows from %s)\n", NROW(new), script))
}

# Colour maps (ColorBrewer data, as matplotlib's named maps): returns f(x in [0,1]) -> colours
mpl_cmap <- function(name) {
  rev_map <- grepl("_r$", name); base <- sub("_r$", "", name)
  cols <- switch(base,
    RdBu = RColorBrewer::brewer.pal(11, "RdBu"),
    viridis = viridisLite::viridis(256),
    RColorBrewer::brewer.pal(9, base))
  if (rev_map) cols <- rev(cols)
  ramp <- grDevices::colorRamp(cols, space = "rgb")
  function(x) {
    x <- pmin(pmax(x, 0), 1)
    m <- ramp(x)
    grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  }
}
mpl_colours <- function(name, n = 256) mpl_cmap(name)(seq(0, 1, length.out = n))

# matplotlib scatter `s` (points^2) -> ggplot point size (mm)
mpl_s_to_size <- function(s) sqrt(s) / ggplot2::.pt

# matplotlib "tab10" colour cycle
TAB10 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
TAB20 <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
           "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")

# Print a data.frame like pandas' to_string (for console reports only).
print_df <- function(df, ...) print(as.data.frame(df), row.names = FALSE, ...)

# Command-line flags: --name / --name=value / --name value
cli_args <- function() commandArgs(trailingOnly = TRUE)
cli_flag <- function(name, args = cli_args()) any(args == paste0("--", name))
cli_value <- function(name, default = NULL, args = cli_args()) {
  key <- paste0("--", name)
  i <- which(args == key)
  if (length(i) && i[1] < length(args)) return(args[i[1] + 1])
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", key, "="), "", hit[1]))
  default
}
