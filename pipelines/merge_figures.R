# merge_figures.R
# Combines existing pipeline PNGs into composite article figures (figures/, 200 dpi).
#   Figure 1 — Temporal incidence & regional distribution
#   Figure 2 — Trend analysis (Mann-Kendall + Sen's slope)
#   Figure 3 — Climate-disease correlations (Spearman)
#   Figure 4 — Spatial autocorrelation (Moran scatterplot)
#   Figure 5 — Patient demographics (age + sex, all three diseases)
#   Figure 6 — Combined map: mean incidence class of the three diseases on one map of Brazil
#   Figure 7 — Dengue clinical (ELISA, IgM, hospitalization, serotype)
#   Figure 8 — Chikungunya and Zika clinical
#   Figure S1 — Annual incidence, one combined map per year (FigureS1_<year>) + overview grid
#   Figure S2 — LISA clusters (A) and Mann-Kendall trends (B), each one combined map
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
quiet_library("png", "grid")

OUT_DIR <- OUTPUT_DIR
SAVE_TO <- FIGURES_DIR
DPI <- 200

load_img <- function(path) png::readPNG(path, native = TRUE)
img_h <- function(img) dim(img)[1]
img_w <- function(img) dim(img)[2]
ratio <- function(img) img_w(img) / img_h(img)       # width / height

# Draw an image inside the current viewport, keeping its aspect ratio (centred),
# with an optional bold panel letter just above its top-left corner.
draw_img <- function(img, label = NULL) {
  vw <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE)
  vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  ia <- ratio(img)
  if (vw / vh > ia) { h <- vh; w <- vh * ia } else { w <- vw; h <- vw / ia }
  pushViewport(viewport(width = unit(w, "in"), height = unit(h, "in")))
  grid.raster(img, interpolate = TRUE)
  if (!is.null(label) && nzchar(label))
    grid.text(label, x = unit(-0.01, "npc"), y = unit(1.03, "npc"), hjust = 0, vjust = 0,
              gp = gpar(fontsize = 16, fontface = "bold", fontfamily = PLOT_FONT))
  popViewport()
}

# Composite figure: a grid layout (heights/widths are relative), `cells` = list of
# list(img, row, col, label); `gap` is the spacing between rows/columns (fraction of a cell);
# `extra` is a function drawn last (e.g. side labels).
composite <- function(name, fw, fh, heights, widths, cells, hgap = 0.04, wgap = 0.03, extra = NULL, left = 0.3) {
  top <- 0.4
  path <- file.path(SAVE_TO, name)
  ragg::agg_png(path, width = fw + left, height = fh + top, units = "in", res = DPI, background = "white")
  grid.newpage()
  nr <- length(heights); nc <- length(widths)
  # interleave gaps between rows / columns
  hs <- as.vector(rbind(heights, rep(hgap * mean(heights), nr)))[-(2 * nr)]
  ws <- as.vector(rbind(widths, rep(wgap * mean(widths), nc)))[-(2 * nc)]
  pushViewport(viewport(x = unit(left, "in"), y = unit(0, "npc"), width = unit(fw, "in"), height = unit(fh, "in"),
                        just = c("left", "bottom"),
                        layout = grid.layout(length(hs), length(ws), heights = unit(hs, "null"), widths = unit(ws, "null"))))
  for (cl in cells) {
    rows <- 2 * cl$row - 1; cols <- 2 * range(cl$col) - 1
    pushViewport(viewport(layout.pos.row = rows, layout.pos.col = cols[1]:cols[2]))
    draw_img(cl$img, cl$label)
    popViewport()
  }
  popViewport()
  if (!is.null(extra)) extra(fh, top, left)
  grDevices::dev.off()
  trim_png(path)
  cat(sprintf("  Saved -> %s\n", name))
}

# Crop the white border (matplotlib bbox_inches="tight"), keeping a small pad.
trim_png <- function(path, pad = 12) {
  img <- png::readPNG(path, native = TRUE)
  h <- dim(img)[1]; w <- dim(img)[2]
  m <- matrix(unclass(img), nrow = h, ncol = w, byrow = TRUE)   # nativeRaster data are row-major
  ink <- m != -1L                                                # -1L = opaque white
  rows <- which(rowSums(ink) > 0); cols <- which(colSums(ink) > 0)
  if (!length(rows)) return(invisible(path))
  r <- max(1, min(rows) - pad):min(h, max(rows) + pad)
  cc <- max(1, min(cols) - pad):min(w, max(cols) + pad)
  sub <- m[r, cc, drop = FALSE]
  out <- structure(as.vector(t(sub)), dim = c(nrow(sub), ncol(sub)), class = "nativeRaster", channels = attr(img, "channels"))
  png::writePNG(out, path, dpi = DPI)
  invisible(path)
}
cell <- function(img, row, col, label = NULL) list(img = img, row = row, col = col, label = label)

# ── Figures 1-3: two stacked panels ──────────────────────────────────────────
two_rows <- function(a, b, name) {
  i1 <- load_img(file.path(OUT_DIR, a)); i2 <- load_img(file.path(OUT_DIR, b))
  fw <- 18
  fh <- fw * (1 / ratio(i1) + 1 / ratio(i2)) * 0.5 + 0.3
  composite(name, fw, fh, heights = c(img_h(i1) / img_w(i1), img_h(i2) / img_w(i2)), widths = 1,
            cells = list(cell(i1, 1, 1, "A"), cell(i2, 2, 1, "B")))
}

main <- function() {
  dir.create(SAVE_TO, showWarnings = FALSE)
  two_rows("py_01_annual_incidence.png", "py_03_boxplot_regional.png", "Figure1_incidence_regional.png")
  two_rows("py_04_mannkendall_tau.png", "py_05_sens_slope.png", "Figure2_trend_analysis.png")
  two_rows("py_06_spearman_summary.png", "py_08_lag_profile_region.png", "Figure3_climate_correlations.png")

  # ── Figure 4 ──
  im <- load_img(file.path(OUT_DIR, "py_10_moran_scatterplot.png"))
  composite("Figure4_spatial_autocorrelation.png", 16, 16 / ratio(im), 1, 1, list(cell(im, 1, 1)))

  # ── Figure 5: age violins (row A) + sex (row B) ──
  diseases <- c("Dengue", "Chikungunya", "Zika")
  ages <- lapply(diseases, function(d) load_img(file.path(OUT_DIR, sprintf("%s_04_age_group_violin.png", d))))
  sexes <- lapply(diseases, function(d) load_img(file.path(OUT_DIR, sprintf("%s_03_sex.png", d))))
  fw <- 18
  fh <- fw * (1 / (ratio(ages[[1]]) * 3) + 1 / (ratio(sexes[[1]]) * 3)) + 0.5
  cells <- c(lapply(1:3, function(i) cell(ages[[i]], 1, i, if (i == 1) "A")),
             lapply(1:3, function(i) cell(sexes[[i]], 2, i, if (i == 1) "B")))
  composite("Figure5_demographics.png", fw, fh, heights = c(1 / ratio(ages[[1]]), 1 / ratio(sexes[[1]])), widths = rep(1, 3),
            cells = cells, hgap = 0.05)

  # ── Figure 6: combined three-disease map ──
  combo <- load_img(file.path(OUT_DIR, "geo_07_combined_map.png"))
  composite("Figure6_combined_map.png", 8, 8 / ratio(combo), 1, 1, list(cell(combo, 1, 1)))

  # ── Figure S1: one map per year (the year is tagged on each map) + an overview grid ──
  yrs <- sort(as.integer(sub("^geo_02_temporal_([0-9]{4})\\.png$", "\\1",
                             list.files(OUT_DIR, pattern = "^geo_02_temporal_[0-9]{4}\\.png$"))))
  ty <- lapply(yrs, function(y) load_img(file.path(OUT_DIR, sprintf("geo_02_temporal_%d.png", y))))
  for (i in seq_along(yrs))
    composite(sprintf("FigureS1_%d.png", yrs[i]), 8, 8 / ratio(ty[[i]]), 1, 1, list(cell(ty[[i]], 1, 1)))
  nc <- 4; nr <- ceiling(length(yrs) / nc)
  composite("FigureS1_overview.png", 6 * nc, 6 * nr / ratio(ty[[1]]), rep(1, nr), rep(1, nc),
            lapply(seq_along(yrs), function(i) cell(ty[[i]], (i - 1) %/% nc + 1, (i - 1) %% nc + 1)))

  # ── Figure S2: LISA clusters (A) and Mann-Kendall trends (B) ──
  lisa <- load_img(file.path(OUT_DIR, "geo_05_lisa_clusters.png"))
  trend <- load_img(file.path(OUT_DIR, "geo_06_trend_map.png"))
  composite("FigureS2_lisa_trend.png", 16, 8 / ratio(lisa) + 0.3, 1, rep(1, 2),
            list(cell(lisa, 1, 1, "A"), cell(trend, 1, 2, "B")))

  # ── Figure 7: Dengue clinical 2x2 ──
  d7 <- lapply(c("ELISA", "IgM", "hospitalization", "serotype"), function(s)
    load_img(file.path(OUT_DIR, sprintf("Dengue_%s_%s.png", if (s == "serotype") "03" else "02", s))))
  fw <- 14
  composite("Figure7_dengue_clinical.png", fw, fw / ratio(d7[[1]]), heights = c(1, 1), widths = c(1, 1),
            cells = list(cell(d7[[1]], 1, 1, "A"), cell(d7[[2]], 1, 2, "B"), cell(d7[[3]], 2, 1, "C"), cell(d7[[4]], 2, 2, "D")))

  # ── Figure 8: Chikungunya (3 panels) over Zika (4 panels) on a 12-column grid ──
  chik <- lapply(c("classification", "criteria", "evolution"), function(s) load_img(file.path(OUT_DIR, sprintf("Chikungunya_03_%s.png", s))))
  zika <- lapply(c("classification", "criteria", "evolution", "pregnancy"), function(s) load_img(file.path(OUT_DIR, sprintf("Zika_03_%s.png", s))))
  fw <- 18
  fh <- fw * 2 / (4 * ratio(chik[[1]])) + 0.4
  cells <- c(lapply(1:3, function(i) cell(chik[[i]], 1, c((i - 1) * 4 + 1, i * 4), if (i == 1) "A")),
             lapply(1:4, function(i) cell(zika[[i]], 2, c((i - 1) * 3 + 1, i * 3), if (i == 1) "B")))
  side <- function(fh, top, left) {
    lab <- function(txt, y, col) grid.text(txt, x = unit(left * 0.35, "in"), y = unit(y * fh, "in"), rot = 90,
                                           gp = gpar(fontsize = 11, fontface = "bold", col = col, fontfamily = PLOT_FONT))
    lab("Chikungunya", 0.74, "#E67E22"); lab("Zika", 0.27, "#2980B9")
  }
  composite("Figure8_chikv_zika_clinical.png", fw, fh, heights = c(1, 1), widths = rep(1, 12), cells = cells,
            hgap = 0.06, wgap = 0.0, extra = side, left = 0.7)

  cat(sprintf("\nAll figures saved to: %s\n", SAVE_TO))
}

if (sys.nframe() == 0L) main()
