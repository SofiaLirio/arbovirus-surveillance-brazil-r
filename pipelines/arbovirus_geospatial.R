# ================================================================================
#  GEOSPATIAL ANALYSIS — NORMALIZED DB version
# ================================================================================
#  Global Moran's I, LISA clusters and Mann-Kendall trends per state, and every map
#  figure: mean incidence (Figure 6), annual incidence per year (Figure S1), LISA and
#  trend maps (Figure S2) -- all drawn as one map with the three diseases as symbols.
#  Reads the normalized DB and writes to <repo>/outputs/.
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
source(path_root("pipelines", "data_loader.R"))
quiet_library("ggplot2", "sf")

ALPHA <- 0.05
# Fixed seed so Moran / LISA permutation p-values (and LISA categories) are reproducible
PERM_SEED <- 12345

Q_MAP <- c("1" = "HH", "2" = "LH", "3" = "LL", "4" = "HL")

compute_mean_inc <- function(annual) {
  lapply(annual, function(ann) ann[, .(Mean_Inc = pop_weighted_mean(.SD)), by = Location_Name][order_c(Location_Name)])
}

states_with <- function(gdf_states, df) {
  out <- merge(gdf_states, as.data.frame(df), by = "Location_Name", all.x = TRUE)
  out[order_c(out$Location_Name), ]
}

run_global_moran <- function(gdf_states, mean_inc) {
  res <- list()
  for (virus in names(mean_inc)) {
    g <- states_with(gdf_states, mean_inc[[virus]]); g <- g[!is.na(g$Mean_Inc), ]
    w <- knn_weights(g)
    m <- moran_global(g$Mean_Inc, w, permutations = 999, seed = PERM_SEED)
    res[[virus]] <- list(I = m$I, EI = m$EI, z_sim = m$z_sim, p_sim = m$p_sim, Significant = m$p_sim < ALPHA)
    cat(sprintf("  [%-12s] I=%.4f p=%.4f %s\n", virus, m$I, m$p_sim, if (m$p_sim < ALPHA) "Significant" else "NS"))
  }
  res
}

run_lisa <- function(gdf_states, mean_inc) {
  lisa <- list()
  for (virus in names(mean_inc)) {
    g <- states_with(gdf_states, mean_inc[[virus]]); g <- g[!is.na(g$Mean_Inc), ]
    w <- knn_weights(g)
    lm <- moran_local(g$Mean_Inc, w, permutations = 999, seed = PERM_SEED)
    quad <- unname(Q_MAP[as.character(lm$q)])
    g$LISA_Cat <- ifelse(lm$p_sim < ALPHA, quad, "NS")
    g$Local_I <- lm$Is; g$LISA_p <- lm$p_sim; g$LISA_Quadrant <- lm$q
    fdr <- add_fdr(sf::st_drop_geometry(g), p_col = "LISA_p")   # 27 local tests = one family per virus
    g$q_value <- fdr$q_value; g$Significant_FDR <- fdr$Significant_FDR
    g$LISA_Cat_FDR <- ifelse(g$Significant_FDR, quad, "NS")
    lisa[[virus]] <- g
    counts <- table(g$LISA_Cat); counts <- counts[sort_c(names(counts))]
    cat(sprintf("  [%-12s] %s\n", virus, paste(sprintf("%s=%d", names(counts), as.integer(counts)), collapse = " | ")))
  }
  lisa
}

run_mann_kendall <- function(annual) {
  lapply(annual, function(ann) {
    rows <- list()
    for (uf in sort_c(unique(ann$Location_Name))) {
      series <- ann[Location_Name == uf][order(Year)]$Annual_Inc_100k
      if (length(series) < 4) next
      res <- mk_original_test(series, alpha = ALPHA)
      rows[[length(rows) + 1]] <- data.table(Location_Name = uf, Trend = res$trend, Tau = res$Tau, p_value = res$p,
                                             Sen_Slope = res$slope, Significant = res$p < ALPHA)
    }
    as.data.table(add_fdr(rbindlist(rows)))
  })
}

# ---- symbol maps: the three diseases on one map ---------------------------------------------------
# Every map figure is one map of Brazil with one symbol per disease and UF (circle = dengue,
# triangle = chikungunya, square = Zika), filled by a category. UF polygons stay neutral; states
# too small for three symbols get them offshore with a leader line.
MAP_SHAPE <- c(Dengue = 21, Chikungunya = 24, Zika = 22)
MAP_SIZE  <- c(Dengue = 2.6, Chikungunya = 2.4, Zika = 2.4)     # the triangle reads larger at equal size
MAP_XLIM <- c(-74.5, -28.0); MAP_YLIM <- c(-34.2, 5.8)
SYMBOL_DX <- 1.05                                                # spacing of the three symbols (degrees)
OFFSHORE <- data.frame(Abbr = c("RN", "PB", "PE", "AL", "SE", "ES", "RJ"),
                       ax = c(-31.4, -31.4, -31.4, -31.4, -31.4, -36.4, -40.2),
                       ay = c(-4.9, -6.6, -8.3, -10.0, -11.7, -20.2, -24.8))
GO_ANCHOR <- c(-50.9, -17.1)   # GO's symbols step south-west so DF (inside GO) keeps its own row
MAP_KEY <- paste("Circle = dengue, triangle = chikungunya, square = Zika; symbols for RN, PB, PE, AL, SE, ES and RJ",
                 "are drawn offshore with a leader line to the state.")

# Fill categories. Incidence classes: viridis (colour-blind-safe, light -> dark in grayscale).
# LISA / trend: hues separated in lightness as well (CIELAB L* 25 / 46 / 71 / 84 / 100 and
# 38 / 68 / 100), so the categories stay distinct in grayscale; non-significant = white.
CLASS_FILL <- setNames(viridisLite::viridis(5, direction = -1), 1:5)
CLASS_LABELS <- c("1 (lowest)", "2", "3", "4", "5 (highest)")
LISA_FILL <- c(HH = "#7f0000", LL = "#2171b5", HL = "#fc9272", LH = "#a8d8f8", NS = "#ffffff")
LISA_LABELS <- c(HH = "HH (high-high cluster)", LL = "LL (low-low cluster)", HL = "HL (high-low outlier)",
                 LH = "LH (low-high outlier)", NS = "Not significant")
TREND_FILL <- c(increasing = "#b2182b", decreasing = "#6baed6", "no trend" = "#ffffff")
TREND_LABELS <- c(increasing = "Increasing (p < 0.05)", decreasing = "Decreasing (p < 0.05)", "no trend" = "No significant trend")

# Where each UF's three symbols sit: the polygon centroid, or the offshore / shifted anchor.
symbol_anchors <- function(gdf_states) {
  cen <- as.data.table(sf::st_drop_geometry(gdf_states))[, .(Location_Name, Abbr, cx = centroid_x, cy = centroid_y)]
  cen[, `:=`(ax = cx, ay = cy, Offshore = FALSE)]
  cen[OFFSHORE, on = "Abbr", `:=`(ax = i.ax, ay = i.ay, Offshore = TRUE)]
  cen[Abbr == "GO", `:=`(ax = GO_ANCHOR[1], ay = GO_ANCHOR[2])]
  cen[Abbr == "DF", ay := cy + 0.25]
  cen
}
# d: one row per UF x disease (Location_Name, Virus, ...); adds the symbol position x, y.
place_symbols <- function(anchors, d) {
  sym <- merge(anchors, d, by = "Location_Name")
  sym[, Virus := factor(Virus, names(MAP_SHAPE))]
  sym[, `:=`(x = ax + (as.integer(Virus) - 2) * SYMBOL_DX, y = ay)]
  sym[order_c(Location_Name, as.integer(Virus))]
}
# Natural-breaks (k = 5) classes of `col`, computed per disease over every row of d.
jenks_classes <- function(d, col) {
  sym <- list(); bnd <- list()
  for (virus in unique(d$Virus)) {
    s <- d[Virus == virus & !is.na(get(col))]
    v <- s[[col]]; nb <- natural_breaks(v, 5)
    sym[[virus]] <- s[, Class := digitize_right(v, nb$bins) + 1L]
    bnd[[virus]] <- data.table(Virus = virus, Class = seq_len(nb$k), Lower = c(min(v), nb$bins[-nb$k]), Upper = nb$bins)
  }
  list(rows = rbindlist(sym), bounds = rbindlist(bnd))
}
# "Class bounds ..." note printed on the map, plus one formatted string per disease for the captions.
bounds_note <- function(bounds) {
  f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
  bt <- bounds[, .(s = paste(sprintf("%s–%s", f0(Lower), f0(Upper)), collapse = "  |  ")), by = Virus]
  list(text = paste0("Class bounds, cases per 100,000 (class 1 → 5)\n", paste(sprintf("%s:  %s", bt$Virus, bt$s), collapse = "\n")),
       by_virus = bt)
}

# The map itself. sym: placed symbols with a factor column Fill whose levels are names(fill_values).
# note: text in the bottom-right ocean; tag: panel identifier in the top-left corner (e.g. the year).
symbol_map <- function(gdf_states, sym, fill_values, fill_name, fill_labels, note = NULL, tag = NULL) {
  uf <- unique(sym[, .(Abbr, cx, cy, ax, ay, Offshore)])
  lead <- uf[Offshore == TRUE]
  lab <- uf[, .(Abbr, lx = ifelse(Offshore, ax + SYMBOL_DX * 1.75, ax), ly = ifelse(Offshore, ay, ay - 0.95), hj = ifelse(Offshore, 0, 0.5))]
  p <- ggplot() +
    geom_sf(data = gdf_states, fill = "#ececec", colour = "#9a9a9a", linewidth = 0.18) +
    geom_segment(data = lead, aes(cx, cy, xend = ax - SYMBOL_DX * 1.6, yend = ay), colour = "#7a7a7a", linewidth = 0.25) +
    geom_point(data = lead, aes(cx, cy), size = 0.5, colour = "#7a7a7a") +
    geom_point(data = sym, aes(x, y, shape = Virus, fill = Fill, size = Virus), colour = "#333333", stroke = 0.35,
               show.legend = TRUE) +   # keep legend keys for categories absent from a map (e.g. a quiet year)
    geom_text(data = lab, aes(lx, ly, label = Abbr, hjust = hj), size = 6 / .pt, colour = "#333333", family = PLOT_FONT)
  if (!is.null(note))
    p <- p + annotate("text", x = -46.5, y = -33.6, label = note, hjust = 0, vjust = 0, size = 6.4 / .pt, colour = "#333333",
                      family = PLOT_FONT, lineheight = 1.1)
  if (!is.null(tag))
    p <- p + annotate("text", x = -73.6, y = 5.2, label = tag, hjust = 0, vjust = 1, size = 15 / .pt, fontface = "bold",
                      colour = "#333333", family = PLOT_FONT)
  p + scale_shape_manual(values = MAP_SHAPE, name = "Disease") +
    scale_size_manual(values = MAP_SIZE, guide = "none") +
    scale_fill_manual(values = fill_values, drop = FALSE, name = fill_name, labels = fill_labels) +
    guides(shape = guide_legend(order = 1, override.aes = list(fill = "#bdbdbd", size = 3)),
           fill = guide_legend(order = 2, override.aes = list(shape = 22, size = 4, colour = "#333333"))) +
    coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE, datum = NA) +
    theme_void(base_size = 9, base_family = PLOT_FONT) +
    theme(plot.background = element_rect(fill = "white", colour = NA), legend.position = "inside",
          legend.position.inside = c(0.015, 0.03), legend.justification = c(0, 0), legend.box = "vertical",
          legend.box.just = "left", legend.spacing.y = unit(6, "pt"), legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7.5), legend.key.size = unit(11, "pt"),
          legend.background = element_rect(fill = "white", colour = NA))
}
save_map <- function(p, fname) {
  save_png(p, file.path(OUTPUT_DIR, fname), 8, 7.4, dpi = 200)
  cat(sprintf("    Saved -> %s\n", fname))
}
note_counts <- function(fname, sym, col, labels) {   # e.g. "HH 8 | LL 4 | ..." per disease
  for (virus in levels(sym$Virus)) {
    n <- table(factor(sym[Virus == virus][[col]], names(labels)))
    caption_note(fname, "States per category", paste(sprintf("%s %d", names(n), as.integer(n)), collapse = " | "), virus)
  }
}

# ---- Figure 6: mean annual incidence ----
# list(symbols = one row per UF x disease with class and symbol position, bounds = class limits per disease)
combined_map_data <- function(gdf_states, mean_inc) {
  d <- rbindlist(lapply(names(mean_inc), function(v) mean_inc[[v]][, .(Location_Name, Mean_Inc, Virus = v)]))
  jc <- jenks_classes(d, "Mean_Inc")
  list(symbols = place_symbols(symbol_anchors(gdf_states), jc$rows), bounds = jc$bounds)
}

plot_combined_map <- function(gdf_states, mean_inc, moran_global) {
  md <- combined_map_data(gdf_states, mean_inc)
  sym <- md$symbols[, Fill := factor(Class, 1:5)]
  bn <- bounds_note(md$bounds)
  p <- symbol_map(gdf_states, sym, CLASS_FILL, "Mean annual incidence\n(natural-breaks class)", CLASS_LABELS, note = bn$text)
  fname <- "geo_07_combined_map.png"
  for (virus in names(moran_global)) {
    mr <- moran_global[[virus]]
    caption_note(fname, "Global Moran's I (p, 999 permutations)", sprintf("%.3f (p = %.3f)", mr$I, mr$p_sim), virus)
  }
  for (virus in bn$by_virus$Virus) caption_note(fname, "Class bounds, cases per 100,000 (class 1 -> 5)", bn$by_virus[Virus == virus]$s, virus)
  caption_note(fname, "Key", paste("Fill = natural-breaks (Jenks, k = 5) class of the population-weighted mean annual incidence,",
                                   "computed separately for each disease.", MAP_KEY))
  save_map(p, fname)
}

# ---- Figure S1: annual incidence, one map per year ----
# Classes are computed per disease over all state-years pooled, so a colour means the same
# incidence range in every year.
temporal_map_data <- function(gdf_states, annual) {
  d <- rbindlist(lapply(names(annual), function(v) annual[[v]][, .(Location_Name, Year, Annual_Inc_100k, Virus = v)]))
  jc <- jenks_classes(d, "Annual_Inc_100k")
  list(symbols = place_symbols(symbol_anchors(gdf_states), jc$rows), bounds = jc$bounds)
}

plot_temporal_evolution <- function(gdf_states, annual) {
  md <- temporal_map_data(gdf_states, annual)
  md$symbols[, Fill := factor(Class, 1:5)]
  bn <- bounds_note(md$bounds)
  years <- sort(unique(md$symbols$Year))
  for (yr in years) {
    p <- symbol_map(gdf_states, md$symbols[Year == yr], CLASS_FILL, "Annual incidence\n(natural-breaks class)", CLASS_LABELS,
                    note = bn$text, tag = as.character(yr))
    save_map(p, sprintf("geo_02_temporal_%d.png", yr))
  }
  fname <- "geo_02_temporal_YYYY.png"   # one row per item for the 11 yearly maps (Figure S1)
  for (virus in bn$by_virus$Virus)
    caption_note(fname, "Class bounds, cases per 100,000 (class 1 -> 5; all years pooled)", bn$by_virus[Virus == virus]$s, virus)
  gaps <- CJ(Location_Name = unique(md$symbols$Location_Name), Year = years, Virus = names(MAP_SHAPE))[!md$symbols, on = c("Location_Name", "Year", "Virus")]
  if (nrow(gaps)) caption_note(fname, "No symbol (no SINAN records)",
                               paste(gaps[, sprintf("%s %d (%s)", Location_Name, Year, Virus)], collapse = "; "))
  caption_note(fname, "Key", paste("One map per year; fill = natural-breaks (Jenks, k = 5) class of the annual incidence per",
                                   "100,000, computed per disease over all state-years (2016-2026; 2026 annualised from January-September).",
                                   MAP_KEY))
}

# ---- Figure S2A: LISA clusters ----
plot_lisa_clusters <- function(gdf_states, lisa_results) {
  d <- rbindlist(lapply(names(lisa_results), function(v)
    data.table(Location_Name = lisa_results[[v]]$Location_Name, Virus = v, LISA_Cat = lisa_results[[v]]$LISA_Cat)))
  sym <- place_symbols(symbol_anchors(gdf_states), d)[, Fill := factor(LISA_Cat, names(LISA_FILL))]
  p <- symbol_map(gdf_states, sym, LISA_FILL, "Local Moran's I (LISA)\np < 0.05, 999 permutations", LISA_LABELS)
  fname <- "geo_05_lisa_clusters.png"
  note_counts(fname, sym, "LISA_Cat", LISA_LABELS)
  caption_note(fname, "Key", paste("Fill = LISA category of the population-weighted mean annual incidence (KNN k = 5 weights);",
                                   "HH/LL = high/low clusters, HL/LH = spatial outliers, white = not significant.", MAP_KEY))
  save_map(p, fname)
}

# ---- Figure S2B: Mann-Kendall trends ----
# The Mann-Kendall test labels a state increasing/decreasing only when p < 0.05, so three categories occur.
plot_trend_map <- function(gdf_states, mk_results) {
  d <- rbindlist(lapply(names(mk_results), function(v) mk_results[[v]][, .(Location_Name, Virus = v, Trend = tolower(Trend), Significant)]))
  d[!(Significant %in% TRUE) | !(Trend %in% c("increasing", "decreasing")), Trend := "no trend"]
  sym <- place_symbols(symbol_anchors(gdf_states), d)[, Fill := factor(Trend, names(TREND_FILL))]
  p <- symbol_map(gdf_states, sym, TREND_FILL, "Mann-Kendall trend\n(annual incidence)", TREND_LABELS)
  fname <- "geo_06_trend_map.png"
  note_counts(fname, sym, "Trend", TREND_LABELS)
  caption_note(fname, "Key", paste("Fill = Mann-Kendall trend of the annual incidence per state, 2016-2026 (p < 0.05); white = no",
                                   "significant trend.", MAP_KEY))
  save_map(p, fname)
}

# State-years with no SINAN records: not drawn in the temporal maps (grey polygon) or the py_02
# heatmap (no tile). Logged so the gaps can be reported in the captions.
missing_state_years <- function(gdf_states, annual) {
  rows <- lapply(names(annual), function(virus) {
    ann <- annual[[virus]]
    full <- CJ(Location_Name = sort_c(gdf_states$Location_Name), Year = sort(unique(ann$Year)))
    gap <- full[!ann, on = c("Location_Name", "Year")]
    if (nrow(gap)) gap[, .(Figure = "Figure S1 (temporal maps); supplementary py_02_heatmap_state_year", Virus = virus, Location_Name,
                           Abbr = unname(NAME_TO_ABBR[Location_Name]), Year, Reason = "no SINAN notification records")]
  })
  out <- rbindlist(rows)
  write_csv_py(out, file.path(OUTPUT_DIR, "figure_exclusions_state_years.csv"))
  cat(sprintf("    Saved -> figure_exclusions_state_years.csv (%d state-years without records)\n", nrow(out)))
}

export_stats <- function(mean_inc, mk_results, lisa_results, moran_global) {
  path <- file.path(OUTPUT_DIR, "geo_spatial_stats.xlsx")
  cols <- c("Location_Name", "Region", "Abbr", "Mean_Inc", "LISA_Cat", "Local_I", "LISA_p", "LISA_Quadrant",
            "q_value", "Significant_FDR", "LISA_Cat_FDR")
  sheets <- list(
    Mean_Incidence = rbindlist(lapply(names(mean_inc), function(v) copy(mean_inc[[v]])[, Virus := v])),
    Mann_Kendall = rbindlist(lapply(names(mk_results), function(v) copy(mk_results[[v]])[, Virus := v])),
    LISA_Clusters = rbindlist(lapply(names(lisa_results), function(v) {
      d <- as.data.table(sf::st_drop_geometry(lisa_results[[v]]))[, ..cols]
      d[, Virus := v]
    })),
    Global_Moran = rbindlist(lapply(names(moran_global), function(v) c(list(Virus = v), moran_global[[v]]))))
  write_xlsx_sheets(sheets, path)
  cat("    Saved -> geo_spatial_stats.xlsx\n")
}

main <- function() {
  sep <- strrep("=", 62)
  cat(sprintf("\n%s\n  GEOSPATIAL ANALYSIS — NORMALIZED DB version\n%s\n", sep, sep))
  cat("\n[1/6] State boundaries\n")
  gdf_states <- get_brazil_states()
  cat("\n[2/6] Data loading & annual aggregation\n")
  datasets <- load_pipeline_datasets()
  annual <- compute_annual(datasets); mean_inc <- compute_mean_inc(annual)
  cat("\n[3/6] Global Moran's I\n")
  moran_global <- run_global_moran(gdf_states, mean_inc)
  cat("\n[4/6] Local Moran's I (LISA)\n")
  lisa_results <- run_lisa(gdf_states, mean_inc)
  cat("\n[5/6] Mann-Kendall per state\n")
  mk_results <- run_mann_kendall(annual)
  cat("\n[6/6] Generating maps\n")
  plot_combined_map(gdf_states, mean_inc, moran_global)
  plot_temporal_evolution(gdf_states, annual)
  plot_lisa_clusters(gdf_states, lisa_results)
  plot_trend_map(gdf_states, mk_results)
  cat("\n  Exporting...\n")
  export_stats(mean_inc, mk_results, lisa_results, moran_global)
  missing_state_years(gdf_states, annual)
  write_caption_notes("arbovirus_geospatial.R")
  cat(sprintf("\n%s\n  Done! -> %s\n%s\n\n", sep, OUTPUT_DIR, sep))
}

if (sys.nframe() == 0L) main()
