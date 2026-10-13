# ================================================================================
#  GEOSPATIAL ANALYSIS — NORMALIZED DB version
# ================================================================================
#  Combined three-disease incidence map (Figure 6), temporal maps, global Moran's I,
#  LISA clusters and Mann-Kendall trend maps. Reads the normalized DB and writes to
#  <repo>/outputs/.
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
quiet_library("ggplot2", "patchwork", "sf")

ALPHA <- 0.05
# Fixed seed so Moran / LISA permutation p-values (and LISA categories) are reproducible
PERM_SEED <- 12345

COLORS <- c(Dengue = "#C0392B", Chikungunya = "#E67E22", Zika = "#2980B9")
CMAPS  <- c(Dengue = "Reds", Chikungunya = "Oranges", Zika = "Blues")
LISA_COLORS <- c(HH = "#D7191C", LL = "#2C7BB6", HL = "#FDAE61", LH = "#ABD9E9", NS = "#EEEEEE")
XLIM <- c(-74, -33); YLIM <- c(-34, 6)
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

# ---- map helpers --------------------------------------------------------------------------------
map_theme <- function(base = 9) {
  theme_void(base_size = base, base_family = PLOT_FONT) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.background = element_rect(fill = "white", colour = NA))
}
map_extent <- function() coord_sf(xlim = XLIM, ylim = YLIM, expand = FALSE, datum = NA)
state_labels <- function(g, size = 4.5) {
  d <- sf::st_drop_geometry(g)
  d <- d[!is.na(d$Abbr), ]
  geom_text(data = d, aes(centroid_x, centroid_y, label = Abbr), size = size / .pt, fontface = "bold",
            family = PLOT_FONT, inherit.aes = FALSE)
}
virus_title <- function(virus, size = 11) theme(plot.title = element_text(colour = COLORS[[virus]], face = "bold", size = size, hjust = 0.5))

# ---- combined map: the three diseases on one map --------------------------------------------------
# One symbol per disease and UF (circle = dengue, triangle = chikungunya, square = Zika), filled by the
# natural-breaks (Jenks, k = 5) class of the population-weighted mean annual incidence, computed per
# disease. UF polygons stay neutral; states too small for three symbols get them offshore with a leader.
MAP_SHAPE <- c(Dengue = 21, Chikungunya = 24, Zika = 22)
MAP_SIZE  <- c(Dengue = 2.6, Chikungunya = 2.4, Zika = 2.4)     # the triangle reads larger at equal size
CLASS_FILL <- viridisLite::viridis(5, direction = -1)            # colour-blind-safe, light -> dark in grayscale
MAP_XLIM <- c(-74.5, -28.0); MAP_YLIM <- c(-34.2, 5.8)
SYMBOL_DX <- 1.05                                                # spacing of the three symbols (degrees)
OFFSHORE <- data.frame(Abbr = c("RN", "PB", "PE", "AL", "SE", "ES", "RJ"),
                       ax = c(-31.4, -31.4, -31.4, -31.4, -31.4, -36.4, -40.2),
                       ay = c(-4.9, -6.6, -8.3, -10.0, -11.7, -20.2, -24.8))
GO_ANCHOR <- c(-50.9, -17.1)   # GO's symbols step south-west so DF (inside GO) keeps its own row

# list(symbols = one row per UF x disease with class and symbol position, bounds = class limits per disease)
combined_map_data <- function(gdf_states, mean_inc) {
  cen <- as.data.table(sf::st_drop_geometry(gdf_states))[, .(Location_Name, Abbr, cx = centroid_x, cy = centroid_y)]
  cen[, `:=`(ax = cx, ay = cy, Offshore = FALSE)]
  cen[OFFSHORE, on = "Abbr", `:=`(ax = i.ax, ay = i.ay, Offshore = TRUE)]
  cen[Abbr == "GO", `:=`(ax = GO_ANCHOR[1], ay = GO_ANCHOR[2])]
  cen[Abbr == "DF", ay := cy + 0.25]
  sym <- list(); bnd <- list()
  for (virus in names(mean_inc)) {
    d <- mean_inc[[virus]][!is.na(Mean_Inc)]
    nb <- natural_breaks(d$Mean_Inc, 5)
    sym[[virus]] <- d[, .(Location_Name, Mean_Inc, Virus = virus, Class = digitize_right(Mean_Inc, nb$bins) + 1L)]
    bnd[[virus]] <- data.table(Virus = virus, Class = seq_len(nb$k), Lower = c(min(d$Mean_Inc), nb$bins[-nb$k]), Upper = nb$bins)
  }
  sym <- merge(cen, rbindlist(sym), by = "Location_Name")
  sym[, Virus := factor(Virus, names(MAP_SHAPE))]
  sym[, `:=`(x = ax + (as.integer(Virus) - 2) * SYMBOL_DX, y = ay)]
  list(symbols = sym[order_c(Location_Name, as.integer(Virus))], bounds = rbindlist(bnd))
}

plot_combined_map <- function(gdf_states, mean_inc, moran_global) {
  md <- combined_map_data(gdf_states, mean_inc)
  sym <- md$symbols
  uf <- unique(sym[, .(Abbr, cx, cy, ax, ay, Offshore)])
  lead <- uf[Offshore == TRUE]
  lab <- uf[, .(Abbr, lx = ifelse(Offshore, ax + SYMBOL_DX * 1.75, ax), ly = ifelse(Offshore, ay, ay - 0.95), hj = ifelse(Offshore, 0, 0.5))]
  f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
  bt <- md$bounds[, .(s = paste(sprintf("%s–%s", f0(Lower), f0(Upper)), collapse = "  |  ")), by = Virus]
  bounds_txt <- paste0("Class bounds, cases per 100,000 (class 1 → 5)\n", paste(sprintf("%s:  %s", bt$Virus, bt$s), collapse = "\n"))
  p <- ggplot() +
    geom_sf(data = gdf_states, fill = "#ececec", colour = "#9a9a9a", linewidth = 0.18) +
    geom_segment(data = lead, aes(cx, cy, xend = ax - SYMBOL_DX * 1.6, yend = ay), colour = "#7a7a7a", linewidth = 0.25) +
    geom_point(data = lead, aes(cx, cy), size = 0.5, colour = "#7a7a7a") +
    geom_point(data = sym, aes(x, y, shape = Virus, fill = factor(Class, 1:5), size = Virus), colour = "#333333", stroke = 0.35) +
    geom_text(data = lab, aes(lx, ly, label = Abbr, hjust = hj), size = 6 / .pt, colour = "#333333", family = PLOT_FONT) +
    annotate("text", x = -46.5, y = -33.6, label = bounds_txt, hjust = 0, vjust = 0, size = 6.4 / .pt, colour = "#333333",
             family = PLOT_FONT, lineheight = 1.1) +
    scale_shape_manual(values = MAP_SHAPE, name = "Disease") +
    scale_size_manual(values = MAP_SIZE, guide = "none") +
    scale_fill_manual(values = setNames(CLASS_FILL, 1:5), drop = FALSE, name = "Mean annual incidence\n(natural-breaks class)",
                      labels = c("1 (lowest)", "2", "3", "4", "5 (highest)")) +
    guides(shape = guide_legend(order = 1, override.aes = list(fill = "#bdbdbd", size = 3)),
           fill = guide_legend(order = 2, override.aes = list(shape = 22, size = 4, colour = "#333333"))) +
    coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE, datum = NA) +
    theme_void(base_size = 9, base_family = PLOT_FONT) +
    theme(plot.background = element_rect(fill = "white", colour = NA), legend.position = "inside",
          legend.position.inside = c(0.015, 0.03), legend.justification = c(0, 0), legend.box = "vertical",
          legend.box.just = "left", legend.spacing.y = unit(6, "pt"), legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7.5), legend.key.size = unit(11, "pt"),
          legend.background = element_rect(fill = "white", colour = NA))
  fname <- "geo_07_combined_map.png"
  for (virus in names(moran_global)) {
    mr <- moran_global[[virus]]
    caption_note(fname, "Global Moran's I (p, 999 permutations)", sprintf("%.3f (p = %.3f)", mr$I, mr$p_sim), virus)
  }
  for (virus in bt$Virus) caption_note(fname, "Class bounds, cases per 100,000 (class 1 -> 5)", bt[Virus == virus]$s, virus)
  caption_note(fname, "Key", paste("Circle = dengue, triangle = chikungunya, square = Zika; fill = natural-breaks (Jenks, k = 5) class of",
                                   "the population-weighted mean annual incidence, computed separately for each disease; symbols for",
                                   "RN, PB, PE, AL, SE, ES and RJ are drawn offshore with a leader line to the state."))
  save_png(p, file.path(OUTPUT_DIR, fname), 8, 7.4, dpi = 200)
  cat(sprintf("    Saved -> %s\n", fname))
}

plot_temporal_evolution <- function(gdf_states, annual) {
  slugs <- list(Dengue = c("dengue", "02"), Chikungunya = c("chikv", "03"), Zika = c("zikv", "04"))
  for (virus in names(annual)) {
    ann <- annual[[virus]]
    sl <- slugs[[virus]]
    years <- sort(unique(ann$Year)); n <- length(years)
    ncols <- 4; nrows <- (n + ncols - 1) %/% ncols
    vmax <- unname(quantile(ann$Annual_Inc_100k, 0.97, na.rm = TRUE))
    long <- do.call(rbind, lapply(years, function(yr) {
      g <- states_with(gdf_states, ann[Year == yr, .(Location_Name, Inc = Annual_Inc_100k)])
      g$Year <- yr
      g
    }))
    labd <- sf::st_drop_geometry(long)
    p <- ggplot(long) + geom_sf(aes(fill = Inc), colour = "#666666", linewidth = 0.25 / 2.13) +
      geom_text(data = labd, aes(centroid_x, centroid_y, label = Abbr), size = 4.5 / .pt, fontface = "bold", family = PLOT_FONT) +
      facet_wrap(~Year, ncol = ncols) + map_extent() +
      scale_fill_gradientn(colours = mpl_colours(CMAPS[[virus]]), limits = c(0, vmax), oob = scales::squish,
                           na.value = "#DDDDDD", name = "Incidence / 100k",
                           guide = guide_colourbar(title.position = "bottom", title.hjust = 0.5, barwidth = unit(6.7, "in"),
                                                   barheight = unit(0.25, "in"), frame.colour = "black")) +
      labs(title = virus) + map_theme(9) + virus_title(virus, 13) +
      theme(strip.text = element_text(face = "bold", size = 9, margin = margin(2, 0, 2, 0)),
            legend.position = "bottom", legend.justification = c(0.9, 0), legend.text = element_text(size = 7),
            legend.title = element_text(size = 9))
    fname <- sprintf("geo_%s_temporal_%s.png", sl[2], sl[1])
    caption_note(fname, "Former title", sprintf("%s — Annual Incidence per 100,000 by State", virus))
    save_png(p, file.path(OUTPUT_DIR, fname), ncols * 4.0, nrows * 3.8, dpi = 130)
    cat(sprintf("    Saved -> %s\n", fname))
  }
}

plot_lisa_clusters <- function(lisa_results) {
  lv <- c("HH", "LL", "HL", "LH", "NS")
  lab <- c(HH = "HH — High cluster", LL = "LL — Low cluster", HL = "HL — High outlier",
           LH = "LH — Low outlier", NS = "NS — Not significant")
  panels <- lapply(names(lisa_results), function(virus) {
    g <- lisa_results[[virus]]
    g$LISA_Cat <- factor(g$LISA_Cat, lv)
    ggplot(g) + geom_sf(aes(fill = LISA_Cat), colour = "#555555", linewidth = 0.35 / 2.13, show.legend = TRUE) + state_labels(g) + map_extent() +
      scale_fill_manual(values = LISA_COLORS, labels = lab, limits = lv, drop = FALSE, name = "LISA Category") +
      labs(title = virus) + map_theme(9) + virus_title(virus, 11)
  })
  p <- wrap_plots(panels, nrow = 1, guides = "collect") & theme(legend.position = "right", legend.justification = c(1, 0),
    legend.text = element_text(size = 8), legend.title = element_text(size = 8.5),
    legend.background = element_rect(colour = "grey80", fill = "white"), legend.margin = margin(4, 6, 4, 6))
  caption_note("geo_05_lisa_clusters.png", "Former title", "LISA Cluster Maps — Local Moran's I (999 permutations, alpha = 0.05)")
  caption_note("geo_05_lisa_clusters.png", "Key (former title)",
               "HH = high-high cluster; LL = low-low cluster; HL/LH = spatial outliers; NS = not significant")
  save_png(p, file.path(OUTPUT_DIR, "geo_05_lisa_clusters.png"), 19, 7.5)
  cat("    Saved -> geo_05_lisa_clusters.png\n")
}

TREND_PALETTE <- data.frame(
  trend = c("increasing", "decreasing", "increasing", "decreasing", "no trend", "no trend"),
  sig = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
  col = c("#B2182B", "#2166AC", "#FDDBC7", "#D1E5F0", "#CCCCCC", "#CCCCCC"),
  label = c("Increasing (sig.)", "Decreasing (sig.)", "Increasing (NS)", "Decreasing (NS)", "No trend", "No trend"))

plot_trend_map <- function(gdf_states, mk_results) {
  labels <- unique(TREND_PALETTE$label)
  cols <- setNames(TREND_PALETTE$col[match(labels, TREND_PALETTE$label)], labels)
  panels <- lapply(names(mk_results), function(virus) {
    g <- states_with(gdf_states, mk_results[[virus]][, .(Location_Name, Trend, Significant)])
    g$Trend[is.na(g$Trend)] <- "no trend"
    g$Significant[is.na(g$Significant)] <- FALSE
    key <- match(paste(tolower(g$Trend), g$Significant), paste(TREND_PALETTE$trend, TREND_PALETTE$sig))
    g$lab <- factor(ifelse(is.na(key), "No trend", TREND_PALETTE$label[key]), labels)
    ggplot(g) + geom_sf(aes(fill = lab), colour = "#555555", linewidth = 0.35 / 2.13, show.legend = TRUE) + state_labels(g) + map_extent() +
      scale_fill_manual(values = cols, limits = labels, drop = FALSE, name = NULL) +
      labs(title = virus) + map_theme(9) + virus_title(virus, 11)
  })
  p <- wrap_plots(panels, nrow = 1, guides = "collect") & theme(legend.position = "right", legend.justification = c(1, 0),
    legend.text = element_text(size = 8), legend.background = element_rect(colour = "grey80", fill = "white"),
    legend.margin = margin(4, 6, 4, 6))
  caption_note("geo_06_trend_map.png", "Former title", "Mann-Kendall Trend Classification by State")
  caption_note("geo_06_trend_map.png", "Key (former title)", "Solid = p < 0.05; pale = not significant")
  save_png(p, file.path(OUTPUT_DIR, "geo_06_trend_map.png"), 19, 7.5)
  cat("    Saved -> geo_06_trend_map.png\n")
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
  plot_lisa_clusters(lisa_results)
  plot_trend_map(gdf_states, mk_results)
  cat("\n  Exporting...\n")
  export_stats(mean_inc, mk_results, lisa_results, moran_global)
  missing_state_years(gdf_states, annual)
  write_caption_notes("arbovirus_geospatial.R")
  cat(sprintf("\n%s\n  Done! -> %s\n%s\n\n", sep, OUTPUT_DIR, sep))
}

if (sys.nframe() == 0L) main()
