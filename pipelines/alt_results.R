# ================================================================================
#  ALTERNATIVE-VIEW FIGURES
# ================================================================================
#  Five alternative presentations of the main results, read back from the
#  pipeline outputs (no statistics are recomputed). Written to alt_results/ at 300 dpi:
#    trend01_slope_first_last.png   slope chart, first vs last year, MK-significant states coloured
#    trend03_sen_dotbar.png         sorted Sen-slope dot-and-bar per state
#    spat01_choropleth_facet.png    mean annual incidence small multiples
#    spat02_lisa_facet.png          LISA cluster categories per disease
#    glm02_caterpillar_temp_irr.png per-state temperature IRR (NegBin), sorted
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
quiet_library("ggplot2", "patchwork", "readxl", "sf")

OUT_MAIN <- OUTPUT_DIR
OUT <- path_root("alt_results")
dir.create(OUT, showWarnings = FALSE)

VIRUS_ORDER <- c("Dengue", "Chikungunya", "Zika")
VIRUS_COLOR <- c(Dengue = "#d7263d", Chikungunya = "#1b998b", Zika = "#2e86ab")

theme_alt <- function() {
  theme_mpl(8.5, grid = FALSE) + theme(plot.title = element_text(size = 11, face = "bold"), axis.title = element_text(size = 10),
                                       legend.text = element_text(size = 9))
}
vtitle <- function(v) theme(plot.title = element_text(colour = VIRUS_COLOR[[v]]))
save_alt <- function(p, name, w, h) {
  out <- file.path(OUT, name)
  save_png(p, out, w, h, dpi = 300)
  out
}

# ---- data loaders ----------------------------------------------------------------
annual_incidence <- function() as.data.table(read_excel(file.path(OUT_MAIN, "py_statistical_results.xlsx"), sheet = "Annual_Incidence"))
mann_kendall <- function() as.data.table(read_excel(file.path(OUT_MAIN, "py_statistical_results.xlsx"), sheet = "Mann_Kendall"))
lisa <- function() as.data.table(read_excel(file.path(OUT_MAIN, "geo_spatial_stats.xlsx"), sheet = "LISA_Clusters"))
mean_inc <- function() as.data.table(read_excel(file.path(OUT_MAIN, "geo_spatial_stats.xlsx"), sheet = "Mean_Incidence"))
glm_csv <- function(family, spec) fread(file.path(OUT_MAIN, sprintf("glm_%s_%s.csv", family, spec)))

# matplotlib "symlog" axis (linthresh 10, linscale 1, base 10)
symlog_trans <- function(linthresh = 10) {
  c <- 1 / (1 - 1 / 10)
  fwd <- function(x) ifelse(abs(x) <= linthresh, x * c, sign(x) * linthresh * (c + log10(abs(x) / linthresh)))
  inv <- function(y) ifelse(abs(y) <= linthresh * c, y / c, sign(y) * linthresh * 10^(abs(y) / linthresh - c))
  scales::trans_new("symlog", fwd, inv, breaks = function(lim) { b <- c(0, 10^(1:6)); b[b >= lim[1] & b <= lim[2]] })
}

# ---- figures ----------------------------------------------------------------------
fig_slope_chart <- function(ai, mk) {
  panels <- lapply(VIRUS_ORDER, function(virus) {
    sub <- ai[Virus == virus]
    fy <- min(sub$Year); ly <- max(sub$Year)
    first <- sub[Year == fy, .(Location_Name, y0 = Annual_Inc_100k)]
    last <- sub[Year == ly, .(Location_Name, y1 = Annual_Inc_100k)]
    d <- merge(first, last, by = "Location_Name")
    sig <- mk[Virus == virus & Significant == TRUE]$Location_Name
    d[, hl := Location_Name %in% sig]
    d <- d[order(hl)]
    seg <- rbind(d[, .(Location_Name, hl, x = 0, y = y0)], d[, .(Location_Name, hl, x = 1, y = y1)])
    ggplot(seg, aes(x, y, group = Location_Name)) +
      geom_line(aes(colour = hl, linewidth = hl, alpha = hl)) +
      geom_point(aes(colour = hl, size = hl, alpha = hl)) +
      geom_text(data = d[hl == TRUE], aes(x = 1.02, y = y1, label = Location_Name), inherit.aes = FALSE,
                hjust = 0, size = 8 / .pt, colour = VIRUS_COLOR[[virus]], family = PLOT_FONT) +
      scale_colour_manual(values = c(`TRUE` = VIRUS_COLOR[[virus]], `FALSE` = "#bcbcbc"), guide = "none") +
      scale_linewidth_manual(values = c(`TRUE` = 1.8 / 2.13, `FALSE` = 0.6 / 2.13), guide = "none") +
      scale_size_manual(values = c(`TRUE` = mpl_s_to_size(18), `FALSE` = mpl_s_to_size(6)), guide = "none") +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.6), guide = "none") +
      scale_x_continuous(breaks = c(0, 1), labels = c(fy, ly), expand = c(0, 0)) +
      scale_y_continuous(trans = symlog_trans()) +
      coord_cartesian(xlim = c(-0.1, 1.4), clip = "off") +
      labs(title = virus, x = NULL, y = if (virus == VIRUS_ORDER[1]) "Annual incidence per 100k" else NULL) +
      theme_alt() + vtitle(virus) + theme(panel.grid.major.y = element_line(colour = "#ececec", linewidth = 0.3))
  })
  caption_note("trend01_slope_first_last.png", "Former title",
               "Slope chart \u2014 first vs last year of surveillance, coloured states are MK-significant (p<0.05)")
  save_alt(wrap_plots(panels, nrow = 1), "trend01_slope_first_last.png", 12.5, 5.4)
}

fig_sen_dotbar <- function(mk) {
  key <- c("\u2191 significant" = "#d7263d", "\u2193 significant" = "#2e86ab", ns = "#cccccc")
  panels <- lapply(VIRUS_ORDER, function(virus) {
    sub <- mk[Virus == virus][order(Sen_Slope)]
    sub[, cls := ifelse(!Significant, "ns", ifelse(Sen_Slope > 0, "\u2191 significant", "\u2193 significant"))]
    sub[, Location_Name := factor(Location_Name, Location_Name)]
    # invisible points carrying every class, so all panels share one legend
    all_cls <- data.frame(Sen_Slope = 0, Location_Name = sub$Location_Name[1], cls = names(key))
    ggplot(sub, aes(Sen_Slope, Location_Name, colour = cls)) +
      geom_point(data = all_cls, aes(fill = cls), shape = 21, alpha = 0, show.legend = TRUE) +
      geom_vline(xintercept = 0, linewidth = 0.6 / 2.13) +
      geom_segment(aes(x = 0, xend = Sen_Slope, yend = Location_Name), linewidth = 1.6 / 2.13) +
      geom_point(aes(fill = cls), shape = 21, colour = "white", stroke = 0.6 / 2.13, size = mpl_s_to_size(28)) +
      scale_colour_manual(values = key, limits = names(key), guide = "none") +
      scale_fill_manual(values = key, limits = names(key), name = NULL, drop = FALSE,
                        guide = if (virus == VIRUS_ORDER[1]) guide_legend(override.aes = list(shape = 22, size = 5, colour = NA, alpha = 1))
                                else "none") +
      labs(title = virus, x = "Sen slope (cases per 100k per year)", y = NULL) + theme_alt() + vtitle(virus) +
      theme(axis.text.y = element_text(size = 7.5), panel.grid.major.x = element_line(colour = "#ececec", linewidth = 0.3))
  })
  p <- wrap_plots(panels, nrow = 1, guides = "collect") & theme(legend.position = "top")
  caption_note("trend03_sen_dotbar.png", "Former title", "Mann-Kendall Sen slopes by state \u2014 sorted dot-and-bar")
  save_alt(p, "trend03_sen_dotbar.png", 13, 6.2)
}

fig_choropleth_facet <- function(mean, gdf) {
  panels <- lapply(VIRUS_ORDER, function(virus) {
    g <- merge(gdf, as.data.frame(mean[Virus == virus]), by = "Location_Name", all.x = TRUE)
    ggplot(g) + geom_sf(aes(fill = Mean_Inc), colour = "white", linewidth = 0.4 / 2.13) +
      scale_fill_gradientn(colours = c("#f5f5f5", VIRUS_COLOR[[virus]]), na.value = "#e0e0e0", name = "Mean annual\nincidence /100k",
                           guide = guide_colourbar(barheight = unit(2.2, "in"))) +
      coord_sf(datum = NA) + labs(title = virus) + theme_void(base_family = PLOT_FONT) +
      theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5, colour = VIRUS_COLOR[[virus]]),
            legend.title = element_text(size = 8), legend.text = element_text(size = 8))
  })
  caption_note("spat01_choropleth_facet.png", "Former title", "Choropleth small-multiples \u2014 mean annual incidence per 100k")
  save_alt(wrap_plots(panels, nrow = 1), "spat01_choropleth_facet.png", 13.5, 5.4)
}

fig_lisa_grid <- function(lisa_df, gdf) {
  cat_color <- c(HH = "#b2182b", LL = "#2166ac", HL = "#ef8a62", LH = "#67a9cf", NS = "#dddddd")
  panels <- lapply(VIRUS_ORDER, function(virus) {
    g <- merge(gdf, as.data.frame(lisa_df[Virus == virus, .(Location_Name, LISA_Cat)]), by = "Location_Name", all.x = TRUE)
    g$LISA_Cat[is.na(g$LISA_Cat) | !(g$LISA_Cat %in% names(cat_color))] <- "NS"
    ggplot(g) + geom_sf(aes(fill = LISA_Cat), colour = "white", linewidth = 0.4 / 2.13, show.legend = TRUE) +
      scale_fill_manual(values = cat_color, limits = names(cat_color), drop = FALSE, name = NULL) +
      coord_sf(datum = NA) + labs(title = virus) + theme_void(base_family = PLOT_FONT) +
      theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5, colour = VIRUS_COLOR[[virus]]))
  })
  p <- wrap_plots(panels, nrow = 1, guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 9))
  p <- p & guides(fill = guide_legend(nrow = 1))
  caption_note("spat02_lisa_facet.png", "Former title", "LISA cluster categories \u2014 faceted by disease")
  save_alt(p, "spat02_lisa_facet.png", 13.5, 5.4)
}

fig_caterpillar_per_state <- function(per_state_nb) {
  panels <- lapply(VIRUS_ORDER, function(virus) {
    sub <- per_state_nb[Virus == virus & Term == "Temp_Avg"][order(IRR)]
    sub[, Location_Name := factor(Location_Name, Location_Name)]
    col <- ifelse(sub$Significant, VIRUS_COLOR[[virus]], "#bbbbbb")
    ggplot(sub, aes(IRR, Location_Name)) +
      geom_vline(xintercept = 1, linewidth = 0.6 / 2.13, linetype = "dashed") +
      geom_segment(aes(x = IRR_CI_L, xend = IRR_CI_U, yend = Location_Name), colour = col, linewidth = 1.4 / 2.13) +
      geom_point(shape = 21, fill = col, colour = "white", stroke = 0.6 / 2.13, size = 5 / .pt * 1.5) +
      scale_x_log10() +
      labs(title = virus, x = "Temperature IRR (95% CI) \u2014 NegBin per-state", y = NULL) + theme_alt() + vtitle(virus) +
      theme(axis.text.y = element_text(size = 7.5), panel.grid.major.x = element_line(colour = "#ececec", linewidth = 0.3))
  })
  caption_note("glm02_caterpillar_temp_irr.png", "Former title", "Caterpillar plot \u2014 per-state temperature IRR (NegBin), sorted")
  save_alt(wrap_plots(panels, nrow = 1), "glm02_caterpillar_temp_irr.png", 13, 6.5)
}

main <- function() {
  ai <- annual_incidence(); mk <- mann_kendall(); li <- lisa(); mi <- mean_inc()
  per_state_nb <- glm_csv("negbin", "per_state")
  gdf <- get_brazil_states()
  written <- c(fig_slope_chart(ai, mk), fig_sen_dotbar(mk), fig_choropleth_facet(mi, gdf), fig_lisa_grid(li, gdf),
               fig_caterpillar_per_state(per_state_nb))
  write_caption_notes("alt_results.R")
  for (p in written) cat(sprintf("  - %s\n", basename(p)))
}

if (sys.nframe() == 0L) main()
