# ================================================================================
#  ARBOVIRUS ANALYSIS PIPELINE (NORMALIZED DB version)
# ================================================================================
#  Descriptives, Mann-Kendall trends, Spearman climate correlations and global
#  Moran's I. Reads the normalized DB via load_pipeline_datasets() and writes to
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
quiet_library("ggplot2", "patchwork", "classInt")

ALPHA <- 0.05
# Fixed seed so Moran's I permutation p-values are reproducible run-to-run
PERM_SEED <- 12345

COLORS <- c(Dengue = "#C0392B", Chikungunya = "#E67E22", Zika = "#2980B9")
REGION_COLORS <- c(Norte = "#1B7837", Nordeste = "#762A83", Sudeste = "#E08214", Sul = "#4393C3", "Centro-Oeste" = "#D6604D")

STATE_CENTROIDS <- data.frame(
  Location_Name = c("Rondonia", "Acre", "Amazonas", "Roraima", "Para", "Amapa", "Tocantins", "Maranhao", "Piaui", "Ceara",
                    "Rio Grande do Norte", "Paraiba", "Pernambuco", "Alagoas", "Sergipe", "Bahia", "Minas Gerais",
                    "Espirito Santo", "Rio de Janeiro", "Sao Paulo", "Parana", "Santa Catarina", "Rio Grande do Sul",
                    "Mato Grosso do Sul", "Mato Grosso", "Goias", "Distrito Federal"),
  lon = c(-63.0, -70.5, -64.0, -61.4, -52.0, -51.0, -48.3, -44.5, -42.8, -39.3, -36.7, -36.7, -37.8, -36.6, -37.1, -41.7,
          -44.7, -40.3, -43.2, -48.5, -51.6, -50.0, -53.1, -54.4, -56.1, -49.6, -47.9),
  lat = c(-10.8, -9.0, -3.5, 2.0, -3.4, 1.0, -9.0, -4.9, -7.7, -5.5, -5.8, -7.1, -8.4, -9.7, -10.6, -12.5, -18.5, -20.0,
          -22.9, -22.0, -24.9, -27.2, -30.0, -20.5, -13.0, -15.8, -15.8))

virus_title <- function(virus, size = 11) {
  theme(plot.title = element_text(colour = COLORS[[virus]], face = "bold", size = size, hjust = 0.5))
}

# ==============================================================================
# 2. DESCRIPTIVE — Annual Incidence
# ==============================================================================
plot_descriptive <- function(annual, datasets) {
  panels <- lapply(names(annual), function(virus) {
    ann <- annual[[virus]]; c <- COLORS[[virus]]
    nat <- ann[, .(Inc = pop_weighted_mean(.SD)), by = Year][order(Year)]
    reg <- ann[, .(Inc = pop_weighted_mean(.SD)), by = .(Region, Year)][order_c(Region, Year)]
    lv <- c("Brazil", sort_c(unique(reg$Region)))
    p <- ggplot() +
      geom_area(data = nat, aes(Year, Inc), fill = c, alpha = 0.12) +
      geom_line(data = reg, aes(Year, Inc, colour = factor(Region, lv)), linetype = "dashed", linewidth = 0.4, alpha = 0.75) +
      geom_line(data = nat, aes(Year, Inc, colour = factor("Brazil", lv)), linewidth = 1) +
      geom_point(data = nat, aes(Year, Inc, colour = factor("Brazil", lv)), size = 1.6) +
      scale_colour_manual(values = c(Brazil = c, REGION_COLORS), breaks = lv, name = NULL, drop = FALSE) +
      scale_x_continuous(breaks = scales::breaks_pretty(), labels = function(x) format(x, scientific = FALSE)) +
      labs(title = virus, x = "Year", y = "Incidence per 100,000 inhabitants") +
      theme_mpl(10) + virus_title(virus, 12) +
      theme(axis.title = element_text(size = 9))
    if (virus == names(annual)[1])
      p + theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
                legend.text = element_text(size = 7), legend.background = element_rect(colour = "grey80", fill = "white"),
                legend.key.height = unit(9, "pt"))
    else p + theme(legend.position = "none")
  })
  caption_note("py_01_annual_incidence.png", "Former title", "Annual Arbovirus Incidence in Brazil (per 100,000 inhabitants)")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_01_annual_incidence.png"), 18, 5.5)
  cat("    Saved -> py_01_annual_incidence.png\n")

  panels <- lapply(names(annual), function(virus) {
    ann <- annual[[virus]]
    rm <- ann[, .(m = mean(Annual_Inc_100k, na.rm = TRUE)), by = Location_Name][order(-m)]
    d <- ann[, .(Inc = mean(Annual_Inc_100k)), by = .(Location_Name, Year)]
    d[, Location_Name := factor(Location_Name, rev(rm$Location_Name))]
    ggplot(d, aes(factor(Year), Location_Name, fill = Inc)) +
      geom_tile(colour = "white", linewidth = 0.25) +
      scale_fill_gradientn(colours = mpl_colours("YlOrRd"), na.value = "white", name = "Inc./100k",
                           guide = guide_colourbar(barheight = unit(0.7 * 7, "in") / 1.6)) +
      labs(title = virus, x = "Year", y = NULL) + theme_mpl(10, grid = FALSE) + virus_title(virus) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), axis.text.y = element_text(size = 7),
            axis.line = element_blank(), axis.title.x = element_text(size = 9))
  })
  caption_note("py_02_heatmap_state_year.png", "Former title", "Annual Incidence by State (per 100,000 inhabitants)")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_02_heatmap_state_year.png"), 22, 9)
  cat("    Saved -> py_02_heatmap_state_year.png\n")

  order_r <- c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste")
  y_max <- max(vapply(datasets, function(df) unname(quantile(df$Monthly_Inc_100k, 0.99, na.rm = TRUE)), 0))
  panels <- lapply(names(datasets), function(virus) {
    df <- datasets[[virus]][!is.na(Monthly_Inc_100k) & !is.na(Region)]
    ggplot(df, aes(factor(Region, order_r), Monthly_Inc_100k, fill = Region)) +
      geom_boxplot(outlier.shape = 16, outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.4, width = 0.8) +
      scale_fill_manual(values = REGION_COLORS, guide = "none") +
      coord_cartesian(ylim = c(0, y_max * 1.05)) +
      labs(title = virus, x = NULL, y = "Monthly Incidence / 100k") + theme_mpl(10) + virus_title(virus) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8), axis.title.y = element_text(size = 9))
  })
  caption_note("py_03_boxplot_regional.png", "Former title", "Monthly Incidence Distribution by Region")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_03_boxplot_regional.png"), 18, 5.5)
  cat("    Saved -> py_03_boxplot_regional.png\n")

  all_data <- rbindlist(lapply(names(datasets), function(v) datasets[[v]][!is.na(Monthly_Inc_100k)][, Virus := v]))
  p <- ggplot(all_data, aes(factor(Virus, names(COLORS)), Monthly_Inc_100k, fill = Virus)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.4, width = 0.8) +
    scale_fill_manual(values = COLORS, guide = "none") +
    labs(x = NULL, y = "Monthly Incidence / 100k") +
    theme_mpl(10) + theme(axis.text.x = element_text(size = 11))
  caption_note("py_03b_boxplot_disease_comparison.png", "Former title", "Monthly Incidence Comparison by Disease (per 100,000 inhabitants)")
  save_png(p, file.path(OUTPUT_DIR, "py_03b_boxplot_disease_comparison.png"), 8, 5.5)
  cat("    Saved -> py_03b_boxplot_disease_comparison.png\n")
}

# ==============================================================================
# 3. MANN-KENDALL
# ==============================================================================
run_mann_kendall <- function(annual) {
  out <- list()
  for (virus in names(annual)) {
    ann <- annual[[virus]]
    rows <- list()
    for (uf in sort_c(unique(ann$Location_Name))) {
      grp <- ann[Location_Name == uf][order(Year)]
      series <- grp$Annual_Inc_100k
      if (length(series) < 4) next
      res <- mk_original_test(series, alpha = ALPHA)
      rows[[length(rows) + 1]] <- data.table(Virus = virus, Location_Name = uf, Region = grp$Region[1],
        Trend = res$trend, Tau = res$Tau, p_value = res$p, S = res$s, Sen_Slope = res$slope,
        Sen_Intercept = res$intercept, Significant = res$p < ALPHA, n_years = length(series))
    }
    df_res <- as.data.table(add_fdr(rbindlist(rows)))
    out[[virus]] <- df_res
    cat(sprintf("  [%-12s] Significant: %d/%d (FDR: %d) | Increasing: %d | Decreasing: %d\n", virus,
                sum(df_res$Significant), nrow(df_res), sum(df_res$Significant_FDR),
                sum(df_res$Trend == "increasing"), sum(df_res$Trend == "decreasing")))
  }
  out
}

# Fill scale with the two legend patches the Python figures show ("p < 0.05" / "Not significant");
# "neg" (negative correlation, blue) is drawn but not listed in the legend.
SIG_LAB <- paste("p <", ALPHA)
sig_legend <- function(c) {
  list(scale_fill_manual(values = setNames(c(c, "#CCCCCC", "#2980B9"), c(SIG_LAB, "Not significant", "neg")),
                         breaks = c(SIG_LAB, "Not significant"), limits = c(SIG_LAB, "Not significant", "neg"), name = NULL))
}
# legend in the lower-right corner of the panel (add after the theme)
legend_br <- theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0),
                   legend.text = element_text(size = 8), legend.background = element_rect(colour = "grey80", fill = "white"))

plot_mann_kendall <- function(mk_results) {
  panels <- lapply(names(mk_results), function(virus) {
    df <- mk_results[[virus]][order(Tau)]; c <- COLORS[[virus]]
    df[, fillc := ifelse(Significant, SIG_LAB, "Not significant")]
    caption_note("py_04_mannkendall_tau.png", "Significant states (p < 0.05)", sprintf("%d/%d", sum(df$Significant), nrow(df)), virus)
    ggplot(df, aes(Tau, factor(Location_Name, Location_Name), fill = fillc)) +
      geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, alpha = 0.7) +
      coord_cartesian(xlim = c(-1.1, 1.1)) + sig_legend(c) +
      labs(title = virus, x = "Kendall's Tau", y = NULL) +
      theme_mpl(10) + virus_title(virus, 10) + theme(axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 9)) +
      legend_br
  })
  caption_note("py_04_mannkendall_tau.png", "Former title", "Mann-Kendall Test — Kendall's Tau by State (coloured bars = p < 0.05)")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_04_mannkendall_tau.png"), 18, 9)
  cat("    Saved -> py_04_mannkendall_tau.png\n")

  all_mk <- rbindlist(mk_results)
  ord <- all_mk[, .(m = mean(Sen_Slope)), by = Location_Name][order(m)]$Location_Name
  all_mk[, Location_Name := factor(Location_Name, ord)]
  all_mk[, Virus := factor(Virus, names(COLORS))]
  p <- ggplot(all_mk, aes(Location_Name, Sen_Slope, fill = Virus)) +
    geom_col(position = position_dodge(width = 0.84), width = 0.84, alpha = 0.85, colour = "white", linewidth = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    scale_fill_manual(values = COLORS, name = NULL) +
    labs(x = NULL, y = "Sen's Slope (cases/100k per year)") +
    theme_mpl(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5), axis.title.y = element_text(size = 9),
                          legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
                          legend.background = element_rect(colour = "grey80", fill = "white"))
  caption_note("py_05_sens_slope.png", "Former title", "Sen's Slope by State — Rate of Change in Incidence")
  save_png(p, file.path(OUTPUT_DIR, "py_05_sens_slope.png"), 13, 6)
  cat("    Saved -> py_05_sens_slope.png\n")
}

# ==============================================================================
# 4. SPEARMAN
# ==============================================================================
LAG_MAP <- data.frame(Variable = c("Temp", "Temp", "Temp", "Rain", "Rain", "Rain"), Lag = c(0L, 1L, 2L, 0L, 1L, 2L),
                      col = c("Temp_Avg", "Temp_Lag1", "Temp_Lag2", "Rain_mm", "Rain_Lag1", "Rain_Lag2"))

run_spearman <- function(datasets) {
  out <- list()
  for (virus in names(datasets)) {
    df <- datasets[[virus]]
    rows <- list()
    for (uf in sort_c(unique(df$Location_Name))) {
      gc <- df[Location_Name == uf & !is.na(Monthly_Inc_100k)]
      for (r in seq_len(nrow(LAG_MAP))) {
        col <- LAG_MAP$col[r]
        sub <- gc[!is.na(gc[[col]])]
        if (nrow(sub) < 8) next
        st <- spearman_test(sub$Monthly_Inc_100k, sub[[col]])
        rows[[length(rows) + 1]] <- data.table(Virus = virus, Location_Name = uf, Region = gc$Region[1],
          Variable = LAG_MAP$Variable[r], Lag = LAG_MAP$Lag[r], Rho = st$rho, p_value = st$p,
          Significant = !is.na(st$p) & st$p < ALPHA, n = nrow(sub))
      }
    }
    # One test family per variable (27 states x 3 lags)
    df_res <- as.data.table(add_fdr(rbindlist(rows), by = "Variable"))
    out[[virus]] <- df_res
    cat(sprintf("  [%-12s] %d pairs tested | %d significant (alpha = %s) | %d after FDR\n", virus, nrow(df_res),
                sum(df_res$Significant), ALPHA, sum(df_res$Significant_FDR)))
  }
  out
}

rdbu_scale <- function(name = "Rho", ...) {
  scale_fill_gradientn(colours = mpl_colours("RdBu_r"), limits = c(-1, 1), name = name, ...)
}

best_abs_row <- function(sub) sub[, .SD[which.max(abs(Rho))], by = Location_Name]

plot_spearman <- function(spearman) {
  panels <- lapply(names(spearman), function(virus) {
    df <- spearman[[virus]]
    s <- df[, .(Mean_Rho = mean(Rho), Prop_Sig = mean(Significant)), by = .(Variable, Lag)]
    s[, label := sprintf("%.2f%s", Mean_Rho, ifelse(Prop_Sig >= 0.5, "*", ""))]
    s[, Var := factor(ifelse(Variable == "Rain", "Rainfall\n(mm)", "Temperature\n(°C)"),
                      c("Temperature\n(°C)", "Rainfall\n(mm)"))]
    ggplot(s, aes(factor(paste("Lag", Lag)), Var, fill = Mean_Rho)) + geom_tile() +
      geom_text(aes(label = label, colour = abs(Mean_Rho) > 0.55), fontface = "bold", size = 9.5 / .pt, family = PLOT_FONT) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
      rdbu_scale(guide = guide_colourbar(barheight = unit(1.4, "in"))) +
      scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
      labs(title = virus, x = NULL, y = NULL) + theme_mpl(10, grid = FALSE) + virus_title(virus, 10) +
      theme(axis.line = element_blank())
  })
  caption_note("py_06_spearman_summary.png", "Former title", "Spearman Correlation: Monthly Incidence × Climate")
  caption_note("py_06_spearman_summary.png", "Key (former title)", "Cells = mean Rho across states; * = ≥ 50% of states with p < 0.05")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_06_spearman_summary.png"), 14, 4.5)
  cat("    Saved -> py_06_spearman_summary.png\n")

  panels <- list()
  for (var in c("Temp", "Rain")) for (virus in names(spearman)) {
    best <- best_abs_row(spearman[[virus]][Variable == var])[order(Rho)]
    best[, `:=`(fillc = ifelse(Rho > 0, "#C0392B", "#2980B9"), al = ifelse(Significant, 1, 0.28),
                xpos = Rho + ifelse(Rho >= 0, 0.03, -0.03), hj = ifelse(Rho >= 0, 0, 1))]
    best[, Location_Name := factor(Location_Name, Location_Name)]
    panels[[length(panels) + 1]] <- ggplot(best, aes(Rho, Location_Name)) +
      geom_col(aes(fill = fillc, alpha = al), width = 0.72, colour = "white", linewidth = 0.15) +
      geom_text(aes(x = xpos, label = paste0("L", Lag), hjust = hj), size = 5.5 / .pt, family = PLOT_FONT) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
      scale_fill_identity() + scale_alpha_identity() + coord_cartesian(xlim = c(-1.15, 1.15)) +
      labs(title = virus, x = sprintf("Spearman's Rho, %s", if (var == "Temp") "temperature" else "rainfall"), y = NULL) +
      theme_mpl(10) + virus_title(virus, 9) + theme(axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 8))
  }
  caption_note("py_07_spearman_state_bestlag.png", "Former title", "Spearman's Rho by State — Lag with Highest Absolute Correlation")
  caption_note("py_07_spearman_state_bestlag.png", "Key (former title)", "Opaque bar = p < 0.05; L0/L1/L2 = optimal lag; top row temperature, bottom row rainfall")
  save_png(wrap_plots(panels, nrow = 2), file.path(OUTPUT_DIR, "py_07_spearman_state_bestlag.png"), 20, 12)
  cat("    Saved -> py_07_spearman_state_bestlag.png\n")

  all_sp <- rbindlist(spearman)
  styles <- c(Dengue = "solid", Chikungunya = "dashed", Zika = "dotdash")
  panels <- lapply(c("Temp", "Rain"), function(var) {
    reg_lag <- all_sp[Variable == var, .(Rho = mean(Rho)), by = .(Region, Lag, Virus)]
    reg_lag[, key := sprintf("%s (%s)", Region, Virus)]
    keys <- unique(reg_lag[order_c(Region, Virus)]$key)
    kr <- sub(" \\(.*$", "", keys); kv <- sub("^.*\\((.*)\\)$", "\\1", keys)
    ggplot(reg_lag, aes(Lag, Rho, colour = key, linetype = key, group = key)) +
      geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3) +
      geom_line(linewidth = 0.55, alpha = 0.82) + geom_point(size = 1.4, alpha = 0.82) +
      scale_colour_manual(values = setNames(unname(REGION_COLORS[kr]), keys), breaks = keys[1:12], name = NULL) +
      scale_linetype_manual(values = setNames(unname(styles[kv]), keys), breaks = keys[1:12], name = NULL) +
      scale_x_continuous(breaks = 0:2, labels = c("Lag 0", "Lag 1", "Lag 2")) + coord_cartesian(ylim = c(-0.65, 0.65)) +
      labs(x = NULL, y = sprintf("Mean Spearman's Rho, %s", if (var == "Temp") "temperature" else "rainfall")) +
      theme_mpl(10) + theme(axis.title.y = element_text(size = 9))
  })
  p <- wrap_plots(panels, nrow = 1, guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 6.5),
    legend.background = element_rect(colour = "grey80", fill = "white"), legend.key.width = unit(22, "pt"))
  p <- p & guides(colour = guide_legend(ncol = 6, byrow = TRUE), linetype = guide_legend(ncol = 6, byrow = TRUE))
  caption_note("py_08_lag_profile_region.png", "Former title", "Lag Profile: Mean Rho by Region and Time Lag")
  save_png(p, file.path(OUTPUT_DIR, "py_08_lag_profile_region.png"), 14, 5)
  cat("    Saved -> py_08_lag_profile_region.png\n")
}

plot_spearman_rho_only <- function(spearman) {
  lag_order <- c("Temp Lag0", "Temp Lag1", "Temp Lag2", "Rain Lag0", "Rain Lag1", "Rain Lag2")
  panels <- lapply(names(spearman), function(virus) {
    df <- copy(spearman[[virus]])
    df[, Col := paste0(Variable, " Lag", Lag)]
    pv <- df[, .(Rho = mean(Rho)), by = .(Location_Name, Col)]
    ord <- pv[, .(m = mean(Rho, na.rm = TRUE)), by = Location_Name][order(-m)]$Location_Name
    pv[, `:=`(Location_Name = factor(Location_Name, rev(ord)), Col = factor(Col, intersect(lag_order, unique(Col))))]
    ggplot(pv, aes(Col, Location_Name, fill = Rho)) + geom_tile(colour = "white", linewidth = 0.25) +
      geom_text(aes(label = sprintf("%.2f", Rho), colour = abs(Rho) > 0.6), size = 6.5 / .pt, family = PLOT_FONT) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
      rdbu_scale(guide = guide_colourbar(barheight = unit(4.5, "in"))) +
      labs(title = virus, x = NULL, y = NULL) + theme_mpl(10, grid = FALSE) + virus_title(virus) +
      theme(axis.line = element_blank(), axis.text.y = element_text(size = 6.5), axis.text.x = element_text(size = 8, angle = 35, hjust = 1))
  })
  caption_note("py_06b_spearman_rho_only.png", "Former title", "Spearman's Rho by State — All Variables and Lags")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_06b_spearman_rho_only.png"), 22, 10)
  cat("    Saved -> py_06b_spearman_rho_only.png\n")
}

plot_spearman_rain_significant <- function(spearman) {
  sig_data <- lapply(spearman, function(df) {
    rain <- df[Variable == "Rain"]
    sig_only <- rain[Significant == TRUE]
    if (!nrow(sig_only)) return(data.table())
    best_abs_row(sig_only)[order(Rho)]
  })
  if (all(vapply(sig_data, nrow, 0L) == 0)) { cat("    (no significant rain-incidence correlations found)\n"); return(invisible()) }
  panels <- lapply(names(sig_data), function(virus) {
    df <- sig_data[[virus]]; c <- COLORS[[virus]]
    caption_note("py_06d_rain_significant_states.png", "Significant states (p < 0.05)", sprintf("%d", nrow(df)), virus)
    if (!nrow(df)) return(ggplot() + annotate("text", x = 0, y = 0, label = "No state with p < 0.05", size = 9 / .pt, family = PLOT_FONT) +
                            theme_void() + labs(title = virus) + virus_title(virus, 10))
    df[, `:=`(fillc = ifelse(Rho > 0, c, "#2980B9"), xpos = Rho + ifelse(Rho >= 0, 0.03, -0.03), hj = ifelse(Rho >= 0, 0, 1))]
    df[, Location_Name := factor(Location_Name, Location_Name)]
    ggplot(df, aes(Rho, Location_Name)) + geom_col(aes(fill = fillc), width = 0.72, colour = "white", linewidth = 0.15) +
      geom_text(aes(x = xpos, label = paste0("L", Lag), hjust = hj), size = 7 / .pt, family = PLOT_FONT) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) + scale_fill_identity() +
      coord_cartesian(xlim = c(-1.15, 1.15)) +
      labs(title = virus, x = "Spearman's Rho", y = NULL) +
      theme_mpl(10) + virus_title(virus, 10) + theme(axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 9))
  })
  caption_note("py_06d_rain_significant_states.png", "Former title", "States with Significant Spearman Correlation: Rainfall × Incidence")
  caption_note("py_06d_rain_significant_states.png", "Key (former title)", "Best significant lag shown (L0/L1/L2); p < 0.05")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_06d_rain_significant_states.png"), 18, 8)
  cat("    Saved -> py_06d_rain_significant_states.png\n")
}

run_spearman_rain_lag2 <- function(datasets) {
  out <- list()
  for (virus in names(datasets)) {
    df <- datasets[[virus]]
    rows <- list()
    for (uf in sort_c(unique(df$Location_Name))) {
      grp <- df[Location_Name == uf][order(Date)]
      rain <- grp$Rain_mm; inc <- grp$Monthly_Inc_100k
      if (length(rain) < 10) next
      rain_t <- head(rain, -2); inc_tp2 <- tail(inc, -2)
      valid <- is.finite(rain_t) & is.finite(inc_tp2)
      if (sum(valid) < 8) next
      st <- spearman_test(inc_tp2[valid], rain_t[valid])
      rows[[length(rows) + 1]] <- data.table(Virus = virus, Location_Name = uf, Region = grp$Region[1], Rho = st$rho,
        p_value = st$p, Significant = !is.na(st$p) & st$p < ALPHA, n = as.integer(sum(valid)))
    }
    out[[virus]] <- as.data.table(add_fdr(rbindlist(rows)))
    cat(sprintf("  [%-12s] Rain Lag-2 | %d states | %d significant (alpha = %s) | %d after FDR\n", virus, length(rows),
                sum(out[[virus]]$Significant), ALPHA, sum(out[[virus]]$Significant_FDR)))
  }
  out
}

plot_spearman_rain_lag2 <- function(rain_lag2) {
  panels <- lapply(names(rain_lag2), function(virus) {
    df <- rain_lag2[[virus]][order(Rho)]; c <- COLORS[[virus]]
    df[, `:=`(fillc = ifelse(Rho > 0, SIG_LAB, "neg"), al = ifelse(Significant, 1, 0.28))]
    df[, Location_Name := factor(Location_Name, Location_Name)]
    caption_note("py_06c_spearman_rain_lag2.png", "Significant states (p < 0.05)", sprintf("%d", sum(df$Significant)), virus)
    ggplot(df, aes(Rho, Location_Name)) +
      geom_col(aes(fill = fillc, alpha = al), width = 0.75, colour = "white", linewidth = 0.15) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) + scale_alpha_identity() +
      sig_legend(c) + coord_cartesian(xlim = c(-1.1, 1.1)) +
      labs(title = virus, x = "Spearman's Rho", y = NULL) +
      theme_mpl(10) + virus_title(virus, 10) + theme(axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 9)) +
      legend_br
  })
  caption_note("py_06c_spearman_rain_lag2.png", "Former title", "Spearman Correlation: Rainfall (month t) × Incidence (month t+2)")
  caption_note("py_06c_spearman_rain_lag2.png", "Key (former title)", "Opaque bar = p < 0.05")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_06c_spearman_rain_lag2.png"), 18, 9)
  cat("    Saved -> py_06c_spearman_rain_lag2.png\n")
}

# ==============================================================================
# 5. SPATIAL
# ==============================================================================
build_spatial <- function(annual) {
  lapply(annual, function(ann) {
    mean_inc <- ann[, .(Mean_Inc = pop_weighted_mean(.SD)), by = Location_Name]
    g <- merge(as.data.table(STATE_CENTROIDS), mean_inc, by = "Location_Name", all.x = TRUE, sort = FALSE)
    g <- g[match(STATE_CENTROIDS$Location_Name, g$Location_Name)]
    g[, Region := unname(UF_REGION[Location_Name])]
    g
  })
}

# State polygons joined to Mean_Inc (same weights as the geospatial stage)
moran_frame <- function(gdf_states, g) {
  out <- merge(gdf_states, as.data.frame(g)[, c("Location_Name", "Mean_Inc")], by = "Location_Name", all.x = TRUE)
  out <- out[order_c(out$Location_Name), ]
  out[!is.na(out$Mean_Inc), ]
}

run_moran <- function(spatial, gdf_states) {
  res <- list()
  for (virus in names(spatial)) {
    gc <- moran_frame(gdf_states, spatial[[virus]])
    if (nrow(gc) < 5) next
    w <- knn_weights(gc)
    m <- moran_global(gc$Mean_Inc, w, permutations = 999, seed = PERM_SEED)
    res[[virus]] <- list(I = m$I, EI = m$EI, z_sim = m$z_sim, p_sim = m$p_sim, Significant = m$p_sim < ALPHA)
    cat(sprintf("  [%-12s] I = %.4f | p = %.4f | %s\n", virus, m$I, m$p_sim, if (m$p_sim < ALPHA) "Significant" else "NS"))
  }
  res
}

# Moran scatterplot (Figure 4). The mean-incidence map is drawn once, for all three diseases, by
# plot_combined_map() in arbovirus_geospatial.R (Figure 6).
plot_spatial <- function(spatial, moran_results, gdf_states) {
  panels <- lapply(names(spatial), function(virus) {
    gc <- moran_frame(gdf_states, spatial[[virus]])
    vals <- gc$Mean_Inc
    if (length(vals) < 5) return(ggplot() + theme_void())
    w <- knn_weights(gc)
    z <- (vals - mean(vals)) / sqrt(mean((vals - mean(vals))^2))
    wz <- spatial_lag(w, z)
    fit <- stats::lm.fit(cbind(1, z), wz)$coefficients
    d <- data.frame(z = z, wz = wz, lab = substr(gc$Location_Name, 1, 3))
    xl <- seq(min(z), max(z), length.out = 100)
    mr <- moran_results[[virus]]
    if (!is.null(mr)) caption_note("py_10_moran_scatterplot.png", "Global Moran's I (p, 999 permutations)", sprintf("%.3f (p = %.3f)", mr$I, mr$p_sim), virus)
    ggplot(d, aes(z, wz)) +
      geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey50", linetype = "dotted", linewidth = 0.3) +
      geom_line(data = data.frame(z = xl, wz = fit[2] * xl + fit[1]), linewidth = 0.6) +
      geom_point(shape = 21, fill = COLORS[[virus]], colour = "white", size = mpl_s_to_size(45), stroke = 0.4, alpha = 0.72) +
      geom_text(aes(label = lab), size = 4.8 / .pt, alpha = 0.72, hjust = 0, vjust = 0, family = PLOT_FONT) +
      labs(title = virus, x = "z (standardised incidence)", y = "Wz (spatial lag)") +
      theme_mpl(10) + virus_title(virus, 10) + theme(axis.title = element_text(size = 9))
  })
  caption_note("py_10_moran_scatterplot.png", "Former title", "Moran Scatterplot — Spatial Autocorrelation")
  save_png(wrap_plots(panels, nrow = 1), file.path(OUTPUT_DIR, "py_10_moran_scatterplot.png"), 16, 5)
  cat("    Saved -> py_10_moran_scatterplot.png\n")
}

# ==============================================================================
# 6. EXPORT
# ==============================================================================
export_results <- function(annual, mk_results, sp_results, moran_results, rain_lag2_results = NULL) {
  path <- file.path(OUTPUT_DIR, "py_statistical_results.xlsx")
  annual_all <- rbindlist(lapply(names(annual), function(v) copy(annual[[v]])[, Virus := v]))
  sheets <- list(Annual_Incidence = annual_all, Mann_Kendall = rbindlist(mk_results), Spearman_Lag = rbindlist(sp_results))
  if (length(moran_results))
    sheets$Global_Moran <- rbindlist(lapply(names(moran_results), function(v) c(list(Virus = v), moran_results[[v]])))
  if (length(rain_lag2_results)) sheets$Rain_Lag2_Incidence <- rbindlist(rain_lag2_results)
  write_xlsx_sheets(sheets, path)
  cat(sprintf("    Saved -> %s\n", basename(path)))
}

main <- function() {
  sep <- strrep("=", 62)
  cat(sprintf("\n%s\n  ARBOVIRUS PIPELINE — NORMALIZED DB version\n%s\n", sep, sep))
  cat("\n[1/5] Data loading from normalized DB\n")
  datasets <- load_pipeline_datasets()
  cat("\n[2/5] Descriptive — Annual incidence\n")
  annual <- compute_annual(datasets)
  plot_descriptive(annual, datasets)
  cat("\n[3/5] Mann-Kendall trend\n")
  mk_results <- run_mann_kendall(annual)
  plot_mann_kendall(mk_results)
  cat("\n[4/5] Spearman (lag 0-2)\n")
  sp_results <- run_spearman(datasets)
  plot_spearman(sp_results)
  cat("\n[4b/5] Spearman — rho-only / sig rain\n")
  plot_spearman_rho_only(sp_results)
  plot_spearman_rain_significant(sp_results)
  cat("\n[4c/5] Rain Lag-2 (explicit)\n")
  rain_lag2_results <- run_spearman_rain_lag2(datasets)
  plot_spearman_rain_lag2(rain_lag2_results)
  cat("\n[5/5] Spatial — Moran\n")
  spatial <- build_spatial(annual)
  gdf_states <- get_brazil_states()
  moran_results <- run_moran(spatial, gdf_states)
  plot_spatial(spatial, moran_results, gdf_states)
  cat("\n  Exporting...\n")
  export_results(annual, mk_results, sp_results, moran_results, rain_lag2_results)
  write_caption_notes("arbovirus_pipeline_en.R")
  cat(sprintf("\n%s\n  Done! -> %s\n%s\n\n", sep, OUTPUT_DIR, sep))
}

if (sys.nframe() == 0L) main()
