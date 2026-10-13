# ================================================================================
#  GRAPHPAD PRISM (.pzfx) -> R PLOTTING PIPELINE
# ================================================================================
#  Reads the three GraphPad Prism files for Brazilian arbovirus surveillance
#  (Dengue 2016-2025, Chikungunya 2017-2025, Zika 2016-2025) and draws the
#  clinical graphs with GraphPad-Prism-like styling.
#
#  Plot recipes per table type
#  ---------------------------
#    Monthly YYYY-YYYY              -> heat-map (years x months)
#    Sex / ELISA / IgM / Hosp /
#    Serotype / Evolution /
#    Criteria / Classification /
#    Pregnancy                      -> dot-line chart of % per category over the years
#                                      (ignored/blank excluded; valid n under each year)
#    Age group YYYY-YYYY            -> violin plot of the reconstructed age distribution
#    Federative unit YYYY-YYYY      -> horizontal bar chart of cumulative totals
#                                      (colored by Brazilian macro-region)
#    General Data (Dengue only)     -> exported to CSV (data dictionary only)
#  Output: PNG files in <repo>/outputs, all prefixed with the virus.
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
quiet_library("ggplot2", "xml2")

INPUT_DIR <- file.path(DB_DIR, "graphpad_data")
FILES <- c(Dengue = "GPad_Dengue_analysis.pzfx", Chikungunya = "GPad_chikungunya_analysis.pzfx", Zika = "GPad_zika_analysis.pzfx")

# Virus signature colors (consistent with the epidemiological pipeline)
VIRUS_COLOR <- c(Dengue = "#C0392B", Chikungunya = "#E67E22", Zika = "#2980B9")
REGION_COLOR <- c(Norte = "#1B7837", Nordeste = "#762A83", Sudeste = "#E08214", Sul = "#4393C3", "Centro-Oeste" = "#D6604D")

# ── GraphPad-Prism look ───────────────────────────────────────────────────────
theme_prism <- function() {
  theme_classic(base_size = 11, base_family = PLOT_FONT) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          axis.title = element_text(face = "bold", size = 11),
          axis.text = element_text(colour = "black", size = 11),
          axis.line = element_line(linewidth = 1.4 / 2.13, colour = "black"),
          axis.ticks = element_line(linewidth = 1.4 / 2.13, colour = "black"),
          axis.ticks.length = unit(5, "pt"),
          legend.text = element_text(size = 9), legend.key = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}
virus_title <- function(virus) theme(plot.title = element_text(colour = VIRUS_COLOR[[virus]]))
space_thousands <- function(x) formatC(trunc(x), format = "d", big.mark = " ")

# ── Ignored / blank records ───────────────────────────────────────────────────
# SINAN "ignored / blank" columns of the GraphPad tables. They are left out of every plot
# (percentages are of records with a known value) and each exclusion is logged by log_exclusion().
MISSING_CATS <- c("Ign/Blank", "Ign", "Blank", "FU - Ign/exterior")

# No in-plot titles: the disease name is shown only on Figure 5 panels (the one clinical figure that
# spans the three diseases); Figures 7/8 and the supplementary plots are named by their captions.
panel_title <- function(fname, virus) if (startsWith(figure_of(fname), "Figure 5")) virus else NULL
note_title <- function(fname, virus, table_title, suffix = "") caption_note(fname, "Former title", sprintf("%s — %s%s", virus, table_title, suffix))

EXCLUSIONS <- list()   # one data.table per plotted table, filled by log_exclusion()

# Records excluded from the plot of one table: one row per year plus an "All" row.
log_exclusion <- function(df, virus, table_title, fname) {
  cats <- setdiff(names(df), "Year")
  miss <- intersect(MISSING_CATS, cats)
  tot <- rowSums(as.matrix(df[cats]), na.rm = TRUE)
  exc <- if (length(miss)) rowSums(as.matrix(df[miss]), na.rm = TRUE) else 0 * tot
  d <- data.table(Year = c(as.character(df$Year), "All"), n_total = c(tot, sum(tot)), n_excluded = c(exc, sum(exc)))
  d[, `:=`(n_valid = n_total - n_excluded, pct_excluded = py_round(100 * n_excluded / n_total, 2))]
  for (cc in c("n_total", "n_excluded", "n_valid")) set(d, j = cc, value = as.integer(d[[cc]]))
  EXCLUSIONS[[length(EXCLUSIONS) + 1]] <<- data.table(
    Figure = figure_of(fname), Output_file = fname,
    Virus = virus, Table = table_title, d, Excluded_categories = paste(miss, collapse = "; "))
}

# One caption sentence per composite figure, from the "All" rows of the log.
exclusion_captions <- function(log) {
  fmt <- function(x) formatC(x, format = "d", big.mark = ",")
  all <- log[Year == "All" & Figure != "Supplementary"]
  all[, Variable := sub("^Classfication", "Classification", sub(" [0-9].*$", "", Table))]   # Prism file typo
  all[!Variable %in% c("ELISA", "IgM"), Variable := tolower(Variable)]
  all[, pct := ifelse(pct_excluded < 0.1, "<0.1%", sprintf("%.1f%%", pct_excluded))]
  all[, Fig := sub("[A-Z]$", "", Figure)]
  vapply(sort_c(unique(all$Fig)), function(f) {
    s <- all[Fig == f][order(match(sub("\\.png$", "", Output_file), names(FIGURE_OF)))]   # composite panel order
    parts <- sprintf("%s, %s %s: n = %s of %s (%s)", sub("^Figure [0-9]+", "", s$Figure), s$Virus, s$Variable,
                     fmt(s$n_excluded), fmt(s$n_total), s$pct)
    what <- if (any(grepl("violin", s$Output_file))) "Percentages and age distributions are" else "Percentages are"
    sprintf("%s. %s calculated on records with a known value; records with ignored/blank values were excluded: %s.",
            f, what, paste(parts, collapse = "; "))
  }, "")
}

# seaborn desaturate(color, prop): scale the HLS saturation
desaturate <- function(col, prop) {
  rgb <- grDevices::col2rgb(col)[, 1] / 255
  mx <- max(rgb); mn <- min(rgb); l <- (mx + mn) / 2
  if (mx == mn) return(col)
  d <- mx - mn
  s <- if (l > 0.5) d / (2 - mx - mn) else d / (mx + mn)
  h <- if (mx == rgb[1]) ((rgb[2] - rgb[3]) / d) %% 6 else if (mx == rgb[2]) (rgb[3] - rgb[1]) / d + 2 else (rgb[1] - rgb[2]) / d + 4
  h <- h / 6; s <- s * prop
  q <- if (l < 0.5) l * (1 + s) else l + s - l * s
  p <- 2 * l - q
  hue <- function(t) { t <- t %% 1; if (t < 1/6) p + (q - p) * 6 * t else if (t < 1/2) q else if (t < 2/3) p + (q - p) * (2/3 - t) * 6 else p }
  grDevices::rgb(hue(h + 1/3), hue(h), hue(h - 1/3))
}

# ==============================================================================
# 1. PARSE THE .pzfx FILES
# ==============================================================================
parse_pzfx <- function(path) {
  doc <- xml_ns_strip(read_xml(path))
  txt <- function(node) if (length(node) && !inherits(node, "xml_missing")) trimws(xml_text(node)) else ""
  lapply(xml_find_all(doc, "/*/Table"), function(t) {
    title <- txt(xml_find_first(t, "./Title"))
    sub <- xml_find_first(t, "./RowTitlesColumn/Subcolumn")
    row_titles <- if (!inherits(sub, "xml_missing")) trimws(xml_text(xml_find_all(sub, "./d"))) else character()
    columns <- lapply(xml_find_all(t, "./YColumn"), function(yc) {
      sc <- xml_find_first(yc, "./Subcolumn")
      vals <- if (!inherits(sc, "xml_missing")) suppressWarnings(as.numeric(trimws(xml_text(xml_find_all(sc, "./d"))))) else numeric()
      list(title = txt(xml_find_first(yc, "./Title")), values = vals)
    })
    list(id = xml_attr(t, "ID"), title = title, row_titles = row_titles, columns = columns)
  })
}

# rows = row titles (named "Year"), one column per series. A repeated series title keeps the
# position of its first occurrence and the values of its last (Python dict semantics).
table_to_df <- function(tbl) {
  out <- list()
  for (cc in tbl$columns) out[[cc$title]] <- cc$values
  df <- as.data.frame(out, check.names = FALSE, optional = TRUE)
  if (!length(out)) df <- data.frame(row.names = seq_along(tbl$row_titles))
  data.frame(Year = tbl$row_titles, df, check.names = FALSE)
}

# ==============================================================================
# 2. PLOT RECIPES
# ==============================================================================
save_plot <- function(p, fname, w, h) {
  save_png(p, file.path(OUTPUT_DIR, fname), w, h, dpi = 150)
  cat(sprintf("    Saved -> %s\n", fname))
}

plot_monthly <- function(df, virus, table_title) {
  month_order <- c("JAN", "FEV", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
  months <- month_order[month_order %in% names(df)]
  years <- df$Year
  long <- data.frame(Year = factor(rep(years, length(months)), rev(years)),
                     Month = factor(rep(months, each = length(years)), months),
                     v = unlist(df[months], use.names = FALSE))
  vmax <- max(long$v, na.rm = TRUE)
  lab <- long[!is.na(long$v) & long$v > 0, ]
  p <- ggplot(long, aes(Month, Year, fill = v)) + geom_raster() +
    scale_fill_gradientn(colours = c("#FFFFFF", VIRUS_COLOR[[virus]]), name = "Cases", labels = space_thousands,
                         guide = guide_colourbar(barheight = unit(0.8 * (0.45 * length(years) + 2) - 0.6, "in"),
                                                 title.position = "right", title.theme = element_text(angle = 90, size = 11))) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Month", y = "Year") +
    theme_prism() + theme(axis.line = element_blank())
  if (length(months) * length(years) <= 13 * 12)
    p <- p + geom_text(data = lab, aes(label = space_thousands(v), colour = v > 0.55 * vmax), size = 7 / .pt, family = PLOT_FONT) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none")
  fname <- sprintf("%s_01b_monthly_heatmap.png", virus)
  note_title(fname, virus, table_title, " (heat-map)")
  save_plot(p, fname, 11, 0.45 * length(years) + 2)
}

# ── Proportion charts (dot-line) ──────────────────────────────────────────────
# Colour-blind-safe categorical order (checked with the dataviz palette validator). A category keeps
# its slot by meaning: its pzfx column order, or CATEGORY_SLOT where one concept has different names.
PROP_PALETTE <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
PROP_SHAPES  <- c(16, 17, 15, 18, 25, 4)   # secondary encoding, so the lines stay distinct in grayscale
CATEGORY_SLOT <- c(Confirmed = 1, Discarded = 2)
# Plot labels for the GraphPad category names (typos / terse codes); source data keeps the raw names.
CATEGORY_LABEL <- c(M = "Male", F = "Female", "Not accomplished" = "Not performed", Descarded = "Discarded",
                    Chikungunya = "Confirmed", "Clinic-epidemiological" = "Clinical-epidemiological",
                    "1º Trimester" = "1st trimester", "2º Trimester" = "2nd trimester", "3º Trimester" = "3rd trimester",
                    "Pregnancy age ignored" = "Pregnant, trimester unknown", "Death by other cause" = "Death (other cause)",
                    "Other causes" = "Death (other cause)", "Death in investigation" = "Death under investigation",
                    "Non pregnant" = "Not pregnant")
# Labels that only apply to one table (keyed like VARIABLE_LABEL); they take precedence over CATEGORY_LABEL.
TABLE_CATEGORY_LABEL <- list(Hospitalization = c(Positive = "Hospitalized", Negative = "Not hospitalized"))
# Legend title naming the variable of each table (first word of the pzfx table title); it identifies
# the panel now that there is no plot title (ELISA and IgM share the same categories).
VARIABLE_LABEL <- c(ELISA = "ELISA result", IgM = "IgM result", Hospitalization = "Hospitalization", Sex = "Sex",
                    Serotype = "Serotype", Evolution = "Outcome", Criteria = "Confirmation criterion",
                    Classification = "Final classification", Classfication = "Final classification",   # Prism file typo
                    Pregnancy = "Pregnancy status")
INK_2 <- "#333333"; GRID_LINE <- "#e1e0d9"

# 547 -> "547", 4078 -> "4.1k", 189618 -> "190k", 4275943 -> "4.3M"
n_compact <- function(n) ifelse(is.na(n), "", ifelse(n >= 1e6, sprintf("%.1fM", n / 1e6),
                                ifelse(n >= 1e4, sprintf("%.0fk", n / 1e3), ifelse(n >= 1e3, sprintf("%.1fk", n / 1e3), sprintf("%d", as.integer(n))))))

# Spread end-of-line label positions (0-100 scale) so neighbours are at least `gap` apart.
spread_labels <- function(y, gap) {
  o <- order(y); v <- y[o]
  for (i in seq_along(v)[-1]) v[i] <- max(v[i], v[i - 1] + gap)
  if (v[length(v)] > 100) { v[length(v)] <- 100; for (i in rev(seq_along(v))[-1]) v[i] <- min(v[i], v[i + 1] - gap) }
  y[o] <- v
  y
}

# Share of records with a known value per category and year: one line per category, direct
# labels at the last year, valid n under each year (ignored/blank columns excluded).
plot_proportion <- function(df, virus, table_title, slug) {
  fname <- sprintf("%s_%s.png", virus, slug)
  log_exclusion(df, virus, table_title, fname)
  cats <- setdiff(names(df), c("Year", MISSING_CATS))
  m <- as.matrix(df[cats])
  totals <- rowSums(m, na.rm = TRUE); totals[totals == 0] <- NA
  pct <- m / totals * 100
  years <- as.integer(df$Year)
  key <- strsplit(table_title, " ")[[1]][1]
  labels <- c(TABLE_CATEGORY_LABEL[[key]], CATEGORY_LABEL[setdiff(names(CATEGORY_LABEL), names(TABLE_CATEGORY_LABEL[[key]]))])
  lab <- ifelse(cats %in% names(labels), unname(labels[cats]), cats)
  slot <- seq_along(cats); hit <- lab %in% names(CATEGORY_SLOT); slot[hit] <- CATEGORY_SLOT[lab[hit]]
  lv <- lab[order(slot)]
  long <- data.frame(Year = rep(years, length(cats)), Cat = factor(rep(lab, each = nrow(df)), lv), pct = as.vector(pct))
  long <- long[!is.na(long$pct), ]
  end <- do.call(rbind, lapply(split(long, long$Cat), function(d) d[which.max(d$Year), ]))
  end$lab <- sprintf("%s %s", end$Cat, ifelse(end$pct < 1, "<1%", sprintf("%.0f%%", end$pct)))
  end$ylab <- spread_labels(end$pct, 4.5)
  cols <- setNames(PROP_PALETTE[sort(slot)], lv); shapes <- setNames(PROP_SHAPES[sort(slot)], lv)
  var_lab <- if (key %in% names(VARIABLE_LABEL)) VARIABLE_LABEL[[key]] else key
  note_title(fname, virus, table_title)
  p <- ggplot(long, aes(Year, pct, colour = Cat, fill = Cat, shape = Cat, group = Cat)) +
    geom_hline(yintercept = seq(0, 100, 25), colour = GRID_LINE, linewidth = 0.3) +
    geom_line(linewidth = 0.9, lineend = "round", linejoin = "round") +
    geom_point(colour = "white", fill = "white", shape = 21, size = 3.6, stroke = 0, show.legend = FALSE) +   # surface ring
    geom_point(size = 2.3, stroke = 0.9) +
    geom_text(data = end, aes(x = Year + 0.22, y = ylab, label = lab), hjust = 0, colour = INK_2, size = 8 / .pt,
              family = PLOT_FONT, inherit.aes = FALSE) +
    scale_colour_manual(values = cols, name = var_lab) + scale_fill_manual(values = cols, name = var_lab) +
    scale_shape_manual(values = shapes, name = var_lab) +
    guides(colour = guide_legend(nrow = ceiling(length(lv) / 4), byrow = TRUE)) +
    scale_x_continuous(breaks = years, labels = sprintf("%d\nn=%s", years, n_compact(totals)), expand = expansion(add = 0.35)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(add = c(1.5, 3))) +
    coord_cartesian(clip = "off") +
    labs(title = panel_title(fname, virus), x = NULL, y = "Records with a known value (%)") +
    theme_prism() + virus_title(virus) +
    theme(axis.text.x = element_text(size = 8.5, lineheight = 1.05), legend.position = "top", legend.justification = "left",
          legend.text = element_text(size = 9), legend.title = element_text(size = 9, face = "bold"),
          legend.key.width = unit(22, "pt"),
          plot.margin = margin(5.5, max(nchar(end$lab)) * 4.6 + 10, 5.5, 5.5))
  save_plot(p, fname, 9, 5.5)
}

# Numeric ranges for each age bin (used to reconstruct a continuous age axis
# from binned counts).  "<1 year" -> [0,1); "80 e +" capped at 100.
AGE_BIN_RANGES <- list("<1 year" = c(0, 1), "1-4" = c(1, 5), "5-9" = c(5, 10), "10-14" = c(10, 15), "15-19" = c(15, 20),
                       "20-39" = c(20, 40), "40-59" = c(40, 60), "60-64" = c(60, 65), "65-69" = c(65, 70),
                       "70-79" = c(70, 80), "80 e +" = c(80, 100))
N_AGE_SAMPLES <- 6000   # samples per year used to build the violin

# Synthetic ages per year: each sample picks a bin with probability proportional to its count,
# then a uniform age inside the bin. 'Ign/Blank' is dropped (age unknown). Seeded (R's RNG).
reconstruct_ages <- function(df, seed = 42) {
  set.seed(seed)
  available <- names(AGE_BIN_RANGES)[names(AGE_BIN_RANGES) %in% names(df)]
  if (!length(available)) return(NULL)
  lows <- vapply(AGE_BIN_RANGES[available], `[`, 0, 1); highs <- vapply(AGE_BIN_RANGES[available], `[`, 0, 2)
  recs <- list()
  for (i in seq_len(nrow(df))) {
    w <- as.numeric(unlist(df[i, available]))
    w <- ifelse(is.finite(w) & w > 0, w, 0)
    if (sum(w) <= 0) next
    idx <- sample.int(length(available), N_AGE_SAMPLES, replace = TRUE, prob = w / sum(w))
    ages <- stats::runif(N_AGE_SAMPLES, lows[idx], highs[idx])
    recs[[length(recs) + 1]] <- data.frame(Year = as.character(df$Year[i]), Age = ages,
                                           n_cases = sum(as.numeric(unlist(df[i, available])), na.rm = TRUE))
  }
  if (!length(recs)) return(NULL)
  list(ages = do.call(rbind, recs), available = available)
}

plot_age_group <- function(df, virus, table_title) {
  rec <- reconstruct_ages(df)
  if (is.null(rec)) { cat(sprintf("    (no age bins found, skipping %s)\n", table_title)); return(invisible()) }
  fname <- sprintf("%s_04_age_group_violin.png", virus)
  log_exclusion(df, virus, table_title, fname)
  d <- rec$ages
  year_order <- unique(d$Year)
  d$Year <- factor(d$Year, year_order)
  meds <- aggregate(Age ~ Year, d, stats::median)
  ns <- unique(d[c("Year", "n_cases")])
  col <- desaturate(VIRUS_COLOR[[virus]], 0.85)
  p <- ggplot(d, aes(Year, Age))
  # seaborn: Scott bandwidth (n^-1/5 * sd) x bw_adjust 0.6, per violin; cut=0; equal areas
  for (yr in year_order) {
    s <- d[d$Year == yr, ]
    bw <- 0.6 * stats::sd(s$Age) * nrow(s)^(-1 / 5)
    p <- p + geom_violin(data = s, fill = col, colour = "#3f3f3f", linewidth = 1.0 / 2.13, bw = bw, trim = TRUE,
                         width = 0.8, quantile.linetype = "dashed", quantiles = c(0.25, 0.5, 0.75))
  }
  p <- p + geom_point(data = meds, shape = 21, fill = "white", colour = "black", size = mpl_s_to_size(22), stroke = 0.9 / 2.13 * 2) +
    geom_text(data = ns, aes(Year, 102, label = paste("n =", space_thousands(n_cases))), vjust = 0, size = 7.5 / .pt,
              colour = "#555555", family = PLOT_FONT) +
    scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 10), expand = c(0, 0)) +
    labs(title = panel_title(fname, virus), x = "Year", y = "Patient age (years)") +
    theme_prism() + virus_title(virus) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_line(linetype = "dotted", colour = scales::alpha("grey40", 0.45), linewidth = 0.35))
  note_title(fname, virus, table_title)
  caption_note(fname, "Key (former in-plot note)", paste("Ages reconstructed from binned counts (jittered uniformly inside each bin);",
                                                         "inner lines = quartiles; white dot = median; n = cases with known age."))
  save_plot(p, fname, max(9, 0.85 * length(year_order) + 3), 6)
}

plot_federative_unit <- function(df, virus, table_title) {
  fname <- sprintf("%s_05_federative_unit_totals.png", virus)
  log_exclusion(df, virus, table_title, fname)
  note_title(fname, virus, table_title)
  cols <- setdiff(names(df), "Year")
  cols <- cols[!startsWith(cols, "FU - ")]
  totals <- sort(colSums(as.matrix(df[cols]), na.rm = TRUE))
  d <- data.frame(State = factor(names(totals), names(totals)), v = unname(totals))
  d$Region <- ifelse(as.character(d$State) %in% names(UF_REGION), UF_REGION[as.character(d$State)], "Other")
  xmax <- max(d$v)
  p <- ggplot(d, aes(v, State, fill = Region)) + geom_col(colour = "white", linewidth = 0.4 / 2.13, width = 0.8) +
    geom_text(aes(x = v + 0.01 * xmax, label = space_thousands(v)), hjust = 0, size = 8 / .pt, family = PLOT_FONT) +
    scale_fill_manual(values = c(REGION_COLOR, Other = "#888888"), breaks = names(REGION_COLOR), limits = c(names(REGION_COLOR), "Other"),
                      name = "Region") +
    scale_x_continuous(limits = c(0, xmax * 1.15), labels = space_thousands, expand = c(0, 0)) +
    labs(x = "Total reported cases (period)", y = NULL) +
    theme_prism() +
    theme(panel.grid.major.x = element_line(linetype = "dotted", colour = scales::alpha("grey40", 0.45), linewidth = 0.35),
          legend.position = "inside", legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0),
          legend.text = element_text(size = 8))
  save_plot(p, fname, 9.5, 8.5)
}

# ==============================================================================
# 3. DISPATCH PER TABLE
# ==============================================================================
PLOT_SLUGS <- c(ELISA = "02_ELISA", IgM = "02_IgM", Hospitalization = "02_hospitalization", Sex = "03_sex",
                Serotype = "03_serotype", Evolution = "03_evolution", Criteria = "03_criteria",
                Classification = "03_classification", Classfication = "03_classification",  # Prism file typo
                Pregnancy = "03_pregnancy")

dispatch_table <- function(tbl, virus) {
  title <- tbl$title
  cat(sprintf("  • %s  (%d rows × %d cols)\n", title, length(tbl$row_titles), length(tbl$columns)))
  df <- table_to_df(tbl)
  if (title == "General Data") {
    # Export the raw master table to CSV; no single graph would be readable
    out <- file.path(OUTPUT_DIR, sprintf("%s_00_general_data.csv", virus))
    write_csv_py(df, out)
    cat(sprintf("    Exported -> %s  (74-column master sheet)\n", basename(out)))
    return(invisible())
  }
  if (startsWith(title, "Monthly")) plot_monthly(df, virus, title)
  else if (startsWith(title, "Age group")) plot_age_group(df, virus, title)
  else if (startsWith(title, "Federative unit")) plot_federative_unit(df, virus, title)
  else {
    # Generic categorical table (ELISA, IgM, Sex, Serotype, Evolution, ...)
    key <- strsplit(title, " ")[[1]][1]
    slug <- if (key %in% names(PLOT_SLUGS)) PLOT_SLUGS[[key]] else paste0("03_", tolower(key))
    plot_proportion(df, virus, title, slug)
  }
}

read_tables <- function(virus) parse_pzfx(file.path(INPUT_DIR, FILES[[virus]]))
find_table <- function(tables, prefix) {
  for (t in tables) if (startsWith(t$title, prefix)) return(table_to_df(t))
  NULL
}

main <- function() {
  sep <- strrep("=", 70)
  cat(sprintf("\n%s\n  GRAPHPAD PRISM (.pzfx) -> R PLOTTING PIPELINE\n%s\n", sep, sep))
  for (virus in names(FILES)) {
    cat(sprintf("\n[%s] reading %s\n", virus, FILES[[virus]]))
    tables <- read_tables(virus)
    cat(sprintf("  Found %d tables\n", length(tables)))
    for (tbl in tables) dispatch_table(tbl, virus)
  }
  log <- rbindlist(EXCLUSIONS)
  write_csv_py(log, file.path(OUTPUT_DIR, "figure_exclusions.csv"))
  writeLines(exclusion_captions(log), file.path(OUTPUT_DIR, "figure_exclusions_captions.txt"), useBytes = TRUE)
  cat(sprintf("\n  Excluded ignored/blank records logged -> figure_exclusions.csv (%d tables), figure_exclusions_captions.txt\n",
              length(EXCLUSIONS)))
  write_caption_notes("graphpad_pipeline.R")
  cat(sprintf("\n%s\n  Done. All outputs written to -> %s\n%s\n\n", sep, OUTPUT_DIR, sep))
}

if (sys.nframe() == 0L) main()
