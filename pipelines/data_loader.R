# ================================================================================
#  DATA LOADER — Normalized Schema Adapter
# ================================================================================
#  Loads the normalized DB (CSV-based 3NF schema) and produces per-virus
#  long-format tables in the shape every pipeline expects.
#
#  Normalized schema:
#    dim_location      : Location_Name, State_Code, Location_ID
#    fact_cases        : Location_ID, Disease, Year, Month, Cases, Monthly_Inc_per100k
#    fact_climate      : Location_ID, Year, Month, Temp_Avg, Rain_mm,
#                        Temp_Lag1, Temp_Lag2, Rain_Lag1, Rain_Lag2
#    fact_population   : Location_ID, Year, Population
#
#  Disease key mapping:  DENV -> Dengue, CHIKV -> Chikungunya, ZIKV -> Zika
#  Filter: Year >= 2016 for all diseases.
#  Sourced after R/utils.R.
# ================================================================================
quiet_library("data.table")
source(path_root("R", "stats.R"))

DB_DIR      <- path_root("db")
OUTPUT_DIR  <- path_root("outputs")
FIGURES_DIR <- path_root("figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)

DISEASE_KEY <- c(DENV = "Dengue", CHIKV = "Chikungunya", ZIKV = "Zika")

UF_REGION <- c(
  "Rondonia" = "Norte", "Acre" = "Norte", "Amazonas" = "Norte", "Roraima" = "Norte", "Para" = "Norte", "Amapa" = "Norte",
  "Tocantins" = "Norte", "Maranhao" = "Nordeste", "Piaui" = "Nordeste", "Ceara" = "Nordeste",
  "Rio Grande do Norte" = "Nordeste", "Paraiba" = "Nordeste", "Pernambuco" = "Nordeste", "Alagoas" = "Nordeste",
  "Sergipe" = "Nordeste", "Bahia" = "Nordeste", "Minas Gerais" = "Sudeste", "Espirito Santo" = "Sudeste",
  "Rio de Janeiro" = "Sudeste", "Sao Paulo" = "Sudeste", "Parana" = "Sul", "Santa Catarina" = "Sul",
  "Rio Grande do Sul" = "Sul", "Mato Grosso do Sul" = "Centro-Oeste", "Mato Grosso" = "Centro-Oeste",
  "Goias" = "Centro-Oeste", "Distrito Federal" = "Centro-Oeste")

YEAR_START <- 2016

load_normalized_tables <- function() {
  list(loc     = fread(file.path(DB_DIR, "normalized_dim_location.csv"), encoding = "UTF-8"),
       cases   = fread(file.path(DB_DIR, "normalized_fact_cases.csv")),
       climate = fread(file.path(DB_DIR, "normalized_fact_climate.csv")),
       pop     = fread(file.path(DB_DIR, "normalized_fact_population.csv")))
}

# list(Dengue=, Chikungunya=, Zika=) of state-month tables, sorted by Location_Name, Year, Month.
load_pipeline_datasets <- function(year_start = YEAR_START, verbose = TRUE) {
  tbl <- load_normalized_tables()
  cases <- tbl$cases[Year >= year_start]
  clim  <- tbl$climate[Year >= year_start]
  pop   <- tbl$pop[Year >= year_start]
  out <- list()
  for (code in names(DISEASE_KEY)) {
    virus <- DISEASE_KEY[[code]]
    d <- cases[Disease == code]
    d[, Cases := as.double(Cases)]
    d <- merge(d, clim, by = c("Location_ID", "Year", "Month"), all.x = TRUE, sort = FALSE)
    d <- merge(d, pop, by = c("Location_ID", "Year"), all.x = TRUE, sort = FALSE)
    d <- merge(d, tbl$loc, by = "Location_ID", all.x = TRUE, sort = FALSE)
    setnames(d, "State_Code", "UF_Code")
    d[, Date := as.Date(sprintf("%d-%02d-01", Year, Month))]
    d[, Region := unname(UF_REGION[Location_Name])]
    d[, Monthly_Inc_100k := Monthly_Inc_per100k]          # legacy alias
    d <- d[order_c(Location_Name, Year, Month)]
    setcolorder(d, c("Location_ID", "Disease", "Year", "Month", "Cases", "Monthly_Inc_per100k", "Temp_Avg", "Rain_mm",
                     "Temp_Lag1", "Temp_Lag2", "Rain_Lag1", "Rain_Lag2", "Population", "Location_Name", "UF_Code",
                     "Date", "Region", "Monthly_Inc_100k"))
    out[[virus]] <- d
    if (verbose) cat(sprintf("  [%-12s] %5d obs | %d states | %d-%d\n", virus, nrow(d), uniqueN(d$Location_Name),
                             min(d$Year), max(d$Year)))
  }
  out
}

# Annual totals and incidence per (state, year). Months with missing case counts are not treated as
# zero: totals are annualised as Total_Cases * 12 / Months_Reported.
compute_annual <- function(datasets) {
  lapply(datasets, function(df) {
    ann <- df[, .(Total_Cases = sum_min1(Cases), Months_Reported = sum(!is.na(Cases)),
                  Population = Population[!is.na(Population)][1]),
              by = .(Location_Name, UF_Code, Region, Year)]
    ann <- ann[order_c(Location_Name, UF_Code, Region, Year)]
    mr <- ifelse(ann$Months_Reported > 0, ann$Months_Reported, NA)
    annualised <- ann$Total_Cases * 12 / mr
    ann[, Annual_Inc_100k := (annualised / Population) * 1e5]
    ann[, Months_Reported := as.integer(Months_Reported)]
    ann
  })
}

# Population-weighted mean of `col`, ignoring rows where it is NA.
pop_weighted_mean <- function(g, col = "Annual_Inc_100k") {
  v <- g[[col]]; ok <- !is.na(v)
  if (!any(ok)) return(NaN)
  w <- g$Population[ok]
  sum(v[ok] * w) / sum(w)
}

# ---- State boundaries (shared by the spatial stages) --------------------------------------------
GEO_CACHE <- file.path(DB_DIR, "brazil_states.geojson")

IBGE_CODE_TO_NAME <- c("11" = "Rondonia", "12" = "Acre", "13" = "Amazonas", "14" = "Roraima", "15" = "Para", "16" = "Amapa",
  "17" = "Tocantins", "21" = "Maranhao", "22" = "Piaui", "23" = "Ceara", "24" = "Rio Grande do Norte", "25" = "Paraiba",
  "26" = "Pernambuco", "27" = "Alagoas", "28" = "Sergipe", "29" = "Bahia", "31" = "Minas Gerais", "32" = "Espirito Santo",
  "33" = "Rio de Janeiro", "35" = "Sao Paulo", "41" = "Parana", "42" = "Santa Catarina", "43" = "Rio Grande do Sul",
  "50" = "Mato Grosso do Sul", "51" = "Mato Grosso", "52" = "Goias", "53" = "Distrito Federal")

NAME_TO_ABBR <- c("Rondonia" = "RO", "Acre" = "AC", "Amazonas" = "AM", "Roraima" = "RR", "Para" = "PA", "Amapa" = "AP",
  "Tocantins" = "TO", "Maranhao" = "MA", "Piaui" = "PI", "Ceara" = "CE", "Rio Grande do Norte" = "RN", "Paraiba" = "PB",
  "Pernambuco" = "PE", "Alagoas" = "AL", "Sergipe" = "SE", "Bahia" = "BA", "Minas Gerais" = "MG", "Espirito Santo" = "ES",
  "Rio de Janeiro" = "RJ", "Sao Paulo" = "SP", "Parana" = "PR", "Santa Catarina" = "SC", "Rio Grande do Sul" = "RS",
  "Mato Grosso do Sul" = "MS", "Mato Grosso" = "MT", "Goias" = "GO", "Distrito Federal" = "DF")

.fetch_state <- function(code) {
  url <- sprintf("https://servicodados.ibge.gov.br/api/v3/malhas/estados/%s?formato=application/vnd.geo+json&qualidade=intermediaria", code)
  tmp <- tempfile(fileext = ".geojson")
  curl::curl_download(url, tmp, quiet = TRUE)
  g <- sf::st_read(tmp, quiet = TRUE)
  g$codarea <- code
  g[, "codarea"]
}

.attach_metadata <- function(gdf) {
  gdf$Location_Name <- unname(IBGE_CODE_TO_NAME[as.character(gdf$codarea)])
  gdf <- gdf[!is.na(gdf$Location_Name), ]
  nm <- sort_c(unique(gdf$Location_Name))
  geoms <- lapply(nm, function(n) {
    g <- sf::st_geometry(gdf[gdf$Location_Name == n, ])
    if (length(g) == 1) g[[1]] else sf::st_union(g)[[1]]
  })
  out <- sf::st_sf(Location_Name = nm, geometry = sf::st_sfc(geoms, crs = sf::st_crs(gdf)))
  out$Region <- unname(UF_REGION[out$Location_Name])
  out$Abbr <- unname(NAME_TO_ABBR[out$Location_Name])
  cen <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(out))))
  out$centroid_x <- cen[, 1]; out$centroid_y <- cen[, 2]
  out
}

# Spatial weights used by every Moran / LISA computation: KNN on polygon centroids
KNN_K <- 5
knn_weights <- function(gdf) knn_listw(cbind(gdf$centroid_x, gdf$centroid_y), k = KNN_K)

# State polygons from db/brazil_states.geojson (downloaded from IBGE if absent).
get_brazil_states <- function() {
  quiet_library("sf")
  suppressMessages(sf::sf_use_s2(FALSE))   # planar centroids in lon/lat, as geopandas/shapely compute them
  if (file.exists(GEO_CACHE)) {
    gdf <- sf::st_read(GEO_CACHE, quiet = TRUE)
    if (nrow(gdf) >= 26) {
      cat(sprintf("  Loaded %d states from cache (%s)\n", nrow(gdf), basename(GEO_CACHE)))
      return(.attach_metadata(gdf))
    }
    cat(sprintf("  Cache invalid (%d rows), re-downloading...\n", nrow(gdf)))
  }
  codes <- names(IBGE_CODE_TO_NAME)
  cat(sprintf("  Downloading %d state boundaries from IBGE...\n", length(codes)))
  frames <- list()
  for (cd in codes) {
    g <- tryCatch(.fetch_state(cd), error = function(e) { cat(sprintf("    FAILED %s: %s\n", cd, conditionMessage(e))); NULL })
    if (!is.null(g)) frames[[length(frames) + 1]] <- g
  }
  if (length(frames) < 26) stop(sprintf("Only %d/27 states downloaded", length(frames)))
  gdf <- do.call(rbind, frames)
  sf::st_write(gdf, GEO_CACHE, driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
  cat(sprintf("  Cached -> %s\n", basename(GEO_CACHE)))
  .attach_metadata(gdf)
}
