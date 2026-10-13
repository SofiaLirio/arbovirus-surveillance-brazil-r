# Project guide

Where things are, and where to look first.

## Run it

```powershell
powershell -ExecutionPolicy Bypass -File run_all.ps1     # Windows
```
```bash
bash run_all.sh                                          # macOS / Linux / Git Bash
```
```r
system("Rscript run_all.R")                             # from R or RStudio
```

This rebuilds everything from `raw/`, which takes about 7 minutes. The folder already contains the
results, so you only need to run it after changing data or code.

- `-SkipEtl` / `--skip-etl` reuses `db/` and reruns only the analyses (about 4 minutes). **Use this on
  a fresh clone:** `raw/` is not in the git repository (see "Raw data" in `README.md` to download it).
- `-SkipInstall` / `--skip-install` skips installing packages.

Needs R 4.5 or newer (https://cloud.r-project.org) and internet the first time, to install the
CRAN packages listed in `R/install_packages.R`. Rtools and Python are not needed.

## How data flows

```
raw/  ──etl/──▶  db/ + clinical_tables/  ──pipelines/──▶  outputs/  figures/  alt_results/  source_data/
```

## Folders

| Folder | What it is | Edit it? |
|---|---|---|
| `raw/` | Original downloads: SINAN cases, INMET weather, IBGE population, state shapes (not in git; `tools/fetch_raw_sources.R`) | No |
| `reference/` | Small fixed inputs: state list, map shapes, GraphPad clinical files, expected hashes | Rarely |
| `etl/` | Turns `raw/` into `db/` and `clinical_tables/` | Only to change data rules |
| `db/` | Clean tables every analysis reads | No, it is regenerated |
| `R/` | Shared R helpers (`utils.R`, `stats.R`) and the package installer | Rarely |
| `pipelines/` | All analysis code | Yes, to change analyses |
| `outputs/` | Result tables and individual plots | No, it is regenerated |
| `figures/` | Final composite figures (Figure 1–8, S1–S2, GLM figures) | No, it is regenerated |
| `alt_results/` | 3 alternative-view figures | No, it is regenerated |
| `source_data/` | The numbers behind each figure, one Excel file per figure | No, it is regenerated |
| `clinical_tables/` | Sex, age, outcome, serotype and similar breakdowns, 2016–2026 | No, it is regenerated |
| `tools/` | Re-download raw data; compare with old archives | Optional |
| `docs/` | Data sources, the R port notes and older write-ups | Reading |

## The most important files

**Results to read first**
- `figures/Figure1_incidence_regional.png` … `Figure8_chikv_zika_clinical.png`: main figures.
- `outputs/py_statistical_results.xlsx`: annual incidence, trends (Mann-Kendall), climate
  correlations (Spearman), Moran's I.
- `outputs/geo_spatial_stats.xlsx`: spatial clusters (LISA) and maps data.
- `outputs/glm_*_results.xlsx`: climate regression models (Poisson, quasi-Poisson, negative binomial).
- `outputs/2026-10-08/additional/`: advanced analyses. These cover nonlinear climate (DLNM), virus
  co-circulation, spatial regression, COVID effect (ITS), forecasting, and epidemic thresholds.
- `source_data/Figure*.xlsx`: exact numbers plotted in each figure.

**Input data**
- `db/normalized_fact_cases.csv`: cases per state × disease × month.
- `db/normalized_fact_climate.csv`: temperature and rain per state × month, with lags.
- `db/normalized_fact_population.csv`: population per state × year.

**Code entry points**
- `run_all.R` (or the `run_all.ps1` / `run_all.sh` wrappers): runs everything.
- `etl/build_db.R`: builds `db/` from `raw/`.
- `pipelines/run_all.R`: the 9 core analysis steps.
- `pipelines/additional/run_additional_analyses.R`: the advanced analyses.
- `pipelines/data_loader.R`: how every script loads data; `R/utils.R` holds shared helpers.

## Common tasks

| I want to… | Do this |
|---|---|
| Rerun everything | `run_all.ps1` / `run_all.sh` |
| Get the newest SINAN data | `Rscript tools/fetch_raw_sources.R --only sinan`, then `run_all` |
| Change how cases are counted | Edit `etl/build_cases_clinical.R` |
| Change an analysis or plot | Edit the script in `pipelines/`, then `run_all -SkipEtl` |
| Find the numbers behind a figure | Open `source_data/<Figure>.xlsx` |
| Check where raw data came from | `docs/DATA_SOURCES.md` and `raw/MANIFEST.csv` |

## Good to know

- **2026 is partial.** It covers January–September (SINAN release of 2026-10-06). Annual
  values for 2026 are scaled up to 12 months.
- **Population.** 2016–2023 come from IBGE projections; 2024–2026 from IBGE official estimates.
- **Clinical figures 5, 7 and 8 end in 2025.** They come from the GraphPad files.
  `clinical_tables/` covers 2016–2026.
- **`docs/RESULTS_SUMMARY.*` describes the older 2016–2025 analysis.** Read the current numbers
  from `outputs/`.
- **Output names keep their history.** Files such as `py_statistical_results.xlsx` are produced by R;
  the `py_` prefix was kept so existing references still work. See `docs/R_PORT_NOTES.md`.
- Paths are all relative, so the folder works from any location.

More detail: `README.md` (method and what changed) and `docs/DATA_SOURCES.md` (data provenance).
