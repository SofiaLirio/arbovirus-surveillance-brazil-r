# Arbovirus Surveillance in Brazil — 2016–2026 (updated, self-contained)

Reproducible statistical and geospatial analysis of **Dengue, Chikungunya, and Zika**
notifications across the **27 Brazilian federative units**, monthly, **2016–2026**.

A reproducible **R** project: the code that turns the raw official files into every table and
figure, the cleaned analysis tables (`db/`), and the generated results. Clone it anywhere and run
one script; no paths need editing. **New here? Start with [`GUIDE.md`](GUIDE.md).**

> **Raw data are not stored in this repository** (`raw/`, about 1.7 GB). On a fresh clone, run with
> `--skip-etl`: it rebuilds every result from the committed `db/`. To rebuild `db/` itself, see
> [Raw data](#raw-data) below.

## Quick start
The project is written in **R** (R ≥ 4.5; https://cloud.r-project.org). No Rtools or Python needed.
```powershell
# Windows
powershell -ExecutionPolicy Bypass -File run_all.ps1
```
```bash
# macOS / Linux / Git Bash
bash run_all.sh
```
```r
# or from any R session / RStudio
system("Rscript run_all.R")      # run_all.R accepts --skip-etl and --skip-install
```
The script installs any missing CRAN packages (`R/install_packages.R`, binary builds; the
versions used are recorded in `R_SESSION_INFO.txt`), then runs:

| Step | Command | Produces |
|---|---|---|
| 1 | `etl/build_db.R` | `db/` + `clinical_tables/` from `raw/` (checked against `reference/expected_hashes.csv`) |
| 2 | `pipelines/run_all.R` | 9 core stages → `outputs/`, `figures/` |
| 3 | `pipelines/additional/run_additional_analyses.R` | `outputs/2026-10-08/additional/` |
| 4 | `pipelines/alt_results.R` | alternative-view figures |
| 5 | `pipelines/build_source_data.R` | `source_data/` (one xlsx per figure) |

Flags: `-SkipEtl` / `--skip-etl` reuse an existing `db/`; `-SkipInstall` / `--skip-install`
skip the package check. Step 1 streams ~5.8 GB of uncompressed SINAN records and ~1.1 GB of
INMET zips (about 3 minutes); the whole run takes about 7 minutes.

The project was ported from Python to R on 2026-10-09. The rebuilt `db/` is byte-identical to the
Python build, and the statistics match the Python results; the few expected differences
(permutation p-values, bootstrap draws, …) are listed in [`docs/R_PORT_NOTES.md`](docs/R_PORT_NOTES.md).

- **Data window:** 2016–2026 (CHIKV now includes 2016; 2026 is season-to-date, Jan–Sep)
- **Spatial unit:** 27 federative units (UF) · **Temporal unit:** monthly
- **Generated:** 2026-10-08 (SINAN 2026 release of 2026-10-06)

The repository ships with the results already generated (`db/`, `clinical_tables/`, `outputs/`,
`figures/`, `source_data/`). `run_all.R --skip-etl` regenerates them from `db/`; the full run first
rebuilds `db/` from `raw/`. Output file names
keep their original names (e.g. `py_statistical_results.xlsx`, `py_01_*.png`) so references in
the papers and older notes stay valid.

## What changed vs the 2026-10-05 results (`arbovirus-surveillance-brazil-2026`)
| Input | Before | Now |
|---|---|---|
| SINAN 2026 | release 2026-09-26 (notifications to Sep 18–20) | release 2026-10-06 (to Sep 25–27) |
| 2026 cases | DENV 460,356 · CHIKV 130,003 · ZIKV 12,682 | DENV 465,323 · CHIKV 131,830 · ZIKV 12,971 |
| Population | copied series | parsed from raw IBGE files (identical values, see below) |
| INMET 2026 | each station counted twice for Jan–Aug (two snapshots in `2026.zip`) | latest snapshot per station only; 10/233 temperature and 16/232 rain cells changed (max 0.05 °C, 17 mm) |

SINAN 2016–2025, all climate before 2026, and population are unchanged, so every difference in
the results comes from the 2026 data.

---

## Folder layout
```
arbovirus-brazil-portable/
├── README.md                 # this file
├── run_all.R                 # one-command rebuild (packages + all steps)
├── run_all.ps1 / run_all.sh  # wrappers that find Rscript and run run_all.R
├── R/                        # shared R code: utils.R, stats.R, install_packages.R
├── R_SESSION_INFO.txt        # R + package versions of the last install
├── raw/                      # RAW official files (not in git; only MANIFEST.csv is committed)
│   ├── sinan/                # DENGBR/CHIKBR/ZIKABR 16–26 .csv.zip (OpenDataSUS)
│   ├── inmet/                # 2016.zip … 2026.zip hourly station data
│   ├── ibge/                 # population: DOU estimates + Projeções 2024
│   ├── geo/                  # IBGE state boundaries (uf_*.geojson)
│   └── MANIFEST.csv          # bytes + sha256 of the exact files the results were built from
├── reference/                # fixed inputs not rebuilt from raw
│   ├── dim_location.csv          # 27 UF ids/names
│   ├── brazil_states.geojson     # merged state polygons used by the maps
│   ├── graphpad_data/            # GraphPad clinical source (.pzfx, through 2025)
│   ├── original_db/              # original 2016–2025 series, for calibration printouts only
│   └── expected_hashes.csv       # sha256 of the db/ + clinical_tables/ these results came from
├── etl/                      # raw → db/ + clinical_tables/ (build_db.R runs it all)
├── pipelines/                # all analysis code
│   ├── run_all.R             # core 9-stage orchestrator
│   └── additional/           # advanced analyses (DLNM, forecast, spatial reg, ITS, …)
├── tools/                    # raw-source re-download + archive comparison (see docs/DATA_SOURCES.md)
├── docs/
│   ├── DATA_SOURCES.md           # raw-data URLs, access date, checksums, Source Data index
│   ├── RESULTS_SUMMARY.md/html   # ORIGINAL 2016–2025 narrative (see note below)
│   ├── README_original.md        # original project README (Python era)
│   ├── R_PORT_NOTES.md           # Python → R port: method mapping, verification, known differences
│   └── alt_results_README.md
│   (generated; rebuilt by run_all)
├── db/                       # normalized analysis inputs
├── clinical_tables/          # clinical/demographic breakdowns 2016–2026
├── outputs/                  # tables + plots (outputs/2026-10-08/additional/ = advanced analyses)
├── figures/                  # composite publication figures (Figure 1–8, S1–S2, GLM figures)
├── alt_results/              # 3 alternative-view figures
└── source_data/              # one xlsx per figure, one sheet per panel
```

## Raw data
`raw/` is not committed (about 1.7 GB of SINAN, INMET and IBGE downloads). `raw/MANIFEST.csv` is:
it lists the size and sha256 of the 99 files the results were built from; the source URLs are in
[`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md). To rebuild `db/` from raw, download them into `raw/`
and run the full pipeline:
```bash
Rscript tools/fetch_raw_sources.R     # all official sources (resumable); --only sinan,inmet,ibge,geo
Rscript run_all.R                     # ETL + everything; etl/check_db.R compares db/ with reference/expected_hashes.csv
```
OpenDataSUS republishes the 2026 SINAN files weekly, so a fresh download is newer than the
**2026-10-06 release** these results use and changes the 2026 counts (and the `check_db.R` hashes).
The 2026-10-08 re-download matched the archived inputs byte for byte for all INMET and IBGE files
and SINAN 2016–2024; the 2025 SINAN zips were re-packed with identical content (details in
`docs/DATA_SOURCES.md`).

## Figures
Main figures are in `figures/` (Figure 1–8), supplementary maps in `FigureS1_<year>.png`
(2016–2026, plus `FigureS1_overview.png`) and `FigureS2_lisa_trend.png`, and the numbers behind each
panel in `source_data/`. Every map uses one design: a single map of Brazil with one symbol per state
and disease (circle = dengue, triangle = chikungunya, square = Zika), filled by the mapped category.
Figures carry no in-plot titles. Text for the captions is generated with them:
- `outputs/figure_caption_values.csv`: statistics, colour keys and class bounds per figure panel;
- `outputs/figure_exclusions.csv` and `figure_exclusions_captions.txt`: records with ignored/blank
  values left out of each clinical figure (n and %), plus ready-to-use caption sentences;
- `outputs/figure_exclusions_state_years.csv`: state-years without SINAN records.

---

## How the data was built (provenance & definitions)

All inputs were rebuilt from the raw official sources and **calibrated to reproduce the original
series** (full detail in `clinical_tables/` companion notes and below).

**Cases** — SINAN notification microdata (DENGBR/CHIKBR/ZIKABR, 2016–2026). Counted by
notification state (`SG_UF_NOT`) and month of notification (`DT_NOTIFIC`). Case definition is
**disease-specific**, matching the original series:
- **Dengue = probable cases** (all notifications minus `Descartado`),
- **Chikungunya / Zika = all notifications**.

Validated against the original `db/` totals across 29 overlapping state-years: **25 within 1%,
28 within 5%** (only CHIKV 2020 differs more, −8.9%, a data-vintage difference). `Monthly_Inc_per100k
= Cases / Population × 100000`.

**Population** — parsed from the raw IBGE files by `etl/build_population.R`:
2016–2023 from the census-consistent *Projeções da População, Revisão 2024*
(`projecoes_2024_tab4_indicadores.xlsx`; IBGE published no estimates for 2022–2023), and
2024–2026 from the official 1-July estimates published in the DOU (`estimativa_dou_<year>`;
2026 Brazil total 214.2 M). The pre-census DOU estimates for 2016–2021 are not used: they overstate
several states (Amapá +11%) and would create a ~2.4 M drop between 2021 and 2022.

**Climate** — INMET automatic weather stations aggregated to state × month (−9999 sentinels
masked; coverage thresholds applied). `Temp_Avg` reproduces the original series well (corr 0.83).
⚠ **Rainfall** is on a different, physically-correct scale (~115 mm/month, realistic seasonality)
than the original `Rain_mm` (~36 mm/month) and is **not continuous** with it — the whole climate
series was rebuilt uniformly for all years, so climate-dependent results below are computed on
this consistent INMET series, **not** comparable in absolute rain scale to the original writeup.

**Clinical/demographic** (`clinical_tables/`) — tidy breakdowns (sex, age, classification,
criterion, evolution, hospitalization, serotype, lab results, pregnancy, UF totals) derived from
the same case population, 2016–2026.

---

## Reproducing the results
Use `run_all.ps1` / `run_all.sh` / `run_all.R` (see Quick start), or run the steps by hand from
any folder:
```bash
Rscript R/install_packages.R                         # once: CRAN packages
Rscript etl/build_db.R                               # raw -> db/ + clinical_tables/
Rscript pipelines/run_all.R                          # core 9 stages -> outputs/ + figures/
Rscript pipelines/additional/run_additional_analyses.R
Rscript pipelines/alt_results.R
Rscript pipelines/build_source_data.R
```
Every script finds the project root by itself, so they also run from RStudio (`source()` the file,
or open it and press *Source*). `etl/check_db.R` (run by `build_db.R`) confirms the rebuilt `db/`
CSVs are byte-identical to the ones the published results used. Permutation tests and the
bootstrap use fixed seeds, so reruns reproduce the numbers. Re-downloading with
`tools/fetch_raw_sources.R` pulls the newest 2026 SINAN release, which changes the 2026 counts
(see `docs/DATA_SOURCES.md`).

---

## Important notes
- **2026 is season-to-date** (cases Jan–Sep; climate Jan–Sep, late months partial; climate present
  for 26/27 states). Annual-incidence values for 2026 are annualised from reported months.
- **CHIKV 2016 is new** vs the original series (which began CHIKV in 2017).
- **`docs/RESULTS_SUMMARY.*` describes the ORIGINAL 2016–2025 analysis.** The figures/tables in
  `outputs/` and `figures/` here are the regenerated 2016–2026 results; where climate (esp. rain)
  is involved, re-read the numbers from these regenerated outputs rather than the old narrative.
- **Clinical *figures*** from the core pipeline (`outputs/Dengue_*`, etc.) are built from the
  GraphPad `.pzfx` source, which runs through **2025**; the **2016–2026 clinical *tables*** are in
  `clinical_tables/`.
