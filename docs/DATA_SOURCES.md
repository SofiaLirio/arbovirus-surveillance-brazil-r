# Data sources and Source Data

Provenance for every raw input of the 2016–2026 arbovirus analysis, and an index of the
journal-style Source Data files. Raw files were re-downloaded from the official sources on
**2026-10-08** and checked against the archived copies the results were built from
(`data/updated_data.rar`, `data/metereologicos.rar`).

## 1. Raw sources

| Dataset | Publisher | Files | URL pattern |
|---|---|---|---|
| Arbovirus notifications (SINAN) | Ministério da Saúde, OpenDataSUS | `DENGBR16–26`, `CHIKBR16–26`, `ZIKABR16–26` `.csv.zip` (33) | `https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SINAN/{Dengue,Chikungunya,Zikavirus}/csv/{DENG,CHIK,ZIKA}BR{yy}.csv.zip` |
| Weather stations (hourly, automatic) | INMET, Banco de Dados Históricos | `2016.zip`–`2026.zip` (11) | `https://portal.inmet.gov.br/uploads/dadoshistoricos/{year}.zip` |
| Population estimates 2016–2021, 2024–2026 | IBGE, Estimativas da População | `estimativa_dou_*`, `estimativa_TCU_2016_*`, `POP20xx_*`, `serie_2001_20xx_TCU` (26) | `https://ftp.ibge.gov.br/Estimativas_de_Populacao/Estimativas_{year}/{file}` |
| Population 2022–2023 | IBGE, Projeções da População, Revisão 2024 | `projecoes_2024_tab1_idade_simples.xlsx`, `projecoes_2024_tab4_indicadores.xlsx` | `https://ftp.ibge.gov.br/Projecao_da_Populacao/Projecao_da_Populacao_2024/{file}` |
| State boundaries | IBGE, API de Malhas v3 | `uf_{code}.geojson` (27) | `https://servicodados.ibge.gov.br/api/v3/malhas/estados/{code}?formato=application/vnd.geo+json&qualidade=intermediaria` |

All of these are Brazilian federal open data, published under the Política de Dados Abertos
(Decreto 8.777/2016). Cite the publisher and the access date.

**Population.** `etl/build_population.R` builds the series from the raw IBGE files:
- **2016–2023:** Projeções 2024 revision, `projecoes_2024_tab4_indicadores.xlsx` column `POP_T`.
  IBGE published no estimates for 2022–2023 (census year). The earlier copied series turned out to
  be exactly this file for 2016–2023 (27/27 states every year).
- **2024–2026:** official DOU 1-July estimates, `estimativa_dou_<year>`. 2024–2025 equal the
  earlier series.
- **Not used:** the pre-census DOU estimates for 2016–2021. They sit 1–1.5% above the projections
  nationally and up to 11% above for single states (Amapá), and would break the series at 2022.

**Local copy:** `raw/` (99 files, about 1.7 GB) holds the exact files the results were built from.
It is not committed to the git repository; only `raw/MANIFEST.csv` is. SINAN 2026 is the **2026-10-06 release**, downloaded 2026-10-08 with
`tools/fetch_raw_sources.py` (now `tools/fetch_raw_sources.R`). `raw/MANIFEST.csv` lists bytes and sha256 for each file. Running
`fetch_raw_sources.R` again replaces the 2026 SINAN files with whatever release is current. After
that, `etl/check_db.R` reports differences for the 2026 case counts until
`reference/expected_hashes.csv` is regenerated.

## 2. Fresh downloads vs the archived inputs

| Result | Files |
|---|---|
| Byte-identical (CRC32) | 64 of 72: all INMET years, all IBGE files, and SINAN 2016–2024 |
| Re-zipped, same CSV content | `DENGBR25`, `CHIKBR25`, `ZIKABR25` (inner CSV CRC identical) |
| Newer data | `DENGBR26` 460,356 → 465,323 records; `CHIKBR26` 130,003 → 131,830; `ZIKABR26` 12,682 → 12,971 |

The 2026 SINAN files are preliminary, and OpenDataSUS republishes them weekly. **These results
use the release of 2026-10-06** (notifications through Sep 25–27, still no October). The earlier
2026-10-05 results used the 2026-09-26 release, which is archived in `updated_data.rar`. The
newer release adds about 1% more records. State the release date in the methods.

## 3. From raw files to `db/`

`etl/build_db.R` runs the whole build (population, cases + clinical tables, climate, checks).
The field decodes are in `etl/build_cases_clinical.R`. In short:
- **Cases:** counted by notification state (`SG_UF_NOT`) and month of `DT_NOTIFIC`.
  DENV counts probable cases (`CLASSI_FIN ≠ 5`). CHIKV and ZIKV count all notifications.
- **Climate:** INMET hourly data. −9999 values are masked. A station-month counts only with ≥ 240 valid
  temperature readings or ≥ 360 valid rain readings. State values are the mean across stations.
  When a year zip holds several cumulative snapshots of one station (`2026.zip` has files ending
  31-08-2026 and 30-09-2026 for each station), only the latest is used.
- **Incidence:** cases / population × 100,000.

## 4. Source Data (`source_data/`)

`pipelines/build_source_data.R` writes one workbook per figure, with one sheet per panel. Each
workbook has a `README` sheet naming the code each panel comes from. Numbers are recomputed with
the pipelines' own functions, then cross-checked against the saved results. All checks pass:
annual incidence, Mann-Kendall, Spearman, global Moran, mean incidence, LISA, and the national monthly series.

| Workbook | Contents |
|---|---|
| `Figure1` | A: national and regional pop-weighted annual incidence. B: box statistics and all 10,364 state-month rows |
| `Figure2` | A: Kendall τ per state. B: Sen's slope |
| `Figure3` | A: mean ρ and % significant by variable × lag. B: regional lag profile. Plus all state-level tests |
| `Figure4` | Moran scatter z / Wz per state, with I, p and the fitted line |
| `Figure5` | A: age-band counts and the plotted violin quartiles. B: sex counts and % |
| `Figure6` | Combined map: mean incidence and natural-breaks class per state × disease, symbol positions, class bounds and Moran's I |
| `Figure7`, `Figure8` | Clinical/laboratory counts and the plotted % (DENV; CHIKV and ZIKV) |
| `GLM_Figures` | glm_fig1–5 |
| `Additional_Figures` | DLNM, co-circulation, spatial regression, ITS, forecast, epidemic threshold |
| `Supplementary_Figures` | S1: incidence per state-year and its natural-breaks class (all years pooled), class bounds. S2: LISA and MK trend classes. py_02, py_03b, py_06b/c, py_07, clinical monthly and state totals, alt_results |
| `Tables` | Bootstrap CIs, global Moran, clinical breakdowns for 2016–2026 |

**Caveats:**
- **Figures 5, 7 and 8 end in 2025.** They are drawn from the GraphPad `.pzfx` tables. The 2016–2026
  breakdowns in `clinical_tables/` use different age bands.
- **The age violins are resampled.** They draw 6,000 synthetic ages per year from the binned counts
  (seed 42). The Source Data gives the bin counts plus the plotted quartiles.
- **2026 is annualised.** All 2026 annual incidence is annualised from January–September.

## 5. Reproduce

```bash
Rscript tools/fetch_raw_sources.R                      # re-download + manifest (resumable)
Rscript tools/compare_with_archive.R --data <folder>    # compare with updated_data.rar / metereologicos.rar
Rscript pipelines/build_source_data.R                   # source_data/*.xlsx + consistency checks
```
