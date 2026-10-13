# Arbovirus Surveillance in Brazil — Results Summary (Normalized DB pipeline)

> **Historical document.** This summary describes the original 2016–2025 analysis; its figure file
> names (e.g. `Figure6a`–`Figure6d`) refer to that version. For the current figures see the
> "Figures" section of the [README](../README.md) and [`FIGURE_CAPTIONS.md`](FIGURE_CAPTIONS.md).

**Run date:** 2026-10-01
**Data window:** 2016–2025 (CHIKV starts 2017; ZIKV and DENV start 2016)
**Spatial unit:** 27 Brazilian federative units
**Temporal unit:** monthly

All scripts in `pipelines/` were executed via `run_all.py` (9/9 stages succeeded).

**Inference conventions used throughout**

* Raw p-values are reported alongside **Benjamini–Hochberg FDR q-values** for every family of 27 state-level tests (Mann-Kendall, Spearman, LISA, per-state GLMs). Findings are only called robust if they survive FDR (q<0.05).
* GLMs include **month fixed effects**, so climate effects are measured against each state's usual seasonal cycle.
* GLM standard errors are **cluster-robust by state** (pooled models) or **HAC / Newey–West, 3 lags** (per-state monthly series).

---

## 1. Data scope

|Disease|Years|States|Monthly obs|Total cases|
|-|-|-|-|-|
|Dengue|2016–2025|27|3,240|16,347,701|
|Chikungunya|2017–2025|27|2,916|2,012,359|
|Zika|2016–2025|27|3,240|546,633|

Mean annual incidence per 100,000 inhabitants (mean over state-months × 12): **Dengue ≈ 624**, **Chikungunya ≈ 131**, **Zika ≈ 40**.

Missing monthly counts (4 CHIKV and 22 ZIKV state-months, at most 4 in any state-year) are not treated as zero: annual totals are annualised as `Total_Cases × 12 / Months_Reported` (`Months_Reported` is exported in the annual-incidence sheet).

---

## 2. Descriptive statistics — 95% CIs (normal-approx vs bootstrap)

Confidence intervals on the mean of monthly counts and incidence (pooled across all states × months), with 1,000-resample percentile bootstrap shown alongside the normal-approximation CIs. Both are computed by `bootstrap_ci.py` (`outputs/bootstrap_ci_95_*.csv`).

### Mean monthly cases

|Disease|Mean|Normal-approx 95% CI|Bootstrap 95% CI|n|
|-|-|-|-|-|
|DENV|5,045.6|(4,153.1 – 5,938.1)|(4,282.1 – 6,036.4)|3,240|
|CHIKV|691.1|(596.4 – 785.7)|(603.4 – 786.1)|2,912|
|ZIKV|169.9|(138.4 – 201.4)|(139.3 – 203.9)|3,218|

### Mean monthly incidence per 100,000

|Disease|Mean|Normal-approx 95% CI|Bootstrap 95% CI|n|
|-|-|-|-|-|
|DENV|51.99|(46.85 – 57.13)|(47.15 – 57.25)|3,240|
|CHIKV|10.91|(9.79 – 12.03)|(9.82 – 11.98)|2,912|
|ZIKV|3.34|(2.94 – 3.73)|(2.95 – 3.74)|3,218|

**Note.** For DENV, the bootstrap CI is shifted ~100 cases higher than the normal-approx CI, reflecting the right skew of the distribution (large epidemic peaks). CHIKV and ZIKV CIs are nearly identical between methods. Both methods treat state-months as independent, so they understate uncertainty for the national mean.

---

## 3. Temporal trend — Mann-Kendall on annual incidence per state

|Disease|Significant (raw p<0.05)|Increasing|Decreasing|Significant after FDR|Smallest q|
|-|-|-|-|-|-|
|Dengue|3 / 27|3|0|**0 / 27**|0.098|
|Chikungunya|3 / 27|3|0|**0 / 27**|0.132|
|Zika|3 / 27|0|3|**0 / 27**|0.286|

**States with raw p<0.05 (none survive FDR):**

* **Dengue ↑** Rio Grande do Sul (τ=0.73, p=0.004, q=0.098, Sen=+82.6/yr), Santa Catarina (τ=0.69, p=0.007, q=0.098, Sen=+108.8/yr), São Paulo (τ=0.51, p=0.049, q=0.44).
* **Chikungunya ↑** Santa Catarina (τ=0.78, q=0.13), Mato Grosso do Sul (τ=0.61, q=0.26), Rio Grande do Sul (τ=0.61, q=0.26).
* **Zika ↓** Alagoas, Distrito Federal, Piauí (all τ=−0.56, p=0.032, q=0.29).

**Interpretation.** Three raw-significant states out of 27 is about what chance alone produces at α=0.05 (~1.35 expected), and none survive FDR correction. With only 9–10 annual points per state, Mann-Kendall has little power. The southern Dengue/Chikungunya increases (RS, SC; q≈0.1) are the strongest candidates and are consistent with the spatial and GLM results, but they are **suggestive, not established**, trends.

---

## 4. Climate correlations — Spearman (lags 0–2 months)

|Disease|Variable|Tests|Significant (raw)|Significant (FDR)|Mean ρ|
|-|-|-|-|-|-|
|Dengue|Temp|81|40|35|+0.17|
|Dengue|Rain|81|10|0|−0.02|
|Chikungunya|Temp|81|43|35|+0.15|
|Chikungunya|Rain|81|19|3|+0.02|
|Zika|Temp|81|27|22|+0.07|
|Zika|Rain|81|14|1|+0.02|

**Interpretation.** The positive temperature correlations survive FDR in most cases (22–35 of 81 tests per disease); the rainfall correlations largely do not (0–3 of 81). These are unadjusted monthly correlations and therefore include the shared seasonal cycle — see Section 6 for the season-adjusted effect. Explicit 2-month rain lag (`rain[t] → incidence[t+2]`): 3–5 raw-significant states per disease, 0–2 after FDR — no national rain-lag signal.

---

## 5. Spatial autocorrelation

Both pipeline stages now use the same weights: KNN (k=5) on state polygon centroids from `db/brazil_states.geojson`, 999 permutations with a fixed seed.

### Global Moran's I

|Disease|I|p_sim|Significant|
|-|-|-|-|
|Dengue|**+0.52**|0.001|✅ strong positive clustering|
|Chikungunya|−0.07|0.42|❌ no spatial structure|
|Zika|−0.01|0.34|❌ no spatial structure|

### LISA clusters (Local Moran's I)

|Disease|Raw p<0.05|After FDR (q<0.05)|
|-|-|-|
|Dengue|8 HH (DF, ES, GO, MT, MS, MG, PR, SP); 4 LL (AP, MA, PA, RN); outliers AC (HL), RJ (LH)|**HH: Goiás, Minas Gerais · LL: Amapá · LH: Rio de Janeiro**|
|Chikungunya|1 HL outlier (Roraima)|none|
|Zika|HH Mato Grosso; LL Paraná; LH Pará|none|

**Interpretation.** Dengue has a clear high-incidence core in the Centro-Oeste/Sudeste (the global I is strongly significant, and GO/MG remain HH clusters after FDR). The Chikungunya and Zika local classifications do not survive multiple-testing correction and should not be interpreted.

---

## 6. Count-data regression (GLM) — pooled national models

Model: `Cases ~ Temp_Avg + Rain_mm + C(Month) + C(Location_Name) + offset(log Population)` — state and month fixed effects, so climate IRRs compare months that are warmer or wetter than usual *for that calendar month in that state*.

### Pooled IRRs, all 27 states, monthly grain

|Disease|Model|Temp IRR per °C (95% CI)|Rain IRR per mm (95% CI)|Dispersion|
|-|-|-|-|-|
|**Dengue**|Poisson (cluster SE)|2.155 (1.656–2.804)\*\*\*|1.005 (1.001–1.008)\*\*|φ = 8,376|
||Quasi-Poisson|2.155 (2.012–2.308)\*\*\*|1.005 (1.002–1.007)\*\*\*|φ = 8,376|
||**NegBin (cluster SE)**|**1.415 (1.175–1.703)\*\*\***|**1.006 (1.001–1.011)\***|α = 1.23|
|**Chikungunya**|Poisson (cluster SE)|1.355 (1.024–1.794)\*|1.008 (0.998–1.018) ns|φ = 1,784|
||Quasi-Poisson|1.355 (1.240–1.481)\*\*\*|1.008 (1.004–1.012)\*\*\*|φ = 1,784|
||**NegBin (cluster SE)**|**1.062 (0.798–1.414) ns**|1.001 (0.997–1.005) ns|α = 1.18|
|**Zika**|Poisson (cluster SE)|0.926 (0.717–1.198) ns|1.001 (0.994–1.008) ns|φ = 790|
||Quasi-Poisson|0.926 (0.813–1.056) ns|1.001 (0.996–1.006) ns|φ = 790|
||**NegBin (cluster SE)**|**0.971 (0.798–1.180) ns**|1.001 (0.997–1.004) ns|α = 1.19|

(\*\*\* p<0.001, \*\* p<0.01, \* p<0.05, ns = not significant)

AIC (pooled): NegBin 52,653 / 38,474 / 34,839 vs Poisson 19.1 M / 3.17 M / 1.19 M (Dengue / Chikungunya / Zika).

### Dispersion overview (per-state models)

|Family|Median dispersion|Range|
|-|-|-|
|Poisson|φ = 412.8|9.7 – 44,556|
|Quasi-Poisson|φ = 412.8|9.7 – 44,556|
|Negative Binomial|α = 0.78|0.29 – 2.41|

### Per-state models (NegBin, month FE, HAC SEs)

|Disease|Temp significant (raw / FDR)|Rain significant (raw / FDR)|
|-|-|-|
|Dengue|17 / 17 of 27|8 / 5 of 27|
|Chikungunya|9 / 7 of 27|6 / 5 of 27|
|Zika|10 / 6 of 27|9 / 2 of 27|

### Interpretation

* **Poisson is inappropriate**: residual dispersion is 10× to 44,000× the Poisson assumption, and its AIC is orders of magnitude worse than NegBin's. Even with cluster-robust SEs, its point estimates are dominated by the largest epidemic months (Dengue 2.16 vs NegBin 1.42).
* **Quasi-Poisson** shares Poisson's point estimates and its SEs still ignore within-state autocorrelation, so its CIs are too narrow; it is reported for comparison only.
* **Dengue** (NegBin, the preferred model): a month 1 °C warmer than usual for that state and calendar month is associated with **~1.4× the case rate (+41%, 95% CI +18% to +70%)**; +10 mm of rain → ~6% higher rate (p≈0.02). The temperature effect is significant per-state in 17/27 states after FDR.
* **Chikungunya and Zika**: once seasonality is removed and SEs account for clustering, **neither has a significant pooled temperature or rainfall effect**. The earlier positive Chikungunya and negative Zika temperature associations were driven by the shared seasonal cycle and by understated SEs.
* **Seasonality matters**: the raw Spearman temperature correlations (Section 4) partly reflect the common seasonal cycle; the month-FE models isolate anomalies.

### Model history

* An earlier version of this summary reported NegBin temperature IRRs of 2.10 (Dengue) and 1.32 (Chikungunya) from Nelder–Mead fits that never converged.
* NegBin α is now estimated by profile likelihood (as in R's `MASS::glm.nb`); every one of the 165 NegBin models converges, and the `converged` flag is in each diagnostics sheet.
* The current specification adds month fixed effects and robust SEs. Compared with the converged-but-unadjusted model, the Dengue estimate is essentially unchanged (1.41) while the Chikungunya (1.08→1.06) and Zika (0.89→0.97) effects lose significance.

---

## 7. Methodological notes

1. **Column rename:** the original `Annual_Incidence` column in the normalized DB was actually monthly incidence per 100 k inhabitants. Renamed to `Monthly_Inc_per100k` in `db/normalized_fact_cases.csv`. Values unchanged.
2. **CI 95%:** normal-approximation and bootstrap CIs are both computed by `bootstrap_ci.py` and exported to `outputs/bootstrap_ci_95_*.csv`.
3. **Year window:** all diseases truncated to ≥2016 for comparability; CHIKV naturally begins 2017.
4. **Seasonality and autocorrelation:** GLMs include month fixed effects; SEs are cluster-robust by state (pooled) or HAC with 3 lags (per-state). Quasi-Poisson uses Pearson-scaled SEs with t-distribution p-values/CIs (as R's `quasipoisson`).
5. **Multiple testing:** BH-FDR q-values within each family of 27 state tests (per virus; per virus × variable for Spearman; per virus × term for per-state GLMs).
6. **Spatial weights:** KNN k=5 on polygon centroids, shared by both spatial stages; permutation tests use a fixed seed.
7. **Missing months:** annual totals are annualised over reported months rather than treating missing months as zero.

---

## 8. Output files

### Generated artefacts (`outputs/`, and composite figures in `figures/`)

**Statistical tables**

|File|Contents|
|-|-|
|`outputs/py_statistical_results.xlsx`|Annual incidence (incl. `Months_Reported`), Mann-Kendall, Spearman, Global Moran, Rain Lag-2 — with `q_value` / `Significant_FDR`|
|`outputs/geo_spatial_stats.xlsx`|Mean incidence, MK, LISA clusters (raw and `LISA_Cat_FDR`), Global Moran|
|`outputs/ci_comparison.xlsx`|Normal-approx vs bootstrap CIs side-by-side|
|`outputs/glm_poisson_results.xlsx`|Poisson — 3 specs + diagnostics|
|`outputs/glm_quasi_poisson_results.xlsx`|Quasi-Poisson — 3 specs + diagnostics|
|`outputs/glm_negbin_results.xlsx`|Negative Binomial — 3 specs + diagnostics (incl. `converged`)|
|`outputs/glm_*_per_state.csv`|Per-state GLM coefficients incl. month effects (1,134 rows each)|
|`outputs/glm_*_lagged.csv`|Per-state GLMs with Temp/Rain lag-0/1/2 (1,458 rows each)|
|`outputs/glm_*_pooled_fe.csv`|Pooled national GLM with state + month fixed effects (120 rows each)|

Every GLM row carries a `cov_type` column (`cluster`, `HAC` or `quasi`).

**Plots** (`outputs/*.png`)

* `py_01`–`py_10` — main pipeline plots
* `geo_01`–`geo_06` — geospatial maps
* `Dengue_*`, `Chikungunya_*`, `Zika_*` — GraphPad-converted clinical/demographic plots

**Composite figures** (`figures/Figure*.png`) — publication-ready panels.
Figure 6 (maps) is available as four standalone figures (`Figure6a_choropleth_mean`,
`Figure6b_jenks_classification`, `Figure6c_temporal_maps`, `Figure6d_lisa_trend`) and as
`Figure6_combined` (panels A–D). GLM figures `figures/glm_fig1`–`glm_fig5` are produced by
`plot_glm_results.py`.

---

## 9. Quick takeaways

1. **Dengue dominates** by case volume (16.3 M cases), is the only arbovirus with significant **spatial clustering** (Moran's I=0.52; GO/MG high-incidence core after FDR), and has a robust temperature signal: **~1.4× the case rate per +1 °C above the seasonal norm** (95% CI 1.18–1.70).
2. **Chikungunya and Zika show no robust climate association** once seasonality and within-state correlation are accounted for.
3. **No state-level trend survives multiple-testing correction.** Southern increases in Dengue and Chikungunya (Rio Grande do Sul, Santa Catarina) are the strongest candidates (q≈0.1) but need a longer series to confirm.
4. **Negative Binomial is the appropriate model family** (Poisson over-dispersion 10×–44,000×); use the NegBin pooled and per-state outputs, with their robust SEs and FDR q-values, for inference.
