# Python → R port notes

The project was ported from Python to R on **2026-10-09**. Every script was rewritten in R
(same folders, same base names; `graphpad_to_python_pipeline.py` became `graphpad_pipeline.R`)
and the Python code, `requirements.txt` and `.venv/` were removed. The original Python version is
preserved in `arbovirus-brazil-portable.rar` next to this folder.

Output file and sheet names are unchanged, including the historical `py_` prefix
(`py_statistical_results.xlsx`, `py_01_annual_incidence.png`, …), so the papers, `source_data/`
and older notes keep pointing at the right files.

## How the port was checked

Every table the R code writes was compared, cell by cell, with the Python results it replaces.

| Output | Result |
|---|---|
| `db/*.csv`, `clinical_tables/clinical_demographic_long.csv` | **Byte-identical** (sha256 in `reference/expected_hashes.csv`; `etl/check_db.R` reports OK) |
| Annual incidence, Mann-Kendall (τ, S, p, Sen slope/intercept, trend), Spearman (ρ, p, q), rain lag-2, global Moran's I, mean incidence, Local I and quadrants | Identical (≤ 1e-9) |
| Co-circulation, COVID ITS, epidemic thresholds, `Dengue_00_general_data.csv` | Identical (≤ 1e-9; the CSV is byte-identical) |
| Poisson / quasi-Poisson GLMs | Identical to ~1e-10 (coefficients, SEs, p, AIC) |
| Negative-binomial GLMs | Within ~1e-6 (α is the optimum of a 1-D search with tolerance 1e-6) |
| Spatial regression | Coefficients within 1e-6, AIC / log-likelihood identical |
| DLNM curves | Within 2e-6 |
| SARIMA forecast | 95 of 99 series within 1e-4; 3 series differ by 1–3 % (see below) |
| Natural-breaks (Jenks) classes | Identical |
| Source Data workbooks | Same sheets, columns and rows; the 7 built-in consistency checks pass |

Figures are redrawn with ggplot2 to the same layout, sizes, colours and annotations. They are
equivalent, not pixel-identical.

## Expected differences

1. **Anything that uses random numbers** gives slightly different values because R's random-number
   generator is not numpy's (seeds are kept: 12345 for Moran/LISA, 42 elsewhere, so R reruns are
   reproducible):
   - permutation p-values and `z_sim` of global Moran's I (e.g. Chikungunya p 0.312 → 0.285);
   - LISA p-values; 2 of 81 LISA categories that sit right at p ≈ 0.05 change class;
   - residual Moran p-values in the spatial regression;
   - bootstrap confidence limits (≤ 1.2 %; means and normal CIs are identical);
   - the synthetic ages behind the age violins (quartiles differ by ≤ 1.4 years).
2. **DLNM confidence bands.** The point curves match, but the bands differ by up to ~10 %. The Python
   version built a rank-deficient basis (patsy `cr` plus intercept), statsmodels fell back to Newton
   and inverted a near-singular Hessian, so its standard errors were numerically unstable (some
   Zika points even had SE = 0). The R version uses a full-rank natural-spline basis spanning the
   same space (identical fit) and the observed-information covariance; its bands are the reliable ones.
3. **SARIMA forecasts.** R reproduces the statsmodels model exactly (same state-space form,
   approximate-diffuse start, 26-month likelihood burn-in, CSS starting values; the log-likelihood
   matches to all printed digits). The likelihood is very flat along the AR/MA ridge, so for 3 of 99
   series (Chikungunya–Acre, Zika–Sul, Zika–Sergipe) the two optimisers stop at slightly different
   points; MAE differs by 1.3–2.6 % there. National and regional forecasts match.

## Method mapping

| Python | R |
|---|---|
| pandas I/O, groupby | `data.table`; Kahan summation and numpy rounding reproduced in `R/utils.R` so the ETL is byte-identical |
| `to_csv(index=False)` | `write_csv_py()` in `R/utils.R` (CRLF, Python float repr, `True`/`False`) |
| xlsx (openpyxl) | `readxl` / `writexl` |
| `pymannkendall.original_test` | `mk_original_test()` in `R/stats.R` |
| `scipy.stats.spearmanr` | `spearman_test()` in `R/stats.R` |
| `esda.Moran`, `esda.Moran_Local` | `moran_global()`, `moran_local()` in `R/stats.R` |
| `libpysal` KNN weights | `spdep::knearneigh` (k = 5, planar centroids from `sf` with s2 off) |
| `mapclassify.NaturalBreaks` | `classInt::classIntervals(style = "jenks")` |
| statsmodels GLM (IRLS, HAC, cluster, NB2 profile likelihood) | `pipelines/glm_common.R`: the same IRLS stopping rule, observed-information sandwiches and Brent search over α |
| patsy `cr()` cross-basis (DLNM) | `splines::ns()` with patsy's knots (quantiles of the distinct exposure values) × linear lag basis |
| `spreg` OLS / ML_Lag / ML_Error | `lm`, `spdep::lm.RStests`, `spatialreg::lagsarlm` / `errorsarlm`; spreg's AIC conventions |
| statsmodels OLS HAC (ITS) | `sandwich::NeweyWest(prewhite = FALSE, adjust = FALSE)` |
| statsmodels SARIMAX | `pipelines/additional/adv_forecast.R` (KFAS filter, L-BFGS-B; series fitted in parallel) |
| matplotlib / seaborn / geopandas plots | `ggplot2`, `patchwork`, `sf`; composites with `grid` + `png` |
| `requests` downloads | `curl` |

## Behaviour kept on purpose

The port reproduces what the Python code did, including a few choices worth knowing about:

- `compute_annual` annualises over absent months (`Total × 12 / Months_Reported`), and months with
  no cases have no row, so a few low-incidence state-years are scaled up.
- Climate lags are positional (previous row of the same state), not calendar-based, across gaps.
- The co-circulation p-value is a naive t-test using the first series' length as n, and ignores
  the selection of the best of 13 lags.
- The spatial regression uses an unweighted mean of annual incidence, whereas Moran / LISA use a
  population-weighted mean.
- `cocirculation_heatmap.png` holds line and bar charts, not a heatmap; the ITS figure shows the
  raw series without the fitted line; the epidemic threshold has no endemic channel.

## Environment

R 4.6.1 with the CRAN packages in `R/install_packages.R` (versions in `R_SESSION_INFO.txt`).
Binary packages are used on Windows and macOS, so Rtools is not needed. On Linux, `sf` needs the
system libraries GDAL, GEOS, PROJ and udunits.
