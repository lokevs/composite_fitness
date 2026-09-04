# The flaws of fitness functions in changing environments

Code and data pipeline for the manuscript *The flaws of fitness functions in
changing environments* (von Schmalensee et al., in revision). Running the
scripts below reproduces every figure and table of the paper.

Data and fitted models are archived on Zenodo: [10.5281/zenodo.22297476](https://doi.org/10.5281/zenodo.22297476)

## Requirements

- R 4.5.1 (base graphics only)
- Packages: `brms` (2.22.0) with `rstan` (2.32.7) and `bayestestR` for the
  thermal performance models, `maps`, `lubridate`

## Setup

```
git clone <repository>
cd <repository>
Rscript fetch_data.R
```

`fetch_data.R` downloads the nine data and model files (about 615 MB) from the
Zenodo record into `data/` and `models/`, verifying each against the MD5
checksum reported by Zenodo, and creates `output/`. Files already present with
the correct checksum are skipped. `Rscript fetch_data.R sandbox` uses the test
record on sandbox.zenodo.org instead.

## Run order

Run each script from the repository root, e.g. `Rscript R/comparisons.R`.
`R/functions.R` holds the shared functions and is sourced by the other scripts.

| Step | Script | Reads | Writes | Runtime |
|---|---|---|---|---|
| 1 | `R/thermal_performance_modelling.R` | `data/thermal_performance.txt`, `models/*.rds` | `data/parameter_estimates_*.txt`; Fig S2, S7 | 5 min (hours if the models must be refitted) |
| 2 | `R/temperature_data.R` | `data/era5_thermal_regimes.txt` | Fig 2 (map) | 1 min |
| 3 | `R/fitness_predictions.R` | steps 1-2 inputs, `data/ovi_rate_conversion.txt` | `data/simulated_fitness.txt`, `data/simulated_fitness_controlled_settings.txt`; Fig 2, 5, S4, S8, S9 | 85 min |
| 4 | `R/comparisons.R` | `data/simulated_fitness.txt`, `data/era5_thermal_regimes.txt` | `data/time_to_twofold_summary.txt`; Fig 3, 4, S3, S5, S6 | 1 min |
| 5 | `R/sensitivity_analyses.R` | `data/era5_thermal_regimes.txt` | `data/*_summary.txt`; Fig S10, S11 | 2 min |
| 6 | `R/conceptual_figures.R` | `data/thermal_performance.txt` | Fig 1, S1 | < 1 min |

The `brm()` calls in step 1 load the fitted models from `models/` when present
(as downloaded by `fetch_data.R`); without them the models are refitted with
the same seeds, which requires a working Stan toolchain.
The three `parameter_estimates_*.txt` files are derived from these models and
are included in the record only so that `R/functions.R`, which reads them when
sourced, can be loaded before step 1 recreates them.

## Figures

All figures are written to `output/` as multi-page PDFs named after the
manuscript figure they belong to (`Fig_<n>_<content>.pdf`); panels were
assembled from these pages.

| Figure | Files |
|---|---|
| 1 | `Fig_1_intro_figure_1` to `_6` |
| 2 | `Fig_2_fitness_components`, `Fig_2_composite_fitness`, `Fig_2_map`, `Fig_2_temperature_time_series`, `Fig_2_average_temperature_time_series` |
| 3 | `Fig_3_site_specific_relative_fitness_diff`, `Fig_3_regional_fitness_diff_map_ari/geo`, `Fig_3_regional_fitness_diff_density_ari/geo`, `Fig_3_time_to_twofold` |
| 4 | `Fig_4_example_sites_hourly` |
| 5 | `Fig_5_trait_heatmaps` |
| S1 | `Fig_S1_conceptual_figure_components`, `Fig_S1_conceptual_figure_predictions` |
| S2 | `Fig_S2_fitted_LRS/development/growth_reaction_norms` |
| S3 | `Fig_S3_fitness_diff_histograms` |
| S4 | `Fig_S4_fitness_diff_heatmaps` |
| S5 | `Fig_S5_growth_rate_diff` |
| S6 | `Fig_S6_site_specific_mean_fitness` |
| S7 | `Fig_S7_pp_check_LRS/dev/growth` |
| S8 | `Fig_S8_body_size_to_max_fec_conversion` |
| S9 | `Fig_S9_LRS_to_ovi_conversion` |
| S10 | `Fig_S10_lambda_approximation_validation`, `Fig_S10_fertility_window_sensitivity` |
| S11 | `Fig_S11_extrapolation_uncertainty` |

Three figures scatter points with unseeded random jitter and therefore differ
cosmetically between runs: `Fig_2_composite_fitness`,
`Fig_3_site_specific_relative_fitness_diff` and `Fig_S9_LRS_to_ovi_conversion`.
