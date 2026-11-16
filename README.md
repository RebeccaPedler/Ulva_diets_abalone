# Ulva_diets_abalone

# Does Ulva pay? A probabilistic cost–benefit analysis of macroalgal dietary supplementation in farmed greenlip abalone (*Haliotis laevigata* Donovan) 

**Manuscript in preparation**

---

## Repository structure

```
Ulva_diets_abalone/
├── code/
│   ├── 00_helpers.R                     # Helper function for posterior summary tables
│   ├── 01_clean_data.R                  # Data inspection, cleaning, ICC check and tank-level dataset
│   ├── 02_bayesian_weight_models.R      # Bayesian models of final (log) weight
│   ├── 03_bayesian_SGR_models.R         # Bayesian models of specific growth rate (SGR)
│   ├── 04_bayesian_mortality_models.R   # Effect of mortality on growth, and weekly mortality hazard models
│   ├── 05_economic_models.R             # Economic modelling (break-even Ulva meal price, days saved)
│   └── 06_bayesian_heavy_metal_models.R # Heavy metal concentrations by diet and tissue
├── data/
│   ├── individual_abalone_data.csv      # Individual-level abalone measurements
│   ├── per_capita_feed.csv              # Daily feed, water temperature, running tank counts and per-capita feed
│   ├── tank_data.csv                    # Tank-level dataset (12 rows), written by 01_clean_data.R
│   ├── heavy_metals.csv                 # Heavy metal concentrations in foot and viscera, by tank
│   └── metadata.csv                     # Data dictionary for all data files
├── models/                              # Saved brms model fits (.rds), written by scripts 01-04 and 06
│                                  
├── figures/                             # Figures from scripts 01-05 (data summaries, diagnostics, results)
├── outputs/                             # Results from 05_economic_models.R saved as csv files                       
└── README.md
```

---

## Data

### `individual_abalone_data.csv`
Individual-level records for 93,724 abalone across 12 tanks (4 tanks per diet: control, Ulva, wakame). Each row represents one abalone detected by drone imagery at trial end. Length is measured from the detected polygon and weight is derived from length using a fixed weight-length relationship. Tank-level variables (start and end counts, density, biomass, mortality, per-capita feed and diet cost) are repeated on every row for that tank. A full data dictionary is provided in `metadata.csv`.

### `per_capita_feed.csv`
Daily records of water temperature (minimum, maximum and average), feed delivered per tank (kg) and estimated live abalone counts for each tank, across the 131-day trial (10 November 2025 to 20 March 2026). Feed is delivered at a fixed daily ration per tank, with no feed delivered on weekly non-feeding days (these are recorded as `NA` in `feed_per_tank_kg`). Per-capita feed is derived each day as feed delivered divided by the estimated live count for that tank, and summed across the trial to produce the cumulative `per_capita_feed` value reported in `individual_abalone_data.csv` and `tank_data.csv`.

### `tank_data.csv`
Tank-level dataset (one row per tank) created by `01_clean_data.R` after filtering implausible values. Contains tank mean log weight, mean weight and mean length, starting conditions, mortality, per-capita feed and diet cost. This is the dataset used for the weight, SGR and mortality-growth models.

### `heavy_metals.csv`
Arsenic, cadmium, lead, mercury and iron concentrations (mg/kg) in foot muscle and viscera, with one sample of each tissue per tank. Values below the detection limit are recorded as `<0.01` (lead only) and are censored in `06_bayesian_heavy_metal_models.R`.

### `metadata.csv`
Data dictionary for all data files.

---

## Code

Scripts are numbered and are intended to be run in order. Scripts 01-04 and 06 each contain a `refit` switch near the top: set `refit <- TRUE` to refit models from scratch and overwrite the saved `.rds` files, or `refit <- FALSE` to load the saved fits from `models/`. Scripts 02-06 rely on the packages loaded in `01_clean_data.R` and on objects created earlier in the same R session, so they are best run in sequence within one session.

### `00_helpers.R`

Defines `posterior_summary_table()`and used by scripts 02, 03, 04 and 06.

### `01_clean_data.R`

Loads the individual-level and per-capita feed data, cleans and inspects, and compares starting conditions and moderators across diets. Fits an intercept-only model with a tank random effect to estimate the intra-class correlation, which supports aggregating to tank level. Writes `data/tank_data.csv`.

### `02_bayesian_weight_models.R`

Bayesian Gaussian models of tank mean log final weight, with diet and per-capita feed as predictors. Compares default, weakly informative and wide priors, and sensitivity models adding starting density or starting weight, using leave-one-out cross-validation. The primary model (`fit_weight_final.rds`) includes diet and per-capita feed. A hierarchical individual-level model with a tank random effect (`fit_hier.rds`) is fitted as a sensitivity check.

### `03_bayesian_SGR_models.R`

Calculates tank-level SGR (% per day) from tank mean log weight, starting weight and trial length, and fits Bayesian Gaussian models with diet, per-capita feed and starting weight. Includes prior sensitivity, a prior-only predictive check, and refits dropping tanks flagged by Pareto-k diagnostics. The primary model (`fit_sgr_final.rds`) is the growth input to `05_economic_models.R`.

### `04_bayesian_mortality_models.R`

Tests whether mortality influenced growth, then fits weekly mortality hazard models to assess the effects of water temperature, diet, starting density and starting weight. Produces the survival and water temperature figure.

### `05_economic_models.R`

R script for building break-even Ulva meal price distributions and a probability-of-saving curve across a range of Ulva meal prices, using the primary SGR model (`fit_sgr_final.rds`). Models days to harvest weight, feed cost and non-feed operating costs for control and Ulva diets, and tests sensitivity to harvest target weight and feed share of operating expenditure. This script loads the pre-fitted model from `models/`, so `03_bayesian_SGR_models.R` must be run first (with `refit <- TRUE`) to generate the required `.rds` file.

### `06_bayesian_heavy_metal_models.R`

Fits one Bayesian model per metal (arsenic, cadmium, lead, mercury, iron) to log concentration, with diet, tissue and their interaction, a tank random effect, and censoring for values below the detection limit. 

---

## Outputs

`outputs/` contains:

| File | Generated by |
|---|---|
| `ulva_wakame_premium_curve.csv` | `05_economic_models.R` — diet premium vs. algal meal price, Ulva and wakame (Part A) |
| `ulva_sgr_economic_summary.csv` | `05_economic_models.R` — SGR, days-saved and break-even summary at the actual meal price (Part B) |
| `ulva_sgr_saving_curve.csv` | `05_economic_models.R` — probability-of-positive-saving curve across Ulva meal price (Part C) |
| `probability_break_evens.csv` | `05_economic_models.R` — highest Ulva meal price at each profitability confidence level (Part C) |
| `days_by_target_df.csv` | `05_economic_models.R` — days-saved sensitivity to harvest target weight (Part D) |
| `ulva_breakeven_by_confidence.csv` | `05_economic_models.R` — break-even Ulva meal price at 50%, 95% and 99% confidence

---

## Citation

Citation details will be added upon submission. In the meantime, if you use these data or code, please contact the corresponding author.

## Contact

Rebecca Pedler
rebecca.pedler@yumbah.com
College of Science and Engineering, Flinders University
Yumbah Aquaculture, Port Lincoln SA
