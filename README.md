# Ipsilateral hip-knee arthroplasty multi-state analysis

Reproducible R code for the patient-side multi-state analysis of ipsilateral hip-knee arthroplasty order and revision risk developed for Lorenzo Carvisiglia's PhD thesis.

The repository contains analysis code only. Registry data are not distributed because they are subject to data-access and privacy restrictions. Users with appropriately structured source data can reproduce the workflow by providing their own local data files and, if necessary, adapting the file names and worksheet names in `config.R`.

## Analysis overview

The code implements the analysis from the first eligible ipsilateral primary arthroplasty (`P1`) and distinguishes:

- direct receipt of the second ipsilateral primary arthroplasty (`P2_direct`);
- revision before the second primary arthroplasty (`Rpre`);
- second primary arthroplasty after a previous revision;
- revision after direct `P2`;
- death as a competing event.

Transition-specific Cox models use a clock-reset time scale, age-stratified baseline hazards, patient-clustered robust standard errors, and time-varying implant-order effects where required. Model-based absolute transition probabilities are standardized to a common origin-specific covariate distribution. Additional code covers patient-level bootstrap uncertainty, overlap weighting, sensitivity analyses, joint-specific and implant-specific revision analyses, and validation checks.

## Requirements

A recent R installation is required. The scripts use CRAN packages including `dplyr`, `tidyr`, `stringr`, `purrr`, `readxl`, `lubridate`, `survival`, `broom`, `openxlsx`, `ggplot2`, `scales`, and `tibble`. Some scripts may use only a subset of these packages.

Install the main dependencies with:

```r
install.packages(c(
  "dplyr", "tidyr", "stringr", "purrr", "readxl", "lubridate",
  "survival", "broom", "openxlsx", "ggplot2", "scales", "tibble"
))
```

## Data setup

By default, the project expects three Excel workbooks under `data/raw/`:

```text
data/raw/
  hip.xlsx
  knee.xlsx
  registry.xlsx
```

The default worksheet names are defined in `config.R`:

- hip workbook: `Esporta foglio di lavoro`
- knee workbook: `Esporta foglio di lavoro`
- registry interventions: `INTERVENTI`
- registry deaths: `CON DECESSO COME EVENTO`

These defaults match the source files used for the thesis analysis but the paths and worksheet names can be changed without editing the analysis scripts. Either edit `config.R` locally or set the corresponding environment variables documented in that file.

The cleaning code also expects the original RIPO variable names used in the thesis data extracts. The raw data themselves are not included in this repository.

## Directory structure

```text
.
├── analysis/
│   ├── 01_clean_data.R
│   ├── 02_descriptive_tables.R
│   ├── 03_build_multistate_data.R
│   ├── 04_fit_transition_models.R
│   ├── 05_transition_probabilities.R
│   └── ...
├── data/
│   └── README.md
├── output/
├── config.R
├── config.example.R
├── run_core.R
└── run_full.R
```

`data/raw/` and generated files under `output/` are ignored by Git.

## Running the analysis

Run scripts from the repository root. For the main analysis:

```r
source("run_core.R")
```

This runs data cleaning, descriptive summaries, construction of the multi-state data, final Cox models, transition probabilities, the supporting `Rpre -> death` model check, and validation.

For the complete set of sensitivity, bootstrap, weighting, secondary, and validation analyses:

```r
source("run_full.R")
```

The full workflow includes 500 patient-level bootstrap replicates and can therefore take substantially longer than the core workflow. No parallelization is used.

## Script order

The numbered scripts are intended to be run in order. The most important dependencies are:

1. `01_clean_data.R` creates the cleaned patient-side cohort.
2. `03_build_multistate_data.R` creates the transition-specific analysis objects.
3. `04_fit_transition_models.R` fits the final transition-specific Cox models.
4. `05_transition_probabilities.R` estimates standardized origin-specific transition probabilities.
5. Later scripts use these saved objects for bootstrap inference, sensitivity analyses, overlap weighting, and secondary analyses.

The complete-pathway analysis combines the `P1 -> P2_direct` and `P2_direct -> Rpost` clock-reset components by convolution. The probability-update validation scripts compare the piecewise-exponential implementation with a product-integral calculation.

## Reproducibility notes

- All paths are relative to the repository root unless overridden in `config.R` or through environment variables.
- No user-specific absolute paths are required.
- No cluster, SLURM, task-array, or parallel-computing infrastructure is used.
- The administrative censoring date defaults to `2021-12-31`, matching the thesis analysis.
- Patient-level resampling is used in the bootstrap because a patient can contribute both sides.
- The propensity-score centre threshold is set to 100 patient-side units, matching the final thesis analysis.

## Data availability

The RIPO registry data used for the thesis cannot be redistributed through this repository. Researchers wishing to reproduce the numerical results need authorized access to the corresponding source data.

## Citation

If you use this code, please cite the associated thesis and/or publication when available. A `CITATION.cff` file is included for repository citation.
