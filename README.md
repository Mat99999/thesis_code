# Temporal Dynamics of Development Aid in Africa: Comparing China and World Bank Projects

**Author:** Mattias Antar | **University:** Linköping University |
**Supervisors:** Adel Daoud & Connor Jerzak | **Year:** 2025/2026

> This repository contains the complete replication code, data
> processing pipelines, and analytical datasets used to generate the
> findings for my Master's thesis.

## 1. About

### Context

This thesis investigates: To what extent do World Bank and Chinese aid
projects causally influence neighborhood-level wealth trajectories
across African countries? As China emerges as a major donor with a
distinct "non-interference" approach compared to the World Bank's
conditionality-based model, understanding their comparative
effectiveness at the local level is of importance to understand
development aid to the continent.

### Methodology

The study overcomes the scarcity of local economic data by
leveraging satellite-derived International Wealth Index (IWI) estimates, provided in tabular form at the DHS-cluster level, to construct a high-resolution longitudinal panel. To estimate causal
effects, the analysis employs a Difference-in-Differences (DiD)
framework. Crucially, it contrasts the conventional Two-Way Fixed
Effects (TWFE) model with the robust de Chaisemartin & d'Haultfoeuille
(dCdH) estimator to correct for biases introduced by staggered treatment
timing and heterogeneous effects.

## 2. Data Description

The analysis relies on panel datasets located in
the `Data/Archive_enriched/` directory. These files integrate high-precision geocoded aid data with satellite-based wealth measures and contextual covariates at the DHS-cluster level, where each DHS cluster is interpreted as representing a neighborhood.

### Input Data

-   **Treatment:** Geocoded aid projects from **AidData** (World Bank
    v1.4.2 and Global Chinese Development Finance v1.1.1). Only projects with precision codes 1–3 (exact locations, buffered locations, or administrative-level centroids) were included to ensure that treatment exposure can be meaningfully defined at the same local neighborhood scale as the satellite-derived wealth outcome.

-   **Outcome:** The **International Wealth Index (IWI)**, a continuous
    asset-based measure of household wealth (0-100) derived from
    satellite imagery using deep learning models trained on DHS survey
    data (Pettersson et al., 2023).

### Covariates

To ensure robustness, the analysis included the following covariates in separate DiD runs:

1.  **`log_avg_pop_dens`**: Log of average population density (Proxy for
    urbanization/demand).

2.  **`log_3yr_pre_conflict_deaths`**: Log of conflict-related deaths in the previous 3 years (UCDP data). 

3.  **`log_disasters`**: Log count of natural disasters.

4.  **`election_year`**: Binary indicator for national election years
    (1 if a national election occurred, 0 otherwise).

5.  **`political_stability`**: World Bank Governance Indicator
    (Perceptions of political stability/violence).

6.  **`leader_birthplace`**: Binary indicator (1 if the cluster is in the national leader's home region).

### Panel Structure

The data is split into 23 distinct CSV files, where each file represents
a unique Donor-Sector Panel. This allows for granular analysis of
specific aid types.

**File Naming Convention:** `InputData_{FUNDER}_{SECTOR_CODE}_DiD_enriched.csv`

- `InputData_ch_110_DiD_enriched.csv` → China (`ch`) aid in the Education (`110`) sector
- `InputData_wb_120_DiD_enriched.csv` → World Bank (`wb`) aid in the Health (`120`) sector

## 3. Repository Structure

The repository is organized into four main components that mirror the
analytical workflow: Input, Estimation, Robustness, and Synthesis.

```
thesis_code/
├── Data/
│   └── Archive_enriched/              # 23 enriched panel datasets
│       ├── InputData_{funder}_{sector}_DiD_enriched.csv (×23)
│       └── sector_group_names.csv     # Sector code lookup
├── R/
│   ├── TWFE/                          # Two-Way Fixed Effects models
│   │   ├── TWFE.R                     # Baseline (no covariates)
│   │   └── TWFE_with_covariates.R     # With 6 covariates
│   ├── dCdH/
│   │   └── dcdh_with_covariates.R     # Robust dCdH estimator (with covariates)
│   ├── Robustness/
│   │   └── Robustness.R               # Spatial spillover checks
│   └── Plots and Tables/
│       ├── THESISTABLES.R             # Generate HTML tables
│       ├── plotandtables.R            # Generate figures
│       └── Output/                    # Final HTML tables (Tables 1-4, A1-A7)
├── Stata/
│   ├── dCdH.do                        # Baseline dCdH (no covariates)
│   └── dcdh_results/
│       └── dcdh_combined_all_panels.csv  # Combined Stata results
└── output/
    ├── dCdH/
    │   └── dcdh_results.csv           # R dCdH with covariates (23 panels)
    ├── TWFE/
    │   ├── without_covariates/
    │   │   └── twfe_combined_all_panels.csv
    │   └── with_covariates/
    │       └── twfe_with_cov_combined_all_panels.csv
    └── spillover/                     # Robustness analysis output
```

### Key Output Files

| Object | Output File | Description |
|-------|-------------|-------------|
| dCdH  | `Stata/dcdh_results/dcdh_combined_all_panels.csv` | (Stata, no cov) |
| dCdH  | `output/dCdH/dcdh_results.csv` |  |
| TWFE (no cov) | `output/TWFE/without_covariates/twfe_combined_all_panels.csv` | (R, no cov) |
| TWFE (with cov) | `output/TWFE/with_covariates/twfe_with_cov_combined_all_panels.csv` | (R, with cov) |
| Tables | `R/Plots and Tables/Output/*.html` | Plots and Tables |

## 4. Methodology & Workflow

The analysis follows a reproducible pipeline, including robustness checks.

### Step 1: Baseline Estimation (Stata)

**Script:** `Stata/dCdH.do`

-   **Purpose:** Runs the dCdH DiD model without covariates.
-   **Details:** Loops through all panel datasets and applies the dCdH estimator without covariates.
-   **Output:** Combined results in `Stata/dcdh_results/dcdh_combined_all_panels.csv`

### Step 2: TWFE and dCdH with Covariates (R)

This step establishes the baseline for comparison and tests if results
hold when controlling for the 6 covariates.

-   **TWFE (Baseline & Covariates Adjusted):**
    -   `R/TWFE/TWFE.R`: Standard Two-Way Fixed Effects event-study model
    -   `R/TWFE/TWFE_with_covariates.R`: TWFE with 6 covariates

-   **dCdH (Robustness):**
    -   `R/dCdH/dcdh_with_covariates.R`: dCdH estimator with 6 covariates using `DIDmultiplegtDYN`

### Step 3: Spatial Robustness (R)

**Script:** `R/Robustness/Robustness.R`

-   **Purpose:** Provides additional diagnostic checks for spatial spillovers in the context of the thesis.
-   **Diagnostic Tests:**
    1.  **Distance-Based Dose-Response:** Tests if IWI gradients fade with distance from projects
    2.  **Exclusion Buffer:** Re-estimates effects while removing control units within 10-30km of treated sites
    3.  **Moran's I:** Tests for spatial autocorrelation in model residuals

### Step 4: Synthesis (R)

**Scripts:** `R/Plots and Tables/`

-   **`THESISTABLES.R`**: Uses `gt` package to create HTML summary tables (Tables 1-4, A1-A7)
-   **`plotandtables.R`**: Generates visualizations and tables
-   **Output Location:** `R/Plots and Tables/Output/`

## 5. Installation & Requirements

To replicate this analysis, you will need **R** and **Stata** installed.
The project uses the `here` package to manage relative paths.

### R Dependencies

```r
install.packages(c(
  "fixest", 
  "DIDmultiplegtDYN", 
  "tidyverse", 
  "janitor",
  "sf", 
  "spdep", 
  "rnaturalearth", 
  "gt", 
  "here", 
  "pacman"
))
```

### Stata Dependencies

```stata
ssc install did_multiplegt_dyn, replace
ssc install event_plot, replace
ssc install estout, replace
ssc install fs, replace
```

## 6. Pipeline Instructions

To reproduce the thesis results, recommended execution order:

1.  **Run dCdH Baseline (Stata):**
    -   Open `Stata/dCdH.do` and execute

2.  **Run TWFE & Diagnostics (R):**
    -   Run `R/TWFE/TWFE.R` (Baseline)
    -   Run `R/TWFE/TWFE_with_covariates.R`

3.  **Run dCdH with Covariates (R):**
    -   Run `R/dCdH/dcdh_with_covariates.R`

4.  **Run Robustness Checks (R):**
    -   Run `R/Robustness/Robustness.R`

5.  **Compile Final Outputs (R):**
    -   Run `R/Plots and Tables/THESISTABLES.R` to generate tables
    -   Run `R/Plots and Tables/plotandtables.R` to generate figures
