# Archive_enriched Data Documentation

This folder (`Archive_enriched`) contains the core analytical datasets used for the analysis. Each CSV file represents a specific **Panel** defined by a unique combination of **Funder** (Donor) and **Aid Sector**. 

These datasets are "enriched" with a wide array of covariates including geolocational data, conflict indicators, political variables, and satellite-derived wealth estimates.

## 1. File Logic & Naming Convention

The datasets follow a strict naming convention:

`InputData_{FUNDER}_{SECTOR_CODE}_DiD_enriched.csv`

*   **Funder**: The source of the development finance.
    *   `ch`: China (Chinese Development Finance)
    *   `wb`: World Bank
*   **Sector Code**: A numeric code representing the specific sector of aid (e.g., Education, Health, Transport).
    *   These codes map to the definitions in `sector_group_names.csv`.
    *   *Example*: `110` corresponds to "Education".
*   **DiD**: Indicates the data is structured for Difference-in-Differences analysis.
*   **enriched**: Indicates the presence of additional pre-treatment and time-varying covariates.

**Example**: 
*   `InputData_ch_110_DiD_enriched.csv` contains data for analyzing the impact of **Chinese** aid in the **Education** sector.

### Sector Code Mapping (Reference)
*Refer to `sector_group_names.csv` for the complete list. Common codes include:*
*   `110`: Education
*   `120`: Health
*   `140`: Water Supply and Sanitation
*   `150`: Government and Civil Society
*   `160`: Other Social Infrastructure and Services
*   `210`: Transport and Storage
*   `230`: Energy Generation and Supply
*   `310`: Agriculture, Forestry and Fishing
*   ... and others.

## 2. Data Structure (How to Read)

Each CSV file is a standard comma-separated values file.
*   **Unit of Observation**: The dataset is at the **Cluster-Period** level. Each row represents a specific Demographic and Health Survey (DHS) cluster observed at a specific time point (or associated with a specific treatment window).
*   **Format**: The data serves a panel structure where clusters are tracked.
*   **Total Rows**: Varies by file, typically tens of thousands of rows depending on the number of available clusters and time periods covered.

## 3. Variable Dictionary

Below is a detailed breakdown of the variables present in the datasets.

### Identification & Location
*   **`dhs_id`** *(Numerical/Integer)*: The unique identifier for the DHS cluster. This allows for tracking the same geographical location across different data cuts.
*   **`country`** *(Character)*: The name of the country.
*   **`iso3`** *(Character)*: The 3-letter ISO country code (e.g., "ZAF" for South Africa).
*   **`adm2`** *(Character)*: The code or name for the second-level administrative division (e.g., District or Municipality).
*   **`lat`** *(Numerical/Float)*: Latitude of the cluster centroid.
*   **`lon`** *(Numerical/Float)*: Longitude of the cluster centroid.

### Treatment Variables (The "Panel Logic")
The core logic for the impact evaluation relies on these variables:
*   **`treated`** *(Numerical/Binary)*: The primary treatment indicator.
    *   `1`: The cluster received a project from the specific **Funder** in the specific **Sector** within the analysis window.
    *   `0`: The cluster did not receive such a project (serving as a control).
*   **`year_group`** *(Character)*: The specific time period or "wave" of the analysis (e.g., "2011:2013").
*   **`time`** *(Numerical)*: A broad time index (e.g., year) often used to align the panel.
*   **`log_treated_other_funder_n`** *(Numerical)*: Log count of projects received from the *other* funder (e.g., if the file is Chinese aid, this controls for World Bank aid).
*   **`log_other_sect_n`** *(Numerical)*: Log count of projects received in *other sectors* from the same funder (controls for bundled aid).
*   **`log_ch_loan_proj_n`** *(Numerical)*: Log count of Chinese loan projects (specific control for financing type).

### Outcome Variables (Wealth & Development)
*   **`iwi_est_post_oda`** *(Numerical/Float)*: **International Wealth Index (IWI)** estimate *after* the Official Development Assistance (ODA) period. This is often the primary outcome variable.
*   **`iwi_{START_YEAR}-{END_YEAR}`** (e.g., `iwi_1990-1992`, `iwi_2011-2013`) *(Numerical/Float)*: Historical IWI estimates for the cluster averaged over 3-year windows. These allow for examining pre-treatment trends (parallel trends assumption).
*   **`log_pc_nl_pre_oda`** *(Numerical/Float)*: Log of Nighttime Lights (NL) per capita *before* the aid period. A proxy for baseline economic activity.

### Geographic & Physical Covariates
*   **`log_avg_min_to_city`** *(Numerical/Float)*: Log of the average travel time (in minutes) to the nearest city. Proxy for market access/remoteness.
*   **`log_avg_pop_dens`** *(Numerical/Float)*: Log of average population density.
*   **`log_dist_km_to_gold`**, **`..._gems`**, **`..._dia`**, **`..._petro`** *(Numerical/Float)*: Log distance (in km) to various natural resources (Gold, Gemstones, Diamonds, Petroleum). Used to control for the "resource curse" or resource-driven development.

### Conflict & Risk Covariates
*   **`log_3yr_pre_conflict_deaths`** *(Numerical)*: Log count of conflict-related deaths in the 3 years prior to the observation. Measures baseline instability.
*   **`log_disasters`** *(Numerical)*: Log count (or intensity) of natural disasters.

### Political & Institutional Covariates
*   **`leader_birthplace`** *(Numerical/Binary)*: `1` if the cluster is in the birth region of the national leader, `0` otherwise. Controls for political favoritism.
*   **`election_year`** *(Numerical/Binary)*: `1` if an election occurred in the observation year.
*   **`unsc_aligned_us`** / **`unsc_non_aligned_us`** *(Numerical/Binary)*: Indicators of the country's voting alignment with the US in the UN Security Council. Proxy for geopolitical alignment.
*   **`country_gini`** *(Numerical/Float)*: The Gini coefficient for the country, measuring income inequality.

### Governance Indicators (World Bank WGI)
These are standard governance indicators, typically normally distributed indices (centered around 0, ranges ~ -2.5 to 2.5):
*   **`corruption_control`**: Control of Corruption.
*   **`gov_effectiveness`**: Government Effectiveness.
*   **`political_stability`**: Political Stability and Absence of Violence/Terrorism.
*   **`reg_quality`**: Regulatory Quality.
*   **`rule_of_law`**: Rule of Law.
*   **`voice_accountability`**: Voice and Accountability.

### Other
*   **`image_file_5k_3yr`** *(Character)*: Relative path to the satellite image tile associated with this cluster (used for deep learning models).
*   **`first_year_group`** *(Character/Numerical)*: The first year grouping observed for this unit.

---
**Note**: When reading these files for analysis, ensure your data pipeline treats missing values (NaN/NA) appropriately, as some historical wealth metrics or lag variables might be missing for specific years or disconnected clusters.
