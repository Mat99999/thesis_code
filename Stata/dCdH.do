************************************************************
* dCdH Analysis Script - Full 1990–2013 Enriched Archive
* Final Version - Uses relative paths via locals
************************************************************

clear all
clear mata
set more off

* --- Directories (using relative paths) -----------------------
* Set the base path relative to current working directory
* Users should set this to their thesis_code folder
local base_path "."
local data_dir "`base_path'/Data/Archive_enriched"
local output_dir "`base_path'/Stata/dcdh_results"
capture mkdir "`output_dir'"
capture mkdir "`output_dir'/graphs"
cd "`data_dir'"

* --- Required packages -----------------------------------
ssc install did_multiplegt_dyn, replace
ssc install event_plot, replace
ssc install fs, replace
ssc install estout, replace

* --- STEP 2: Sector name lookup ---------------------------
import delimited "`data_dir'/sector_group_names.csv", clear
forvalues i = 1/`=_N' {
    global sector_`=ad_sector_codes[`i']' "`=ad_sector_names[`i']'"
}
clear

* --- STEP 3: Identify all enriched CSVs -------------------
fs "InputData_*_enriched.csv"
local enriched_files `r(files)'
local nfiles : word count `enriched_files'

display "============================================================="
display "Found " `nfiles' " enriched panel files."
display "============================================================="

* --- STEP 4: Loop over all enriched input files -----------
foreach filename in `enriched_files' {

    * Extract funder/sector from filename
    local parts = subinstr("`filename'", "InputData_", "", 1)
    local parts = subinstr("`parts'", "_DiD_enriched.csv", "", 1)
    local funder_abbr = usubstr("`parts'", 1, strpos("`parts'", "_")-1)
    local sector_code = usubstr("`parts'", strpos("`parts'", "_")+1, .)

    local funder_full "China"
    if "`funder_abbr'" == "wb" local funder_full "World Bank"
    local sector_full "${sector_`sector_code'}"

    di ""
    di "-------------------------------------------------------------"
    di as text "Processing:", "`funder_full'", "Sector", "`sector_code'", "(`sector_full')"
    di "-------------------------------------------------------------"

    * --- STEP 4A: Import panel data -------------------------
    import delimited "`filename'", clear
    quietly duplicates drop dhs_id time, force

    * QC: required columns
    capture confirm variable dhs_id
    if _rc { di as error "Missing dhs_id in `filename' - skipping"; continue }
    capture confirm variable time
    if _rc { di as error "Missing time in `filename' - skipping"; continue }
    capture confirm variable treated
    if _rc { di as error "Missing treated in `filename' - skipping"; continue }
    capture confirm variable iwi_est_post_oda
    if _rc { di as error "Missing iwi_est_post_oda in `filename' - skipping"; continue }

    * Optional: keep rows where period == year_group
    capture confirm variable period
    local has_period = ( _rc == 0 )
    capture confirm variable year_group
    local has_yeargroup = ( _rc == 0 )
    if `has_period' & `has_yeargroup' {
        quietly keep if period == year_group
    }

    * Rename for clarity
    rename dhs_id id
    rename time t
    rename treated D
    rename iwi_est_post_oda Y

    * Ensure numeric time
    capture confirm numeric variable t
    if _rc { destring t, replace force }

    * --- STEP 4B: Declare panel structure -------------------
    xtset id t

    * --- STEP 4C: Run dynamic DiD (dCdH) -------------------
    did_multiplegt_dyn Y id t D, effects(4) placebo(4) cluster(id)

    * --- STEP 4D: Export coefficient table -----------------
    esttab using "`output_dir'/DiD_Table_`funder_abbr'_`sector_code'.csv", ///
        replace se r2 wide not noobs

    * --- STEP 4E: Event-study plot --------------------------
    local graph_title "Effect of Aid: `funder_full' - `sector_full' Sector"
    event_plot e(estimates)#e(variances), default_look ///
        graph_opt( xtitle("Event Time") ///
                   ytitle("Effect on IWI") ///
                   title("`graph_title'") ///
                   xlabel(-4(1)3) ) ///
        stub_lag(Effect_#) stub_lead(Placebo_#) together

    graph export "`output_dir'/graphs/DiD_Plot_`funder_abbr'_`sector_code'.png", replace

    * --- STEP 4F: QC Summary -------------------------------
    di "COMPLETED"
    di "-------------------------------------------------------------"
    
}

di ""
di "All analyses complete! Graphs and tables saved in `output_dir'/"
************************************************************
