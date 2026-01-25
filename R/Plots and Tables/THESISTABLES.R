######################### THESIS TABLES#################################
################################################################################
# Table 1: Parallel Trends Assumption (TWFE DiD)
# Author: Mattias Antar
# Purpose: Combine TWFE DiD results, display pre-treatment (lead) estimates,
#          test for parallel trend violations, and export academic-style table
################################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, janitor, gt, readr, stringr, here)

# ------------------------------------------------------------------------------
# 1. Paths (using here::here for portability)
# ------------------------------------------------------------------------------
twfe_path <- here::here("output", "TWFE", "without_covariates", "tables")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

# ------------------------------------------------------------------------------
# 2. Read all TWFE files
# ------------------------------------------------------------------------------
files <- list.files(twfe_path, pattern = "^twfe_results_.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "TWFE files.\n")

# Helper: read one file and extract donor + sector code
read_one_twfe <- function(f) {
  df <- suppressMessages(read_csv(f, show_col_types = FALSE))
  fname <- basename(f)
  m <- str_match(fname, "^twfe_results_([A-Za-z]+)_([0-9]+)\\.csv$")
  donor <- m[, 2]
  sector_code <- suppressWarnings(as.numeric(m[, 3]))

  df %>%
    mutate(
      funder_abbr = donor,
      sector_code = sector_code
    )
}

twfe_all <- purrr::map_dfr(files, read_one_twfe)

cat("\n[OK] Combined TWFE results preview:\n")
print(head(twfe_all, 10))

# ------------------------------------------------------------------------------
# 3. Read sector lookup and join
# ------------------------------------------------------------------------------
sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE) %>%
  clean_names()

if (!all(c("ad_sector_codes", "ad_sector_names") %in% names(sector_lookup))) {
  stop("[ERROR] Sector lookup missing expected columns.")
}

twfe_all <- twfe_all %>%
  left_join(sector_lookup, by = c("sector_code" = "ad_sector_codes")) %>%
  mutate(
    funder = case_when(
      funder_abbr == "wb" ~ "World Bank",
      funder_abbr == "ch" ~ "China",
      TRUE ~ funder_abbr
    )
  )

# ------------------------------------------------------------------------------
# 4. Keep only pre-treatment event time = -1 and exclude wb_330
# ------------------------------------------------------------------------------
parallel_df <- twfe_all %>%
  filter(event_time == -1, !(funder_abbr == "wb" & sector_code == 330)) %>%
  mutate(
    ci_95 = sprintf("[%.2f, %.2f]", conf.low, conf.high),
    estimate_fmt = sprintf("%.2f", estimate),
    # Violation test: CI excludes 0 ⇒ TRUE, else FALSE
    parallel_violation = if_else(conf.low > 0 | conf.high < 0, "Yes", "No")
  ) %>%
  select(funder, ad_sector_names, estimate_fmt, ci_95, parallel_violation) %>%
  arrange(funder, ad_sector_names)

cat("\n[OK] Parallel Trends Table Preview:\n")
print(parallel_df)

# ------------------------------------------------------------------------------
# 5. Format with gt (academic table)
# ------------------------------------------------------------------------------
parallel_table <- parallel_df %>%
  gt(groupname_col = "funder") %>%
  tab_header(
    title = md("**Table 1. Parallel Trends Assumption (TWFE DiD)**"),
    subtitle = md("Pre-treatment coefficients (event time = −1) with 95% confidence intervals.")
  ) %>%
  cols_label(
    ad_sector_names = "Sector",
    estimate_fmt = html("Estimate (β̂)"),
    ci_95 = html("[95% Confidence Interval]"),
    parallel_violation = html("Parallel Trend<br>Violation")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
    locations = cells_body()
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  tab_footnote(
    footnote = "A coefficient statistically different from zero indicates a potential violation of the parallel trends assumption.",
    locations = cells_title(groups = "subtitle")
  )

# ------------------------------------------------------------------------------
# 6. Save to THESISTABLES directory
# ------------------------------------------------------------------------------
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
output_file <- file.path(save_path, "Table1_ParallelTrends_TWFE.html")

gtsave(parallel_table, output_file)
cat("\n[OK] Table saved to:\n", output_file, "\n")

########### TABLE 2 ##################
################################################################################
# Table 2: Parallel Trends Assumption (dCdH DiD) — Clean Numeric Version (no stars)
# Author: Mattias Antar
# Purpose: Extract Placebo_1 coefficients from dCdH results, remove significance
#          stars, and render academic-style table identical to Table 1
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths (using here::here for portability)
# ------------------------------------------------------------------------------
dcdh_path <- here::here("Stata", "dcdh_results")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

# ------------------------------------------------------------------------------
# 2. Read all dCdH files
# ------------------------------------------------------------------------------
files <- list.files(dcdh_path, pattern = "^DiD_Table_.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "dCdH files.\n")

if (length(files) == 0) stop("[ERROR] No dCdH files found — check path or filenames.")

read_one_dcdh <- function(f) {
  df <- suppressMessages(read_csv(f, col_names = FALSE, show_col_types = FALSE))
  fname <- basename(f)
  m <- str_match(fname, "^DiD_Table_([A-Za-z]+)_([0-9]+)\\.csv$")
  donor <- m[, 2]
  sector_code <- suppressWarnings(as.numeric(m[, 3]))
  df <- df %>%
    filter(!is.na(X1)) %>%
    mutate(
      funder_abbr = donor,
      sector_code = sector_code
    )
  return(df)
}

dcdh_all <- purrr::map_dfr(files, read_one_dcdh)

cat("\n[OK] Combined dCdH raw preview:\n")
print(head(dcdh_all, 10))

# ------------------------------------------------------------------------------
# 3. Extract Placebo_1 rows
# ------------------------------------------------------------------------------
placebo_df <- dcdh_all %>%
  filter(str_detect(X1, "Placebo_1")) %>%
  transmute(
    funder_abbr,
    sector_code,
    coef_raw = X2,
    se_raw = X3
  )

cat("\n[OK] Placebo rows preview:\n")
print(head(placebo_df))

# ------------------------------------------------------------------------------
# 4. Parse numeric values (remove stars)
# ------------------------------------------------------------------------------
placebo_df <- placebo_df %>%
  mutate(
    # remove all asterisks before parsing
    coef_clean = str_remove_all(coef_raw, "\\*"),
    estimate = as.numeric(str_extract(coef_clean, "-?\\d+\\.?\\d*")),
    se = as.numeric(str_extract(se_raw, "\\d+\\.?\\d*")),
    conf.low = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se,
    ci_95 = sprintf("[%.2f, %.2f]", conf.low, conf.high),
    estimate_fmt = sprintf("%.3f", estimate)
  )

# ------------------------------------------------------------------------------
# 5. Join sector names
# ------------------------------------------------------------------------------
sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE) %>%
  clean_names()

if (!all(c("ad_sector_codes", "ad_sector_names") %in% names(sector_lookup))) {
  stop("[ERROR] Sector lookup missing required columns.")
}

placebo_df <- placebo_df %>%
  left_join(sector_lookup, by = c("sector_code" = "ad_sector_codes")) %>%
  mutate(
    funder = case_when(
      funder_abbr == "wb" ~ "World Bank",
      funder_abbr == "ch" ~ "China",
      TRUE ~ funder_abbr
    )
  )

# ------------------------------------------------------------------------------
# 6. Exclude WB Trade & Tourism (330)
# ------------------------------------------------------------------------------
placebo_df <- placebo_df %>%
  filter(!(funder_abbr == "wb" & sector_code == 330))

# ------------------------------------------------------------------------------
# 7. Violation flag
# ------------------------------------------------------------------------------
placebo_df <- placebo_df %>%
  mutate(
    parallel_violation = if_else(conf.low > 0 | conf.high < 0, "Yes", "No")
  ) %>%
  select(funder, ad_sector_names, estimate_fmt, ci_95, parallel_violation) %>%
  arrange(funder, ad_sector_names)

cat("\n[OK] Cleaned Parallel Trends (dCdH) summary preview:\n")
print(placebo_df)

# ------------------------------------------------------------------------------
# 8. Format as gt table
# ------------------------------------------------------------------------------
table2 <- placebo_df %>%
  gt(groupname_col = "funder") %>%
  tab_header(
    title = md("**Table 2. Parallel Trends Assumption (dCdH DiD)**"),
    subtitle = md("Placebo coefficients (pre-treatment) with 95% confidence intervals.")
  ) %>%
  cols_label(
    ad_sector_names = "Sector",
    estimate_fmt = html("Estimate (β̂)"),
    ci_95 = html("[95% Confidence Interval]"),
    parallel_violation = html("Parallel Trend<br>Violation")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
    locations = cells_body()
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  tab_footnote(
    footnote = "A coefficient statistically different from zero indicates a potential violation of the parallel trends assumption.",
    locations = cells_title(groups = "subtitle")
  )

# ------------------------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------------------------
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
output_file <- file.path(save_path, "Table2_ParallelTrends_dCdH.html")
gtsave(table2, output_file)
cat("\n[OK] Table saved to:\n", output_file, "\n")

############ TABLE A1 ##############
################################################################################
# Table A1: TWFE Event-Study Difference-in-Differences Results
# Author: Mattias Antar
# Purpose: Combine TWFE DiD event-study coefficients (post-treatment)
#          and display estimates with 95% confidence intervals in publication style
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths (using here::here for portability)
# ------------------------------------------------------------------------------
twfe_path <- here::here("output", "TWFE", "without_covariates", "tables")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

# ------------------------------------------------------------------------------
# 2. Read TWFE CSV files
# ------------------------------------------------------------------------------
files <- list.files(twfe_path, pattern = "^twfe_results_.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "TWFE files.\n")

if (length(files) == 0) stop("[ERROR] No TWFE files found — check path.")

read_one_twfe <- function(f) {
  df <- suppressMessages(read_csv(f, show_col_types = FALSE))
  fname <- basename(f)
  m <- str_match(fname, "^twfe_results_([A-Za-z]+)_([0-9]+)\\.csv$")
  donor <- m[, 2]
  sector_code <- suppressWarnings(as.numeric(m[, 3]))
  df %>%
    mutate(
      funder_abbr = donor,
      sector_code = sector_code
    )
}

twfe_all <- purrr::map_dfr(files, read_one_twfe)

cat("\n[OK] Combined TWFE results preview:\n")
print(head(twfe_all, 10))

# ------------------------------------------------------------------------------
# 3. Read sector lookup
# ------------------------------------------------------------------------------
sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE) %>%
  clean_names()

if (!all(c("ad_sector_codes", "ad_sector_names") %in% names(sector_lookup))) {
  stop("[ERROR] Sector lookup missing required columns.")
}

# ------------------------------------------------------------------------------
# 4. Keep post-treatment event times (1, 2, 3) and compute formatted CIs
# ------------------------------------------------------------------------------
twfe_all <- twfe_all %>%
  filter(event_time %in% c(1, 2, 3)) %>%
  left_join(sector_lookup, by = c("sector_code" = "ad_sector_codes")) %>%
  mutate(
    funder = case_when(
      funder_abbr == "wb" ~ "World Bank",
      funder_abbr == "ch" ~ "China",
      TRUE ~ funder_abbr
    ),
    estimate_fmt = sprintf("%.3f", estimate),
    ci_95 = sprintf("[%.2f, %.2f]", conf.low, conf.high)
  ) %>%
  filter(!(funder_abbr == "wb" & sector_code == 330)) %>%
  select(funder, ad_sector_names, event_time, estimate_fmt, ci_95) %>%
  arrange(funder, ad_sector_names, event_time)

cat("\n[OK] Post-treatment estimates preview:\n")
print(head(twfe_all))

# ------------------------------------------------------------------------------
# 5. Pivot to wide format for easier reading (t + 1 – t + 3)
# ------------------------------------------------------------------------------
twfe_wide <- twfe_all %>%
  pivot_wider(
    names_from = event_time,
    values_from = c(estimate_fmt, ci_95),
    names_glue = "t{event_time}_{.value}"
  ) %>%
  rename(
    `Sector` = ad_sector_names,
    `Estimate t + 1 (β̂)` = t1_estimate_fmt,
    `95% CI (t + 1)` = t1_ci_95,
    `Estimate t + 2 (β̂)` = t2_estimate_fmt,
    `95% CI (t + 2)` = t2_ci_95,
    `Estimate t + 3 (β̂)` = t3_estimate_fmt,
    `95% CI (t + 3)` = t3_ci_95
  )

cat("\n[OK] Wide-format preview:\n")
print(head(twfe_wide))

# ------------------------------------------------------------------------------
# 6. Create academic-style table
# ------------------------------------------------------------------------------
table_a1 <- twfe_wide %>%
  gt(groupname_col = "funder") %>%
  tab_header(
    title = md("**Table A1. Results from the TWFE Event-Study Difference-in-Differences Model**"),
    subtitle = md("Post-treatment average treatment effects (t + 1 to t + 3) with 95% confidence intervals.")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
    locations = cells_body()
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  cols_width(everything() ~ px(120)) %>%
  tab_footnote(
    footnote = "Each coefficient represents the estimated effect on local wealth (IWI) relative to the pre-treatment period.",
    locations = cells_title(groups = "subtitle")
  )

# ------------------------------------------------------------------------------
# 7. Save table
# ------------------------------------------------------------------------------
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
output_file <- file.path(save_path, "TableA1_TWFE_EventStudy.html")
gtsave(table_a1, output_file)
cat("\n[OK] Table A1 saved to:\n", output_file, "\n")

############# TABLE A2 #############
################################################################################
# Table A2: Results from the dCdH Difference-in-Differences Model (Fixed Parsing)
# Author: Mattias Antar
# Purpose: Correctly extract post-treatment effects (Effect_1–3) from Excel-style
#          encoded CSVs, compute 95% CIs, and render in publication style.
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths (using here::here for portability)
# ------------------------------------------------------------------------------
dcdh_path <- here::here("Stata", "dcdh_results")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

# ------------------------------------------------------------------------------
# 2. Read & clean all dCdH files
# ------------------------------------------------------------------------------
files <- list.files(dcdh_path, pattern = "^DiD_Table_.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "dCdH files.\n")
if (length(files) == 0) stop("[ERROR] No dCdH files found — check path or filenames.")

read_one_dcdh <- function(f) {
  df <- suppressMessages(read_csv(f, col_names = FALSE, show_col_types = FALSE))
  fname <- basename(f)
  m <- str_match(fname, "^DiD_Table_([A-Za-z]+)_([0-9]+)\\.csv$")
  donor <- m[, 2]
  sector_code <- suppressWarnings(as.numeric(m[, 3]))

  # Clean Excel-exported quote patterns like ="Effect_1"
  df <- df %>%
    mutate(across(X1:X3, ~ str_replace_all(., '="|"$', ""))) %>%
    mutate(across(X1:X3, ~ str_replace_all(., '"', ""))) %>%
    mutate(across(X1:X3, ~ str_trim(.))) %>%
    filter(!is.na(X1) & X1 != "") %>%
    mutate(funder_abbr = donor, sector_code = sector_code)
  return(df)
}

dcdh_all <- purrr::map_dfr(files, read_one_dcdh)
cat("\n[OK] Cleaned dCdH raw preview:\n")
print(head(dcdh_all, 10))

# ------------------------------------------------------------------------------
# 3. Extract post-treatment effects (Effect_1–Effect_3)
# ------------------------------------------------------------------------------
effects_df <- dcdh_all %>%
  filter(str_detect(X1, "^Effect_[1-3]$")) %>%
  transmute(
    term = X1,
    coef_raw = X2,
    se_raw = X3,
    funder_abbr,
    sector_code
  )

cat("\n[OK] Post-treatment effect rows preview:\n")
print(head(effects_df, 10))

# ------------------------------------------------------------------------------
# 4. Parse numeric values (remove stars, extract β and SE)
# ------------------------------------------------------------------------------
effects_df <- effects_df %>%
  mutate(
    coef_clean = str_remove_all(coef_raw, "\\*"),
    estimate = as.numeric(str_extract(coef_clean, "-?\\d+\\.?\\d*")),
    se = as.numeric(str_extract(se_raw, "\\d+\\.?\\d*")),
    conf.low = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se,
    ci_95 = sprintf("[%.2f, %.2f]", conf.low, conf.high),
    estimate_fmt = sprintf("%.3f", estimate),
    event_time = as.numeric(str_extract(term, "\\d+"))
  )

cat("\n[OK] Parsed numeric preview:\n")
print(head(effects_df))

# ------------------------------------------------------------------------------
# 5. Join sector lookup
# ------------------------------------------------------------------------------
sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE) %>%
  clean_names()

if (!all(c("ad_sector_codes", "ad_sector_names") %in% names(sector_lookup))) {
  stop("[ERROR] Sector lookup missing required columns.")
}

effects_df <- effects_df %>%
  left_join(sector_lookup, by = c("sector_code" = "ad_sector_codes")) %>%
  mutate(
    funder = case_when(
      funder_abbr == "wb" ~ "World Bank",
      funder_abbr == "ch" ~ "China",
      TRUE ~ funder_abbr
    )
  )

# ------------------------------------------------------------------------------
# 6. Exclude WB Trade & Tourism (330)
# ------------------------------------------------------------------------------
effects_df <- effects_df %>%
  filter(!(funder_abbr == "wb" & sector_code == 330))

# ------------------------------------------------------------------------------
# 7. Reshape to wide format (t + 1 to t + 3)
# ------------------------------------------------------------------------------
dcdh_wide <- effects_df %>%
  select(funder, ad_sector_names, event_time, estimate_fmt, ci_95) %>%
  pivot_wider(
    names_from = event_time,
    values_from = c(estimate_fmt, ci_95),
    names_glue = "t{event_time}_{.value}"
  ) %>%
  rename(
    `Sector` = ad_sector_names,
    `Estimate t + 1 (β̂)` = t1_estimate_fmt,
    `95% CI (t + 1)` = t1_ci_95,
    `Estimate t + 2 (β̂)` = t2_estimate_fmt,
    `95% CI (t + 2)` = t2_ci_95,
    `Estimate t + 3 (β̂)` = t3_estimate_fmt,
    `95% CI (t + 3)` = t3_ci_95
  ) %>%
  arrange(funder, `Sector`)

cat("\n[OK] Wide-format preview:\n")
print(head(dcdh_wide))

# ------------------------------------------------------------------------------
# 8. Format as academic-style gt table
# ------------------------------------------------------------------------------
table_a2 <- dcdh_wide %>%
  gt(groupname_col = "funder") %>%
  tab_header(
    title = md("**Table A2. Results from the dCdH Difference-in-Differences Model**"),
    subtitle = md("Post-treatment average treatment effects (t + 1 to t + 3) with 95% confidence intervals.")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
    locations = cells_body()
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  cols_width(everything() ~ px(120)) %>%
  tab_footnote(
    footnote = "Each coefficient represents the estimated effect on local wealth (IWI) relative to the pre-treatment period, based on the dCdH estimator.",
    locations = cells_title(groups = "subtitle")
  )

# ------------------------------------------------------------------------------
# 9. Save to HTML
# ------------------------------------------------------------------------------
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
output_file <- file.path(save_path, "TableA2_dCdH_EventStudy.html")
gtsave(table_a2, output_file)
cat("\n[OK] Table A2 saved to:\n", output_file, "\n")

############# TABLE A3 #######################

################################################################################
# Table A3: Distance-Based Dose–Response Robustness Test (OLS Model)
# Author: Mattias Antar
# Purpose: Display OLS coefficients for distance bands (dose–response analysis)
#          evaluating spatial spillovers near treated clusters.
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------
band_path <- here::here("output", "spillover", "Table1_band_tidy.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

# ------------------------------------------------------------------------------
# 2. Read and clean data
# ------------------------------------------------------------------------------
band_df <- read_csv(band_path, show_col_types = FALSE) %>%
  clean_names()

cat("\n[OK] Data preview:\n")
print(head(band_df, 10))

# ------------------------------------------------------------------------------
# 3. Keep distance-band coefficients only
# ------------------------------------------------------------------------------
band_filtered <- band_df %>%
  filter(str_detect(term, "dist_band")) %>%
  mutate(
    Distance_Band = str_remove_all(term, "dist_band|\\[|\\)"),
    Estimate = sprintf("%.3f", estimate),
    Std_Error = sprintf("%.3f", std_error),
    t_stat = sprintf("%.2f", statistic),
    p_value_fmt = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Estimate = paste0(Estimate, Significance)
  ) %>%
  select(Distance_Band, Estimate, Std_Error, t_stat, p_value_fmt)

cat("\n[OK] Distance-band coefficients preview:\n")
print(band_filtered)

# ------------------------------------------------------------------------------
# 4. Create gt table (academic style)
# ------------------------------------------------------------------------------
table_a3 <- band_filtered %>%
  gt() %>%
  tab_header(
    title = md("Table A3. Distance-Based Dose–Response Robustness Test (OLS Model)"),
    subtitle = md("Ordinary Least Squares (OLS) estimates comparing IWI levels across distance bands from treated clusters.")
  ) %>%
  cols_label(
    Distance_Band = md("Distance Band (km)"),
    Estimate = html("Estimate (β)*"),
    Std_Error = html("Std. Error"),
    t_stat = html("t Statistic"),
    p_value_fmt = html("p Value")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
    locations = cells_body()
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  tab_footnote(
    footnote = "Dependent variable: post-treatment IWI. The omitted reference category is clusters ≥50 km from treated areas. Significance levels: * p < 0.05, ** p < 0.01, *** p < 0.001.",
    locations = cells_title(groups = "subtitle")
  )

# ------------------------------------------------------------------------------
# 5. Save table
# ------------------------------------------------------------------------------
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
output_file <- file.path(save_path, "TableA3_DistanceBased_DoseResponse.html")
gtsave(table_a3, output_file)
cat("\n[OK] Table A3 saved to:\n", output_file, "\n")

########## TABLE A4 #############

################################################################################
# Table A4: Exclusion Buffer Robustness Test (OLS Model)
# Author: Mattias Antar
# Purpose: Display OLS treatment effect estimates with progressively larger
#          exclusion buffers (0–30 km) around treated clusters.
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------
buffer_path <- here::here("output", "spillover", "Table2_buffer_results.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

# ------------------------------------------------------------------------------
# 2. Read and clean data
# ------------------------------------------------------------------------------
buffer_df <- read_csv(buffer_path, show_col_types = FALSE) %>%
  clean_names()

cat("\n[OK] Data preview:\n")
print(buffer_df)

# ------------------------------------------------------------------------------
# 3. Format and prepare table
# ------------------------------------------------------------------------------
buffer_formatted <- buffer_df %>%
  mutate(
    k_buffer_km = as.character(k_buffer_km), # <- FIXED (use lowercase)
    Estimate = sprintf("%.3f", estimate),
    Std_Error = sprintf("%.3f", std_error),
    t_stat = sprintf("%.2f", statistic),
    p_value_fmt = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Estimate = paste0(Estimate, Significance)
  ) %>%
  select(
    `Exclusion Buffer (km)` = k_buffer_km,
    `Estimate (β̂)` = Estimate,
    `Std. Error` = Std_Error,
    `t Statistic` = t_stat,
    `p Value` = p_value_fmt
  )

cat("\n[OK] Formatted table preview:\n")
print(buffer_formatted)

# ------------------------------------------------------------------------------
# 4. Create gt table (academic style)
# ------------------------------------------------------------------------------
table_a4 <- buffer_formatted %>%
  gt() %>%
  tab_header(
    title = md("**Table A4. Exclusion Buffer Robustness Test**"),
    subtitle = md("Ordinary Least Squares (OLS) estimates excluding control clusters within increasing distances of treated sites.")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
    locations = cells_body()
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  tab_footnote(
    footnote = "Dependent variable: post-treatment International Wealth Index (IWI). Significance levels: * p < 0.05, ** p < 0.01, *** p < 0.001.",
    locations = cells_title(groups = "subtitle")
  )

# ------------------------------------------------------------------------------
# 5. Save table
# ------------------------------------------------------------------------------
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
output_file <- file.path(save_path, "TableA4_ExclusionBuffer_Test.html")
gtsave(table_a4, output_file)
cat("\n[OK] Table A4 saved to:\n", output_file, "\n")


###################### TABLE A5 #######################
################################################################################
# Table A5: Spatial Autocorrelation Diagnostic (Global Moran’s I)
# Author: Mattias Antar
# Purpose: Display the results from the spatial autocorrelation test of OLS
#          residuals using the Global Moran’s I statistic.
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------
moran_path <- here::here("output", "spillover", "Table3_Moran_summary.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

# ------------------------------------------------------------------------------
# 2. Read and clean data
# ------------------------------------------------------------------------------
moran_df <- read_csv(moran_path, show_col_types = FALSE) %>%
  clean_names()

cat("\n[OK] Data preview:\n")
print(moran_df)

# ------------------------------------------------------------------------------
# 3. Format values for presentation
# ------------------------------------------------------------------------------
moran_formatted <- moran_df %>%
  mutate(
    moran_i_fmt = sprintf("%.3f", moran_i),
    expected_fmt = sprintf("%.4f", expected),
    variance_fmt = sprintf("%.6f", variance),
    p_value_fmt = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) %>%
  select(
    `Global Moran’s I` = moran_i_fmt,
    `Expected Value`   = expected_fmt,
    `Variance`         = variance_fmt,
    `p Value`          = p_value_fmt
  )

cat("\n[OK] Formatted Moran's I table preview:\n")
print(moran_formatted)

# ------------------------------------------------------------------------------
# 4. Create gt table (academic style)
# ------------------------------------------------------------------------------
table_a5 <- moran_formatted %>%
  gt() %>%
  tab_header(
    title = md("**Table A5. Spatial Autocorrelation Diagnostic (Global Moran’s I)**"),
    subtitle = md("Test of spatial independence in OLS residuals from the post-treatment model.")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
    locations = cells_body()
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  tab_footnote(
    footnote = "A significant positive Moran’s I indicates that residuals exhibit spatial clustering, implying the presence of unobserved spatial dependence. Significance at p < 0.001 confirms non-random spatial autocorrelation.",
    locations = cells_title(groups = "subtitle")
  )

# ------------------------------------------------------------------------------
# 5. Save table
# ------------------------------------------------------------------------------
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
output_file <- file.path(save_path, "TableA5_SpatialAutocorrelation_Diagnostic.html")
gtsave(table_a5, output_file)
cat("\n[OK] Table A5 saved to:\n", output_file, "\n")

################################################################################
# Table A6: TWFE with Covariates Event-Study Difference-in-Differences Results
# Purpose: Generate a publication-ready horizontal table with robust error handling
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths (using here::here for portability)
# ------------------------------------------------------------------------------
twfe_path <- here::here("output", "TWFE", "with_covariates", "tables")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Read sector lookup
# ------------------------------------------------------------------------------
sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE) %>%
  clean_names()

# ------------------------------------------------------------------------------
# 3. Read & Combine TWFE files
# ------------------------------------------------------------------------------
files <- list.files(twfe_path, pattern = "^TWFE_Covar_Results_.*\\.csv$", full.names = TRUE)
cat("Found", length(files), "TWFE with covariates files.\n")

read_one_twfe_robust <- function(f) {
  # Clean names and unwrap Excel formatting
  df <- suppressMessages(read_csv(f, show_col_types = FALSE)) %>%
    remove_empty("cols") %>%
    clean_names() %>%
    mutate(across(everything(), ~ str_remove_all(.x, '="|"|\\)')))

  fname <- basename(f)
  m <- str_match(fname, "^TWFE_Covar_Results_([A-Za-z]+)_([0-9]+)\\.csv$")

  # Identify columns - use visual_time for event time in with_covariates files
  time_col <- names(df)[str_detect(names(df), "visual_time")][1]
  if (is.na(time_col)) time_col <- names(df)[str_detect(names(df), "time|term|x1")][1]
  est_col <- names(df)[str_detect(names(df), "^estimate$|est|x2")][1]
  low_col <- names(df)[str_detect(names(df), "conf_low|low|conf_l|x3")][1]
  high_col <- names(df)[str_detect(names(df), "conf_high|high|conf_h|x4")][1]

  df %>%
    mutate(
      time_num = as.numeric(.data[[time_col]]),
      funder_name = ifelse(m[, 2] == "ch", "China", "World Bank"),
      ad_sector_codes = as.numeric(m[, 3]),
      est_val = as.numeric(str_remove_all(.data[[est_col]], "\\*")),
      l_val = as.numeric(.data[[low_col]]),
      h_val = as.numeric(.data[[high_col]])
    ) %>%
    select(funder_name, ad_sector_codes, time = time_num, estimate = est_val, conf_low = l_val, conf_high = h_val)
}

twfe_all <- purrr::map_dfr(files, read_one_twfe_robust)

# ------------------------------------------------------------------------------
# 4. Preparation & Formatting (Exclude t0 and Sector Group)
# ------------------------------------------------------------------------------
twfe_final <- twfe_all %>%
  left_join(sector_lookup, by = "ad_sector_codes") %>%
  # Exclude t=0
  filter(time %in% c(-1, 1, 2, 3)) %>%
  mutate(
    time_label = paste0("t", ifelse(time > 0, "+", ""), time),
    cell_value = sprintf("%.3f\n(%.3f, %.3f)", estimate, conf_low, conf_high)
  ) %>%
  # Select only columns matching Table A1 style
  select(funder_name, ad_sector_names, time_label, cell_value) %>%
  pivot_wider(names_from = time_label, values_from = cell_value) %>%
  # Order columns chronologically
  select(funder_name, ad_sector_names, `t-1`, `t+1`, `t+2`, `t+3`) %>%
  arrange(funder_name, ad_sector_names)

# ------------------------------------------------------------------------------
# 5. Create the Publication-Style gt Table
# ------------------------------------------------------------------------------
table_a6 <- twfe_final %>%
  group_by(funder_name) %>%
  gt() %>%
  tab_header(
    title = md("**Table A6: TWFE with Covariates Event-Study Difference-in-Differences Results**"),
    subtitle = "Point Estimates and 95% Confidence Intervals"
  ) %>%
  cols_label(
    ad_sector_names = "Sector",
    `t-1` = "t-1",
    `t+1` = "t+1",
    `t+2` = "t+2",
    `t+3` = "t+3"
  ) %>%
  tab_options(
    table.font.size = px(11),
    row_group.background.color = "#f9f9f9",
    column_labels.font.weight = "bold",
    data_row.padding = px(3)
  ) %>%
  tab_style(
    style = cell_text(align = "center", size = px(10)),
    locations = cells_body(columns = c(`t-1`, `t+1`, `t+2`, `t+3`))
  ) %>%
  cols_width(
    ad_sector_names ~ px(250),
    everything() ~ px(120)
  ) %>%
  tab_source_note(
    source_note = md("*6 Covariates included: population density, conflict, natural disasters, election cycles, institutional stability, and leader birthplace.*")
  )

# ------------------------------------------------------------------------------
# 6. Save Table as HTML
# ------------------------------------------------------------------------------
gtsave(table_a6, file.path(save_path, "Table_A6_TWFE_Covariates.html"))

cat("\n[OK] Table A6 saved successfully as HTML in:", save_path, "\n")
print(table_a6)

################################################################################
# Table 4: Parallel Trends Assumption (dCdH DiD with Covariates)
# Author: Mattias Antar
# Purpose: Generate Table 4 showing pre-treatment (placebo) coefficients
################################################################################

# Packages already loaded at start of script

# ------------------------------------------------------------------------------
# 1. Paths (using here::here for portability)
# ------------------------------------------------------------------------------
dcdh_file <- here::here("output", "dCdH", "dcdh_results.csv")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
save_path <- here::here("R", "Plots and Tables", "Output")

if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Read dCdH Results and Sector Lookup
# ------------------------------------------------------------------------------
if (!file.exists(dcdh_file)) {
  cat("\n[WARNING] dCdH results file not found. Run dcdh_with_covariates.R first.\n")
  cat("Expected path:", dcdh_file, "\n")
} else {
  dcdh_raw <- read_csv(dcdh_file, show_col_types = FALSE) %>% clean_names()
  sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE) %>% clean_names()

  # Map sector codes to names
  sector_map <- setNames(sector_lookup$ad_sector_names, as.character(sector_lookup$ad_sector_codes))
  dcdh_raw$ad_sector_names <- sector_map[as.character(dcdh_raw$sector)]

  # Map funder codes to full names
  dcdh_raw$funder_name <- ifelse(dcdh_raw$funder == "ch", "China", "World Bank")

  # ------------------------------------------------------------------------------
  # 3. Create Table 4: Parallel Trends
  # ------------------------------------------------------------------------------
  placebo_data <- dcdh_raw %>%
    mutate(
      estimate_fmt = sprintf("%.2f", placebo),
      ci_95 = sprintf("[%.2f, %.2f]", placebo_ci_lower, placebo_ci_upper),
      parallel_violation = ifelse(
        (placebo_ci_lower > 0 & placebo_ci_upper > 0) | (placebo_ci_lower < 0 & placebo_ci_upper < 0),
        "Yes", "No"
      )
    ) %>%
    select(funder_name, ad_sector_names, estimate_fmt, ci_95, parallel_violation) %>%
    arrange(funder_name, ad_sector_names)

  table4 <- placebo_data %>%
    gt(groupname_col = "funder_name") %>%
    tab_header(
      title = md("**Table 4. Parallel Trends Assumption (dCdH DiD with Covariates)**"),
      subtitle = md("Placebo coefficients (pre-treatment) with 95% confidence intervals.")
    ) %>%
    cols_label(
      ad_sector_names = "Sector",
      estimate_fmt = "Estimate (β̂)",
      ci_95 = "[95% Confidence Interval]",
      parallel_violation = html("Parallel Trend<br>Violation")
    ) %>%
    cols_align(align = "left", columns = c(ad_sector_names, ci_95, parallel_violation)) %>%
    cols_align(align = "right", columns = estimate_fmt) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_style(
      style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
      locations = cells_body()
    ) %>%
    tab_footnote(
      footnote = "A coefficient statistically different from zero indicates a potential violation of the parallel trends assumption.",
      locations = cells_title(groups = "subtitle")
    )

  gtsave(table4, file.path(save_path, "Table4_ParallelTrends_dCdH_Covariates.html"))
  cat("\n[OK] Table 4 saved to:", save_path, "\n")

  # ------------------------------------------------------------------------------
  # 4. Create Table A7: Treatment Effects (t+1, t+2, t+3)
  # ------------------------------------------------------------------------------
  effects_data <- dcdh_raw %>%
    mutate(
      effect_1_fmt = sprintf("%.3f", effect_1),
      effect_2_fmt = sprintf("%.3f", effect_2),
      effect_3_fmt = sprintf("%.3f", effect_3),
      ci_1 = sprintf("[%.2f, %.2f]", ci_lower_1, ci_upper_1),
      ci_2 = sprintf("[%.2f, %.2f]", ci_lower_2, ci_upper_2),
      ci_3 = sprintf("[%.2f, %.2f]", ci_lower_3, ci_upper_3)
    ) %>%
    select(funder_name, ad_sector_names, effect_1_fmt, effect_2_fmt, effect_3_fmt, ci_1, ci_2, ci_3) %>%
    arrange(funder_name, ad_sector_names)

  table_a7 <- effects_data %>%
    gt(groupname_col = "funder_name") %>%
    tab_header(
      title = md("**Table A7. Results from the dCdH Difference-in-Differences Model with Covariates**"),
      subtitle = md("Post-treatment average treatment effects (t + 1 to t + 3) with 95% confidence intervals.")
    ) %>%
    cols_label(
      ad_sector_names = "Sector",
      effect_1_fmt = "Estimate t + 1 (β̂)",
      effect_2_fmt = "Estimate t + 2 (β̂)",
      effect_3_fmt = "Estimate t + 3 (β̂)",
      ci_1 = "95% CI (t + 1)",
      ci_2 = "95% CI (t + 2)",
      ci_3 = "95% CI (t + 3)"
    ) %>%
    cols_align(align = "left", columns = c(ad_sector_names, ci_1, ci_2, ci_3)) %>%
    cols_align(align = "right", columns = c(effect_1_fmt, effect_2_fmt, effect_3_fmt)) %>%
    cols_width(
      ad_sector_names ~ px(150),
      everything() ~ px(110)
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_style(
      style = cell_borders(sides = "bottom", color = "gray70", weight = px(1)),
      locations = cells_body()
    ) %>%
    tab_source_note(
      source_note = md("*6 Covariates included: population density, conflict, natural disasters, election cycles, institutional stability, and leader birthplace.*")
    ) %>%
    tab_footnote(
      footnote = "Each coefficient represents the estimated effect on local wealth (IWI) relative to the pre-treatment period.",
      locations = cells_title(groups = "subtitle")
    )

  gtsave(table_a7, file.path(save_path, "TableA7_dCdH_Covariates.html"))
  cat("\n[OK] Table A7 saved to:", save_path, "\n")
  print(table_a7)
}

cat("\n=== All thesis tables generated successfully! ===\n")
