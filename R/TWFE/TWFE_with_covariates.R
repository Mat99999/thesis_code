################################################################################
# TWFE DiD (Full Covariates) - FINAL PRODUCTION RUN
# Fixes: Excluded WB 330, Removed V-Line, Forced Integer X-Axis (-1 to 3)
################################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(fixest, tidyverse, stringr, janitor, broom, here)

# --- 1. DEFINE PATHS (using here::here for portability) ---
data_path <- here::here("Data", "Archive_enriched")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
output_path <- here::here("output", "TWFE", "with_covariates")

# --- 2. CREATE DIRECTORIES ---
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(file.path(output_path, "graphs"))) dir.create(file.path(output_path, "graphs"), recursive = TRUE)
if (!dir.exists(file.path(output_path, "tables"))) dir.create(file.path(output_path, "tables"), recursive = TRUE)

# --- 3. LOAD LOOKUP & CONTROLS ---
sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE)
funder_lookup <- tibble(funder_abbr = c("ch", "wb"), funder_full = c("China", "World Bank"))

# The "Big 6" Covariates
controls <- c(
  "log_avg_pop_dens",
  "log_3yr_pre_conflict_deaths",
  "log_disasters",
  "election_year",
  "political_stability",
  "leader_birthplace"
)

# --- 4. MAIN LOOP (PRODUCTION) ---
data_files <- list.files(path = data_path, pattern = "InputData_.*\\.csv", full.names = TRUE)

for (file in data_files) {
  # A. Metadata Extraction
  filename <- basename(file)
  parts <- filename %>%
    str_remove("InputData_") %>%
    str_remove("_DiD.*\\.csv")
  funder_abbr <- str_extract(parts, "^[a-z]+")
  sector_code <- as.numeric(str_extract(parts, "[0-9]+"))

  # --- FIX 1: EXCLUDE WORLD BANK TRADE (330) ---
  if (funder_abbr == "wb" && sector_code == 330) {
    cat("\nSKIPPING: World Bank - Trade and Tourism (330) due to data constraints.")
    next
  }

  funder_full <- funder_lookup$funder_full[funder_lookup$funder_abbr == funder_abbr]
  sector_full <- sector_lookup$ad_sector_names[sector_lookup$ad_sector_codes == sector_code]
  if (length(sector_full) == 0) sector_full <- paste("Sector", sector_code)

  cat(paste0("\nProcessing: ", funder_full, " | ", sector_full, "... "))

  # B. Load & Clean Data
  df <- read_csv(file, show_col_types = FALSE) %>% clean_names()

  df <- df %>%
    mutate(time = as.numeric(time)) %>%
    group_by(dhs_id) %>%
    mutate(treatment_cohort = if_else(any(treated == 1), min(time[treated == 1]), Inf)) %>%
    ungroup() %>%
    mutate(relative_time = if_else(treatment_cohort != Inf, time - treatment_cohort, -Inf)) %>%
    mutate(relative_time = case_when(relative_time <= -5 ~ -5, relative_time >= 5 ~ 5, TRUE ~ relative_time)) %>%
    mutate(relative_time = as.integer(relative_time))

  # C. Dynamic Reference Finding
  available_pre_periods <- df %>%
    filter(relative_time < 0) %>%
    distinct(relative_time) %>%
    pull(relative_time) %>%
    sort()
  if (length(available_pre_periods) == 0) {
    cat("SKIPPING (No pre-periods).\n")
    next
  }
  ref_period <- max(available_pre_periods)

  # D. Run TWFE Model
  formula_str <- paste0("iwi_est_post_oda ~ i(relative_time, ref = ", ref_period, ") + ", paste(controls, collapse = " + "), " | dhs_id + time")
  twfe_model <- feols(as.formula(formula_str), data = df, cluster = ~dhs_id)

  # E. Extract Results & Visual Mapping
  # Mapping: -5 -> -1 | Ref -> 0 | 0 -> 1 | 3 -> 2 | 5 -> 3
  plot_data <- broom::tidy(twfe_model, conf.int = TRUE) %>%
    filter(str_detect(term, "relative_time::")) %>%
    mutate(relative_time = as.numeric(str_remove(term, "relative_time::"))) %>%
    add_row(relative_time = ref_period, estimate = 0, conf.low = 0, conf.high = 0) %>%
    mutate(visual_time = case_when(
      relative_time == -5 ~ -1, relative_time == ref_period ~ 0,
      relative_time == 0 ~ 1, relative_time == 3 ~ 2, relative_time == 5 ~ 3,
      TRUE ~ NA_real_
    )) %>%
    filter(!is.na(visual_time)) %>%
    arrange(visual_time)

  # F. Plotting (No vline, Forced Integers -1 to 3)
  p <- ggplot(plot_data, aes(x = visual_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "steelblue", alpha = 0.2) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 2.5) +
    scale_x_continuous(breaks = seq(-1, 3, 1), labels = c("-1", "0", "1", "2", "3"), limits = c(-1.2, 3.2)) +
    labs(
      x = "Event Time", y = "Effect on IWI", title = paste0(funder_full, ": ", sector_full),
      subtitle = "TWFE Event Study (Big 6 Covariates Included)"
    ) +
    theme_minimal()

  # G. Save Individual Results
  ggsave(filename = file.path(output_path, "graphs", paste0("TWFE_Covar_Plot_", funder_abbr, "_", sector_code, ".png")), plot = p, width = 8, height = 5)
  write_csv(plot_data, file = file.path(output_path, "tables", paste0("TWFE_Covar_Results_", funder_abbr, "_", sector_code, ".csv")))
  cat("Saved.")
}

# --- 5. COMPARISON ANALYSIS (OLD VS NEW) ---
path_old <- here::here("output", "TWFE", "without_covariates", "tables")
path_new <- here::here("output", "TWFE", "with_covariates", "tables")

files_old <- list.files(path_old, pattern = "\\.csv$", full.names = TRUE)
files_new <- list.files(path_new, pattern = "\\.csv$", full.names = TRUE)

get_id <- function(f) str_extract(basename(f), "[a-z]{2}_[0-9]{3}")

comparison_df <- map_dfr(files_old, function(f_old) {
  id <- get_id(f_old)
  if (id == "wb_330") {
    return(NULL)
  } # SKIP WB 330

  f_new <- files_new[str_detect(basename(files_new), id)]
  if (length(f_new) == 0) {
    return(NULL)
  }

  d_old <- read_csv(f_old, show_col_types = FALSE)
  d_new <- read_csv(f_new, show_col_types = FALSE)

  clean_data <- function(df, label) {
    if ("relative_time" %in% names(df)) {
      return(df %>% select(relative_time, estimate, conf.low, conf.high) %>% mutate(model = label))
    }
    if ("term" %in% names(df)) {
      return(df %>% filter(str_detect(term, "relative_time")) %>% mutate(relative_time = as.numeric(str_extract(term, "-?[0-9]+"))) %>% select(relative_time, estimate, conf.low, conf.high) %>% mutate(model = label))
    }
    return(NULL)
  }

  bind_rows(clean_data(d_old, "Old (No Covars)"), clean_data(d_new, "New (With Covars)")) %>% mutate(id = id)
})

# --- 6. PROCESS FINAL OUTPUTS ---
final_table <- comparison_df %>%
  separate(id, into = c("funder_abbr", "sector_code"), sep = "_", remove = FALSE) %>%
  mutate(sector_code = as.numeric(sector_code)) %>%
  left_join(sector_lookup %>% clean_names(), by = c("sector_code" = "ad_sector_codes")) %>%
  mutate(funder = ifelse(funder_abbr == "ch", "China", "World Bank"), panel_name = paste0(funder, " - ", ad_sector_names)) %>%
  filter(relative_time %in% c(1, 2, 3)) %>%
  distinct(panel_name, relative_time, model, .keep_all = TRUE) %>%
  pivot_wider(id_cols = c(panel_name, relative_time), names_from = model, values_from = estimate) %>%
  clean_names() %>%
  mutate(diff = new_with_covars - old_no_covars, change_type = case_when(sign(new_with_covars) != sign(old_no_covars) & abs(new_with_covars) > 0.05 ~ "FLIP", abs(new_with_covars) < abs(old_no_covars) ~ "Shrunk", TRUE ~ "Grew")) %>%
  arrange(panel_name, relative_time)

# --- 7. SAVE FINAL COMPARISON PLOT (UPDATED FOR INTEGER X-AXIS) ---
plot_data_final <- comparison_df %>%
  separate(id, into = c("funder_abbr", "sector_code"), sep = "_", remove = FALSE) %>%
  mutate(sector_code = as.numeric(sector_code)) %>%
  left_join(sector_lookup %>% clean_names(), by = c("sector_code" = "ad_sector_codes")) %>%
  mutate(funder = ifelse(funder_abbr == "ch", "China", "World Bank"), panel_name = paste0(funder, " - ", ad_sector_names)) %>%
  # Apply Visual Mapping to the Comparison Plot Data as well
  mutate(visual_time = case_when(
    relative_time == -5 ~ -1,
    relative_time %in% c(-4, -3, -2) ~ 0, # Assuming most common ref periods mapping to 0
    relative_time == 0 ~ 1,
    relative_time == 3 ~ 2,
    relative_time == 5 ~ 3,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(visual_time))

p_comp <- ggplot(plot_data_final, aes(x = visual_time, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, alpha = 0.4) +
  facet_wrap(~panel_name, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Old (No Covars)" = "gray60", "New (With Covars)" = "red")) +
  # FORCE INTEGER LABELS -1, 0, 1, 2, 3
  scale_x_continuous(breaks = seq(-1, 3, 1), labels = c("-1", "0", "1", "2", "3"), limits = c(-1.2, 3.2)) +
  labs(
    title = "Comparison: Baseline vs. Covariates",
    subtitle = "TWFE",
    x = "Event Time", y = "Estimate (IWI)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(output_path, "graphs", "Final_Comparison_Plot.png"), p_comp, width = 16, height = 12)
write_csv(final_table, file.path(output_path, "tables", "Final_Detailed_Comparison_Table.csv"))

cat("\n================================================================\n")
cat("   FINAL ANALYSIS COMPLETE (WB 330 EXCLUDED / NO V-LINES / INT X-AXIS)\n")
cat("================================================================\n")
options(tibble.print_max = Inf)
print(final_table)

# --- 8. PRE-TREND DIAGNOSTIC (PARALLEL TRENDS CHECK) ---

cat("\n================================================================\n")
cat("   PHASE 2: PARALLEL TRENDS DIAGNOSTIC (Old vs New)\n")
cat("================================================================\n")

# Path to the definitive old file
definitive_old_file <- here::here("output", "TWFE", "tables", "Table_TWFE_Parallel_Trends_Test.csv")

if (file.exists(definitive_old_file)) {
  # 1. Load Definitive Old Results
  old_pre <- read_csv(definitive_old_file, show_col_types = FALSE) %>%
    filter(relative_time == -5) %>%
    select(funder_abbr, sector_code, panel_name, estimate_old = estimate, sig_old = is_significant)

  # 2. Extract New Results from comparison_df (the internal object from Section 5)
  new_pre <- comparison_df %>%
    filter(model == "New (With Covars)", relative_time == -5) %>%
    # Parse ID to match columns in old_pre
    separate(id, into = c("funder_abbr", "sector_code"), sep = "_") %>%
    mutate(sector_code = as.numeric(sector_code)) %>%
    # Determine significance: if the 95% Confidence Interval does not cross zero
    mutate(sig_new = !(conf.low <= 0 & conf.high >= 0)) %>%
    select(funder_abbr, sector_code, estimate_new = estimate, sig_new)

  # 3. Merge and Classify Improvement
  pt_diagnostic <- inner_join(new_pre, old_pre, by = c("funder_abbr", "sector_code")) %>%
    mutate(
      improvement = case_when(
        sig_old == TRUE & sig_new == FALSE ~ "FIXED",
        sig_old == FALSE & sig_new == FALSE ~ "STABLE",
        sig_old == TRUE & sig_new == TRUE ~ "STILL BIASED",
        sig_old == FALSE & sig_new == TRUE ~ "NEW BIAS️"
      ),
      reduction_in_bias = abs(estimate_old) - abs(estimate_new)
    ) %>%
    select(panel_name, estimate_old, sig_old, estimate_new, sig_new, improvement, reduction_in_bias) %>%
    arrange(improvement)

  # 4. Save the table to your TWFECOV tables folder
  write_csv(pt_diagnostic, file.path(output_path, "tables", "Parallel_Trends_Diagnostic_Old_vs_New.csv"))

  # 5. Print Summary results to console
  print(pt_diagnostic)

  cat("\n--- SUMMARY OF PRE-TREND IMPROVEMENT ---\n")
  print(pt_diagnostic %>% count(improvement))
  cat(paste0("\n Diagnostic Table saved to: ", output_path, "/tables/Parallel_Trends_Diagnostic_Old_vs_New.csv\n"))
} else {
  cat("\n ERROR: Definitive old trends file not found at expected path.\n")
}
