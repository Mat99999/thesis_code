################################################################################
# Plotting and tables for thesis
# Mattias Antar, 2025
################################################################################

# ---- Packages -----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, janitor, sf, rnaturalearth, rnaturalearthdata, ggplot2, viridis, here
)

# ---- Paths (using here::here for portability) ---------------------------------
data_path <- here::here("Data", "Archive_enriched")
plot_path <- here::here("output", "plots")
table_path <- here::here("output", "tables")
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Loading and merging input data ===\n")

# ---- Load and Merge All Input Files -------------------------------------------
files <- list.files(data_path, pattern = "InputData_.*_DiD_enriched\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
cat("Files found:", length(files), "\n")

read_one <- function(f) suppressMessages(readr::read_csv(f, show_col_types = FALSE)) %>% clean_names()
raw <- bind_rows(lapply(files, read_one))
cat("Rows loaded:", nrow(raw), "\n")

# ---- Clean and Prepare --------------------------------------------------------
raw <- raw %>%
  mutate(
    time = as.integer(time),
    iso3 = toupper(iso3),
    year_group = as.factor(year_group),
    cluster_id = paste0(iso3, "_", adm2)
  ) %>%
  filter(!is.na(lat), !is.na(lon))

cat("Unique clusters:", n_distinct(raw$cluster_id), "\n")
cat("Countries represented:", n_distinct(raw$iso3), "\n")

# ==============================================================================
# (1) Table & Plot: Number of Clusters by Year Cohort
# ==============================================================================
cluster_by_year <- raw %>%
  group_by(year_group) %>%
  summarise(n_clusters = n_distinct(cluster_id), .groups = "drop") %>%
  arrange(year_group) %>%
  mutate(total_clusters = sum(n_clusters))

write_csv(cluster_by_year, file.path(table_path, "Table1_ClusterCohorts.csv"))
cat("\n--- Cluster cohorts summary ---\n")
print(cluster_by_year)

# ==============================================================================
# (2) Table & Map: Cluster Coverage by Country
# ==============================================================================
cluster_by_country <- raw %>%
  group_by(country, iso3) %>%
  summarise(n_clusters = n_distinct(cluster_id), .groups = "drop") %>%
  arrange(desc(n_clusters))

# Add total rows
total_clusters <- sum(cluster_by_country$n_clusters, na.rm = TRUE)
n_countries <- nrow(cluster_by_country)
summary_tbl <- tibble(
  total_countries = n_countries,
  total_clusters = total_clusters
)

write_csv(summary_tbl, file.path(table_path, "Table3_ClusterSummary.csv"))
write_csv(cluster_by_country, file.path(table_path, "Table2_ClusterCoverage.csv"))

cat("\n--- Cluster coverage summary ---\n")
print(head(cluster_by_country, 10))
cat("\n--- Lowest coverage countries ---\n")
print(tail(cluster_by_country, 10))
cat("\nTotal clusters:", total_clusters, "across", n_countries, "countries\n")

# ---- Map Data -----------------------------------------------------------------
cat("\n=== Preparing map data ===\n")
africa_map <- rnaturalearth::ne_countries(continent = "africa", returnclass = "sf") %>%
  mutate(iso_a3 = ifelse(iso_a3 == "SDS", "SSD", iso_a3))

# Fix geometry issues
africa_map <- sf::st_make_valid(africa_map)

africa_df <- africa_map %>%
  left_join(cluster_by_country, by = c("iso_a3" = "iso3")) %>%
  mutate(
    log_clusters = log1p(n_clusters)
  )

missing <- africa_df %>%
  filter(is.na(n_clusters)) %>%
  pull(name_long)
cat("\nCountries in map without matching data:", length(missing), "\n")
if (length(missing) > 0) print(missing)

# Identify top 5 for labeling
top5 <- cluster_by_country %>% slice_max(n_clusters, n = 5)
cat("\nTop 5 countries by cluster count:\n")
print(top5)

# ---- Compute Label Points -----------------------------------------------------
label_points <- africa_df %>%
  st_point_on_surface() %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  filter(country %in% top5$country)

cat("\nQC — Label coordinates preview:\n")
print(label_points[, c("country", "X", "Y")])

# -----------------------------------------------------------------
# ---- Map Plot (No Country Labels) ---------------------------------------------
fig2 <- ggplot(africa_df) +
  geom_sf(aes(fill = log_clusters), color = "gray40", size = 0.2) +
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "lightgray",
    name = "Log(n Clusters)",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 15,
      barheight = 0.6
    )
  ) +
  coord_sf(default_crs = sf::st_crs(4326)) +
  labs(
    title = "Cluster Coverage Across Africa",
    subtitle = sprintf("Based on ADM2 units across 35 countries", total_clusters, n_countries),
    caption = "Grey areas = not inclued",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(size = 9, hjust = 0.5),
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

ggsave(file.path(plot_path, "Fig2_ClusterCoverage_Map.png"), fig2, width = 8, height = 6, dpi = 300)
cat("\n[OK] Saved Fig2_ClusterCoverage_Map.png (without labels)\n")

# ---- Summary QC Report --------------------------------------------------------
cat("\n=== Quality Control Summary ===\n")
cat("- Total clusters:", total_clusters, "\n")
cat("- Total countries:", n_countries, "\n")
cat("- Missing country joins:", length(missing), "\n")
cat("- Top country by clusters:", top5$country[1], "(", top5$n_clusters[1], "clusters)\n")
cat("- Files saved in:", plot_path, "and", table_path, "\n")

################################################################################
################################################################################
# Plot — Distribution of IWI Across ADM2 Units (2002–2013)
# Mattias Antar, 2025
################################################################################

library(tidyverse)

# ---- Paths (using here::here for portability) ---------------------------------
plot_path <- here::here("output", "plots")
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

# ---- Prepare Data --------------------------------------------------------------
df_iwi_subset <- raw %>%
  select(country, iso3, starts_with("iwi_")) %>%
  pivot_longer(
    cols = starts_with("iwi_"),
    names_to = "period",
    values_to = "iwi"
  ) %>%
  mutate(
    period = str_replace_all(period, "iwi_", ""),
    period = str_replace_all(period, "_", ":"),
    period = factor(period, levels = c(
      "2002:2004", "2005:2007", "2008:2010", "2011:2013"
    ))
  ) %>%
  filter(period %in% levels(period)) %>%
  drop_na(iwi)

# ---- QC Check ------------------------------------------------------------------
cat("Number of ADM2 observations per cohort:\n")
print(df_iwi_subset %>% count(period))

# ---- Density Plot --------------------------------------------------------------
fig_iwi_subset <- ggplot(df_iwi_subset, aes(x = iwi, fill = period)) +
  geom_density(alpha = 0.45, color = NA) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  labs(
    title = "Distribution of IWI Across ADM2 Units (2002–2013)",
    x = "IWI",
    y = "Probability Density",
    fill = "Cohort"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# ---- Save ----------------------------------------------------------------------
ggsave(
  filename = file.path(plot_path, "Fig_Distribution_IWI_2002_2013.png"),
  plot = fig_iwi_subset,
  width = 7, height = 5, dpi = 300
)

cat("\n[OK] Saved plot to:", file.path(plot_path, "Fig_Distribution_IWI_2002_2013.png"), "\n")

#######################################################
################################################################################
# Summary Statistics — IWI Across ADM2 Units (2002–2013)
# Mattias Antar, 2025
# Output: Table_IWI_SummaryStats_2002_2013.csv
################################################################################

library(tidyverse)
library(e1071) # for skewness

# ---- Paths (using here::here for portability) ---------------------------------
table_path <- here::here("output", "tables")
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)

# ---- Prepare Data --------------------------------------------------------------
df_iwi_subset <- raw %>%
  select(country, iso3, starts_with("iwi_")) %>%
  pivot_longer(
    cols = starts_with("iwi_"),
    names_to = "period",
    values_to = "iwi"
  ) %>%
  mutate(
    period = str_replace_all(period, "iwi_", ""),
    period = str_replace_all(period, "_", ":"),
    period = factor(period, levels = c(
      "2002:2004", "2005:2007", "2008:2010", "2011:2013"
    ))
  ) %>%
  filter(period %in% levels(period)) %>%
  drop_na(iwi)

# ---- QC Check ------------------------------------------------------------------
cat("Number of ADM2 observations per cohort:\n")
print(df_iwi_subset %>% count(period))

# ---- Compute Descriptive Statistics --------------------------------------------
iwi_summary_stats <- df_iwi_subset %>%
  group_by(period) %>%
  summarise(
    n_obs = n(),
    mean_iwi = mean(iwi, na.rm = TRUE),
    median_iwi = median(iwi, na.rm = TRUE),
    sd_iwi = sd(iwi, na.rm = TRUE),
    skewness_iwi = e1071::skewness(iwi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_iwi = round(mean_iwi, 2),
    median_iwi = round(median_iwi, 2),
    sd_iwi = round(sd_iwi, 2),
    skewness_iwi = round(skewness_iwi, 2)
  )

# ---- Console Preview -----------------------------------------------------------
cat("\n--- Descriptive Statistics for IWI (2002–2013) ---\n")
print(iwi_summary_stats)

# ---- Save to CSV ---------------------------------------------------------------
output_file <- file.path(table_path, "Table_IWI_SummaryStats_2002_2013.csv")
write_csv(iwi_summary_stats, output_file)

cat("\n[OK] Saved descriptive IWI summary table to:\n", output_file, "\n")

################################################################################
### HÄR NU
############
################################################################################
# Final Harmonized Comparison — TWFE vs dCdH
# Mattias Antar, 2025
################################################################################

library(tidyverse)
library(janitor)
library(readr)
library(stringr)

# ---- Paths (using here::here for portability) ---------------------------------
path_twfe <- here::here("output", "TWFE", "tables")
path_dcdh <- here::here("Stata", "dcdh_results")
sector_map_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
plot_path <- here::here("output", "plots")
table_path <- here::here("output", "tables")

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)

# ---- Helper: Extract Funder + Sector Code -------------------------------------
extract_meta <- function(filename) {
  funder <- case_when(
    str_detect(filename, "ch") ~ "China",
    str_detect(filename, "wb") ~ "World Bank",
    TRUE ~ NA_character_
  )
  sector_code <- as.integer(str_extract(filename, "([0-9]{3})(?=\\.csv$)"))
  tibble(funder = funder, ad_sector_codes = sector_code)
}

# ---- Readers ------------------------------------------------------------------
read_twfe <- function(f) {
  meta <- extract_meta(f)
  df <- suppressMessages(read_csv(f, show_col_types = FALSE)) %>% clean_names()
  df %>%
    rename(
      event_time = event_time,
      estimate   = estimate,
      ci_lower   = conf_low,
      ci_upper   = conf_high
    ) %>%
    mutate(
      funder = meta$funder,
      ad_sector_codes = meta$ad_sector_codes,
      method = "TWFE",
      file_name = basename(f)
    ) %>%
    mutate(across(c(event_time, estimate, ci_lower, ci_upper), as.numeric))
}

read_dcdh <- function(f) {
  meta <- extract_meta(f)
  df <- suppressMessages(read_csv(f, col_names = FALSE, show_col_types = FALSE))
  names(df) <- c("label", "estimate_raw", "se_raw")
  df %>%
    filter(str_detect(label, "Effect|Placebo")) %>%
    mutate(
      label = str_remove_all(label, '=\"|\"|='),
      estimate = as.numeric(str_extract(estimate_raw, "-?[0-9.]+")),
      se = as.numeric(str_extract(se_raw, "[0-9.]+")),
      ci_lower = estimate - 1.96 * se,
      ci_upper = estimate + 1.96 * se,
      event_time = case_when(
        str_detect(label, "Placebo_1") ~ -1,
        str_detect(label, "Effect_1") ~ 1,
        str_detect(label, "Effect_2") ~ 2,
        str_detect(label, "Effect_3") ~ 3,
        TRUE ~ NA_real_
      ),
      funder = meta$funder,
      ad_sector_codes = meta$ad_sector_codes,
      method = "dCdH",
      file_name = basename(f)
    ) %>%
    filter(!is.na(event_time))
}

# ---- Read & Combine -----------------------------------------------------------
files_twfe <- list.files(path_twfe, pattern = "\\.csv$", full.names = TRUE)
files_dcdh <- list.files(path_dcdh, pattern = "\\.csv$", full.names = TRUE)

twfe_all <- bind_rows(lapply(files_twfe, read_twfe))
dcdh_all <- bind_rows(lapply(files_dcdh, read_dcdh))

cat("\nQC — Files Loaded:\n")
cat("- TWFE:", nrow(twfe_all), "rows across", length(files_twfe), "files\n")
cat("- dCdH:", nrow(dcdh_all), "rows across", length(files_dcdh), "files\n")

combined <- bind_rows(twfe_all, dcdh_all) %>%
  filter(!is.na(event_time), !is.na(estimate))

sector_map <- read_csv(sector_map_path, show_col_types = FALSE) %>% clean_names()
combined <- combined %>%
  left_join(sector_map, by = "ad_sector_codes") %>%
  mutate(sector_label = paste0(ad_sector_names, " — ", funder)) %>%
  filter(sector_label != "Trade and Tourism — World Bank")

# ---- QC -----------------------------------------------------------------------
cat("\nQC — Combined Dataset Summary\n")
cat("Unique panels:", length(unique(combined$sector_label)), "\n")
print(
  combined %>%
    summarise(
      min_est = min(estimate, na.rm = TRUE),
      max_est = max(estimate, na.rm = TRUE),
      min_time = min(event_time, na.rm = TRUE),
      max_time = max(event_time, na.rm = TRUE)
    )
)

# ---- ATT Comparison at t+3 ----------------------------------------------------
cat("\nATT Comparison at t+3 (TWFE vs dCdH, with CI and SE)\n")

att_summary <- combined %>%
  filter(event_time == 3) %>%
  group_by(funder, ad_sector_names, method) %>%
  summarise(
    estimate = mean(estimate, na.rm = TRUE),
    se = mean(ci_upper - ci_lower, na.rm = TRUE) / (2 * 1.96),
    ci_lower = mean(ci_lower, na.rm = TRUE),
    ci_upper = mean(ci_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = c(estimate, se, ci_lower, ci_upper),
    names_glue = "{.value}_{method}"
  ) %>%
  # Drop cases with NA in either estimator (e.g., missing dCdH)
  filter(!is.na(estimate_TWFE) & !is.na(estimate_dCdH)) %>%
  mutate(
    diff = estimate_TWFE - estimate_dCdH
  ) %>%
  arrange(desc(abs(diff)))

# Preview in console
print(att_summary, n = 24)

# Save full detailed table with CI and SE
write_csv(att_summary, file.path(table_path, "Table_ATT_Comparison_t3_Full.csv"))
cat("\n[OK] Saved detailed ATT comparison table (with SE & CI) to:", table_path, "\n")

# ---- Keep only complete panels for plotting -----------------------------------
valid_panels <- att_summary$ad_sector_names
combined <- combined %>%
  filter(ad_sector_names %in% valid_panels)

# ---- Reorder Panels by Magnitude (dCdH avg) -----------------------------------
combined <- combined %>%
  group_by(sector_label) %>%
  mutate(mean_effect = mean(estimate[method == "dCdH"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sector_label = fct_reorder(sector_label, mean_effect, .fun = max))

# ---- Faceted Plot: All Panels -------------------------------------------------
cat("\nCreating faceted event-study comparison plot...\n")

fig_faceted <- ggplot(combined, aes(x = event_time, y = estimate, color = method, shape = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = method),
    alpha = 0.15, colour = NA
  ) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2, position = position_dodge(width = 0.3), alpha = 0.6
  ) +
  geom_point(size = 1.8, position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c("TWFE" = "#440154", "dCdH" = "#E66100")) + # dCdH orange
  scale_fill_manual(values = c("TWFE" = "#440154", "dCdH" = "#E66100")) +
  facet_wrap(~sector_label, ncol = 4, scales = "free_y") +
  labs(
    title = "Event-Study Estimates of Aid Impacts on IWI — TWFE vs dCdH",
    subtitle = "Each panel shows a donor–sector combination",
    x = "Event time",
    y = "Estimated Effect on IWI",
    color = "Estimator",
    shape = "Estimator",
    fill = "Estimator"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 8.5, face = "bold"),
    axis.text.x = element_text(size = 8)
  )

ggsave(file.path(plot_path, "Fig_TWFE_dCdH_Faceted_AllPanels.png"),
  fig_faceted,
  width = 10, height = 9, dpi = 300
)
cat("\n[OK] Saved full faceted figure to:", plot_path, "\n")

# ---- Focused Plot: Six Key Panels ---------------------------------------------
selected_panels <- c(
  "Health — China",
  "Education — World Bank",
  "Transport and Storage — China",
  "Water Supply and Sanitation — World Bank",
  "Communications — World Bank",
  "Government and Civil Society — China"
)

fig_focus <- combined %>%
  filter(sector_label %in% selected_panels) %>%
  ggplot(aes(x = event_time, y = estimate, color = method, shape = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = method),
    alpha = 0.15, colour = NA
  ) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2, position = position_dodge(0.3), alpha = 0.6
  ) +
  geom_point(size = 1.8, position = position_dodge(0.3)) +
  scale_color_manual(values = c("TWFE" = "#440154", "dCdH" = "#E66100")) +
  scale_fill_manual(values = c("TWFE" = "#440154", "dCdH" = "#E66100")) +
  facet_wrap(~sector_label, ncol = 3, scales = "free_y") +
  labs(
    title = "Event-Study Estimates of Aid Impacts on IWI — Selected Panels",
    subtitle = "Each panel shows a donor–sector combination",
    x = "Event time (years relative to project start)",
    y = "Estimated Effect on IWI",
    color = "Estimator",
    shape = "Estimator",
    fill = "Estimator"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9)
  )

ggsave(file.path(plot_path, "Fig_TWFE_dCdH_SelectedPanels.png"),
  fig_focus,
  width = 9, height = 6.5, dpi = 300
)
cat("\n[OK] Saved focused six-panel figure to:", plot_path, "\n")
################################################################################

################################################################################
# Parallel Trends Diagnostics — Pre-treatment (t = -1)
# Mattias Antar, 2025
################################################################################

library(tidyverse)
library(readr)

# ---- Paths --------------------------------------------------------------------
table_path <- here::here("output", "tables")
plot_path <- here::here("output", "plots")
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

cat("\nChecking parallel trends assumption (pre-treatment, event_time = -1)...\n")

# ---- Extract pre-treatment coefficients ---------------------------------------
pretrend <- combined %>%
  filter(event_time == -1) %>%
  group_by(funder, ad_sector_names, method) %>%
  summarise(
    estimate = mean(estimate, na.rm = TRUE),
    ci_lower = mean(ci_lower, na.rm = TRUE),
    ci_upper = mean(ci_upper, na.rm = TRUE),
    sig = ifelse(ci_lower > 0 | ci_upper < 0, TRUE, FALSE),
    .groups = "drop"
  )

# ---- Summary statistics by estimator ------------------------------------------
pretrend_summary <- pretrend %>%
  group_by(method) %>%
  summarise(
    mean_pre = mean(estimate, na.rm = TRUE),
    sd_pre = sd(estimate, na.rm = TRUE),
    share_sig = mean(sig, na.rm = TRUE),
    n_panels = n()
  )

cat("\nQC — Summary of pre-treatment coefficients:\n")
print(pretrend_summary)

# ---- Simple pooled t-tests ----------------------------------------------------
t_twfe <- t.test(pretrend$estimate[pretrend$method == "TWFE"])
t_dcdh <- t.test(pretrend$estimate[pretrend$method == "dCdH"])

cat("\n--- TWFE pre-treatment mean test ---\n")
print(t_twfe)
cat("\n--- dCdH pre-treatment mean test ---\n")
print(t_dcdh)

# ---- Save detailed results ----------------------------------------------------
write_csv(pretrend, file.path(table_path, "Table_Pretrend_Check.csv"))
write_csv(pretrend_summary, file.path(table_path, "Table_Pretrend_Summary.csv"))

# Save t-test summaries to text file
sink(file.path(table_path, "Table_Pretrend_TTests.txt"))
cat("Parallel Trends Diagnostic — Pre-treatment (t = -1)\n")
cat("====================================================\n\n")
cat("Summary statistics by estimator:\n")
print(pretrend_summary)
cat("\n\nTWFE mean test:\n")
print(t_twfe)
cat("\n\ndCdH mean test:\n")
print(t_dcdh)
sink()
cat("\n[OK] Saved pretrend summary tables and tests to:", table_path, "\n")

# ---- Visualization 1: Boxplot summary ----------------------------------------
fig_pretrend_box <- ggplot(pretrend, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = 16, outlier.alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Distribution of Pre-treatment (t = -1) Coefficients",
    subtitle = "Parallel Trends Diagnostic — TWFE vs dCdH",
    y = "Estimated Pre-treatment Effect (IWI points)",
    x = "Estimator"
  ) +
  scale_fill_manual(values = c("TWFE" = "#440154", "dCdH" = "#E66100")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave(file.path(plot_path, "Fig_Pretrend_Boxplot.png"),
  fig_pretrend_box,
  width = 6, height = 5, dpi = 300
)
cat("\n[OK] Saved boxplot of pre-treatment coefficients to:", plot_path, "\n")

# ---- Visualization 2: Dot-whisker by panel ------------------------------------
fig_pretrend_dot <- ggplot(pretrend, aes(
  x = reorder(ad_sector_names, estimate),
  y = estimate, color = method
)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.3, position = position_dodge(width = 0.7)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  coord_flip() +
  facet_wrap(~funder, scales = "free_y") +
  scale_color_manual(values = c("TWFE" = "#440154", "dCdH" = "#E66100")) +
  labs(
    title = "Pre-treatment (t = -1) Coefficients by Sector and Funder",
    subtitle = "Parallel Trends Diagnostic — Confidence Intervals shown",
    x = "Sector",
    y = "Pre-treatment Estimate (IWI points)",
    color = "Estimator"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 9, face = "bold")
  )

ggsave(file.path(plot_path, "Fig_Pretrend_DotWhisker.png"),
  fig_pretrend_dot,
  width = 8.5, height = 7, dpi = 300
)
cat("\n[OK] Saved detailed dot-whisker pretrend figure to:", plot_path, "\n")

################################################################################
# Mean Wealth (IWI) by Funder and Sector: 1990–1992 vs 2017–2019 (Faceted Dumbbell)
################################################################################

library(tidyverse)
library(janitor)
library(readr)

# ---- Paths --------------------------------------------------------------------
data_path <- here::here("Data", "Archive_enriched")
plot_path <- here::here("output", "plots")
table_path <- here::here("output", "tables")
sector_map_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Starting corrected mean IWI by funder–sector (both funders) ===\n")

# ---- Load sector map -----------------------------------------------------------
sector_map <- read_csv(sector_map_path, show_col_types = FALSE) %>% clean_names()

# ---- Helper: Extract metadata from file names ----------------------------------
extract_meta <- function(filepath) {
  filename <- basename(filepath)
  funder <- case_when(
    str_detect(filename, "_wb_") ~ "World Bank",
    str_detect(filename, "_ch_") ~ "China",
    TRUE ~ NA_character_
  )
  sector_code <- as.integer(str_extract(filename, "(?<=_)\\d{3}(?=_DiD)"))
  tibble(filepath = filepath, funder = funder, ad_sector_codes = sector_code)
}

# ---- Identify files ------------------------------------------------------------
files <- list.files(data_path, pattern = "InputData_.*_DiD_enriched\\.csv$", full.names = TRUE)
cat("Files found:", length(files), "\n")

meta <- map_dfr(files, extract_meta)
cat("\nExample of extracted metadata:\n")
print(head(meta, 6))

# ---- Read and label each dataset properly --------------------------------------
read_labeled <- function(file, funder, sector_code) {
  df <- suppressMessages(read_csv(file, show_col_types = FALSE)) %>%
    clean_names() %>%
    mutate(funder = funder, ad_sector_codes = sector_code)
}

raw_all <- map2_dfr(meta$filepath, seq_len(nrow(meta)), \(file, i) {
  read_labeled(file, meta$funder[i], meta$ad_sector_codes[i])
})

cat("\n[OK] Loaded all enriched data files with correct funder–sector mapping.\n")
cat("Rows:", nrow(raw_all), "\n")
cat("Distinct funders:", unique(raw_all$funder), "\n")
cat("Distinct sector codes:", n_distinct(raw_all$ad_sector_codes), "\n")

# ---- Join sector names ---------------------------------------------------------
raw_all <- raw_all %>%
  left_join(sector_map, by = "ad_sector_codes") %>%
  filter(!(ad_sector_names == "Trade and Tourism" & funder == "World Bank"))

cat("\n[OK] Added sector names.\n")
cat("Unique panels after filtering:", n_distinct(paste(raw_all$funder, raw_all$ad_sector_names)), "\n")

# ---- Compute mean IWI per period ----------------------------------------------
iwi_change <- raw_all %>%
  group_by(funder, ad_sector_names) %>%
  summarise(
    mean_iwi_1990_1992 = mean(iwi_1990_1992, na.rm = TRUE),
    mean_iwi_2017_2019 = mean(iwi_2017_2019, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(delta_iwi = mean_iwi_2017_2019 - mean_iwi_1990_1992)

cat("\n--- Preview of mean IWI by funder–sector ---\n")
print(iwi_change %>% arrange(funder, ad_sector_names) %>% head(10))
cat("\nFunders in table:", unique(iwi_change$funder), "\n")

# ---- Reshape for plotting ------------------------------------------------------
iwi_long <- iwi_change %>%
  pivot_longer(
    cols = starts_with("mean_iwi_"),
    names_to = "period",
    values_to = "mean_iwi"
  ) %>%
  mutate(
    period = str_replace_all(period, "mean_iwi_", ""),
    period = str_replace_all(period, "_", ":"),
    period = factor(period, levels = c("1990:1992", "2017:2019"))
  )

cat("\n[OK] Prepared long-form data for plotting.\n")
print(iwi_long %>% count(funder, period))

# ---- Faceted Dumbbell Plot with Arrows ----------------------------------------
fig_b <- ggplot(iwi_change, aes(y = reorder(ad_sector_names, mean_iwi_2017_2019))) +
  geom_segment(
    aes(x = mean_iwi_1990_1992, xend = mean_iwi_2017_2019, color = funder),
    linewidth = 1.1, alpha = 0.8,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed") # <-- Added arrowheads
  ) +
  geom_point(aes(x = mean_iwi_1990_1992, color = funder), size = 3, alpha = 0.6) +
  geom_point(aes(x = mean_iwi_2017_2019, color = funder), size = 3.4) +
  facet_wrap(~funder, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("China" = "#E66100", "World Bank" = "#2b83ba")) +
  labs(
    title = "Mean Wealth (IWI) by Funder-Sector",
    subtitle = "1990–1992 (left) → 2017–2019 (right)",
    x = "Mean IWI",
    y = "Sector",
    caption = "Arrows indicate direction of change; mean IWI increased across all sectors during the period."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 9)
  )

print(fig_b)

# Save figure
ggsave(file.path(plot_path, "Fig_IWI_MeanChange_FunderSector_1990_2019_FACETED.png"),
  fig_b,
  width = 9, height = 7, dpi = 300
)
###############################################################################

## Faceted dot plot showing mean wealth (IWI) by sector and donor , wb vs ch 1990 vs 2017

###############################################################################
library(tidyverse)

# --- prepare data --------------------------------------------------------------
iwi_long_mean <- iwi_change %>%
  mutate(
    funder = recode(funder,
      "China" = "China",
      "World Bank" = "World Bank"
    ),
    ad_sector_names = str_trunc(ad_sector_names, 35)
  ) %>%
  pivot_longer(
    cols = starts_with("mean_iwi_"),
    names_to = "period",
    values_to = "mean_iwi"
  ) %>%
  mutate(
    period = recode(
      period,
      "mean_iwi_1990_1992" = "1990–1992",
      "mean_iwi_2017_2019" = "2017–2019"
    ),
    period = factor(period, levels = c("1990–1992", "2017–2019"))
  )

# --- plot (points only, black and white) ---------------------------------------
ggplot(
  iwi_long_mean,
  aes(
    x = mean_iwi,
    y = reorder(ad_sector_names, mean_iwi)
  )
) +
  geom_point(aes(shape = period),
    color = "black",
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~funder, scales = "free_y", ncol = 1) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  labs(
    title = "Mean IWI by Funder–Sector",
    subtitle = "1990–1992 vs 2017–2019",
    x = "Mean IWI",
    y = "Sector"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 9)
  )


############## TWFE TWENTY THREE PANELS INTO ONE PLOT #################
################################################################################
# VISUALIZATION: Create the Faceted Plot with Independent Y-Axes
################################################################################
# This code assumes the 'plot_data' dataframe from the previous step is loaded.

# --- Define save path for the main faceted plot ---
plots_save_path <- here::here("output", "plots")
dir.create(plots_save_path, showWarnings = FALSE)

# --- Create the plot object with the crucial 'scales = "free_y"' addition ---
faceted_event_study_plot_v2 <- ggplot(plot_data, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 2) +

  # KEY CHANGE HERE: `scales = "free_y"` lets each plot have its own y-axis range
  facet_wrap(~ funder_full + ad_sector_names, ncol = 4, scales = "free_y") +
  labs(
    title = "TWFE Event-Study Models by Funder and Sector",
    x = "Event Time",
    y = "Change in IWI"
  ) +
  scale_x_continuous(breaks = -1:3) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

# --- Display the improved plot for inspection ---
print(faceted_event_study_plot_v2)

#--- Save the improved plot ---
ggsave(
  filename = file.path(plots_save_path, "TWFE_Fig_Combined_Event_Study.png"),
  plot = faceted_event_study_plot_v2,
  width = 14,
  height = 10,
  dpi = 300
)
cat("Successfully saved the main faceted plot to:", plots_save_path, "\n")

############# dCdH Model 23 panels into one plot ##############
################################################################################
# 1. SETUP: Load Libraries and Define Paths
################################################################################

# Install pacman if you haven't already, it helps manage packages
if (!require("pacman")) install.packages("pacman")

# Load required libraries
pacman::p_load(tidyverse, stringr)

# --- Define paths to your data ---
dcdh_results_path <- here::here("Stata", "dcdh_results")
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
plots_save_path <- here::here("output", "plots") # General plot save folder

# Ensure the save directory exists
dir.create(plots_save_path, showWarnings = FALSE)

################################################################################
# 2. DATA PREPARATION: Load, Parse, and Clean dCdH Results (CORRECTED)
################################################################################

# --- Load and combine all dCdH result files into a raw data frame ---
dcdh_raw_data <- list.files(path = dcdh_results_path, pattern = "^DiD_Table_.*\\.csv$", full.names = TRUE) %>%
  map_dfr(~ {
    # ---- A. Parse filename for metadata ----
    filename <- basename(.x)
    matches <- str_match(filename, "DiD_Table_([a-z]+)_([0-9]+)\\.csv")
    funder_abbr <- matches[, 2]
    sector_code <- as.numeric(matches[, 3])

    # ---- B. Read and clean the messy CSV ----
    raw_data <- read_csv(.x, col_names = FALSE, show_col_types = FALSE) %>%
      mutate(
        term = str_remove_all(X1, '["=]'),
        estimate_raw = str_remove_all(X2, '["=]'),
        # --- THIS IS THE CORRECTED LINE ---
        std_error_raw = str_remove_all(X3, '[()"=]') # Added '=' to the characters to remove
      ) %>%
      select(term, estimate_raw, std_error_raw) %>%
      filter(str_starts(term, "Effect_|Placebo_"))

    # ---- C. Reshape and calculate values ----
    raw_data %>%
      mutate(
        estimate = as.numeric(str_extract(estimate_raw, "-?[0-9.]+")),
        std_error = as.numeric(std_error_raw),
        event_time = case_when(
          term == "Placebo_1" ~ -1,
          term == "Effect_1" ~ 1,
          term == "Effect_2" ~ 2,
          term == "Effect_3" ~ 3,
          TRUE ~ NA_integer_
        ),
        conf.low = estimate - 1.96 * std_error,
        conf.high = estimate + 1.96 * std_error
      ) %>%
      mutate(funder_abbr = funder_abbr, sector_code = sector_code)
  })

# --- Create a separate tibble for the t=0 anchor points ---
anchor_points <- dcdh_raw_data %>%
  distinct(funder_abbr, sector_code) %>%
  mutate(event_time = 0, estimate = 0, conf.low = 0, conf.high = 0)

# --- Combine the original data with the anchor points ---
dcdh_plot_data <- bind_rows(dcdh_raw_data, anchor_points) %>%
  arrange(funder_abbr, sector_code, event_time)


# --- Load lookup tables and merge with plot data ---
sector_lookup <- read_csv(sector_lookup_path, show_col_types = FALSE) %>%
  select(ad_sector_codes, ad_sector_names)

funder_lookup <- tibble(
  funder_abbr = c("ch", "wb"),
  funder_full = c("China", "World Bank")
)

# Merge and finalize the dataset
final_dcdh_data <- dcdh_plot_data %>%
  filter(!(funder_abbr == "wb" & sector_code == 330)) %>%
  left_join(funder_lookup, by = "funder_abbr") %>%
  left_join(sector_lookup, by = c("sector_code" = "ad_sector_codes")) %>%
  filter(!is.na(event_time))

# Quick check to confirm the fix
glimpse(final_dcdh_data)


################################################################################
# 3. VISUALIZATION: Create the Faceted dCdH Plot
################################################################################

dcdh_faceted_plot <- ggplot(final_dcdh_data, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#d95f02", alpha = 0.25) +
  geom_line(color = "#d95f02", linewidth = 0.8) +
  geom_point(color = "#d95f02", size = 2) +
  facet_wrap(~ funder_full + ad_sector_names, ncol = 4, scales = "free_y") +
  labs(
    title = "dCdH Event-Study Models by Funder and Sector",
    x = "Event Time",
    y = "Change in IWI"
  ) +
  scale_x_continuous(breaks = -1:3) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

# --- Display the plot for inspection ---
print(dcdh_faceted_plot)

# --- Save the plot ---
ggsave(
  filename = file.path(plots_save_path, "dCdH_Fig_Combined_Event_Study.png"),
  plot = dcdh_faceted_plot,
  width = 14,
  height = 10,
  dpi = 300
)

cat("Successfully saved the main dCdH faceted plot to:", plots_save_path, "\n")

########## MAP PLOT FOR SPATIAL SPILLOVER #############
################################################################################
# SETUP FOR MAP: Load Packages and Prep Data
################################################################################
# (Ensuring 'analytic' and 'africa_map' are loaded from the previous step)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, sf, rnaturalearth)

# --- Define save path ---
save_path <- here::here("output", "spillover")
# Ensure the directory exists (it should, but good practice)
dir.create(save_path, showWarnings = FALSE)

# --- Re-run the main model to get residuals ---
files <- list.files(here::here("Data", "Archive_enriched"), pattern = "enriched\\.csv$", full.names = TRUE)
raw <- map_dfr(files, ~ read_csv(.x, show_col_types = FALSE) %>% janitor::clean_names())
unit_level <- raw %>%
  group_by(iso3, adm2) %>%
  summarise(across(c(lat, lon, iwi_2014_2016, iwi_2017_2019, log_avg_min_to_city, log_avg_pop_dens, log_3yr_pre_conflict_deaths), ~ mean(.x, na.rm = TRUE)),
    treated = as.integer(any(treated == 1 & time <= 2012)), .groups = "drop"
  ) %>%
  mutate(iwi_est_post_oda = coalesce(rowMeans(cbind(iwi_2014_2016, iwi_2017_2019), na.rm = TRUE), iwi_2017_2019, iwi_2014_2016), country = iso3)
analytic <- unit_level %>%
  select(iwi_est_post_oda, treated, log_avg_min_to_city, log_avg_pop_dens, log_3yr_pre_conflict_deaths, country, lat, lon) %>%
  drop_na()

m_main <- lm(iwi_est_post_oda ~ treated + log_avg_min_to_city + log_avg_pop_dens + log_3yr_pre_conflict_deaths + country, data = analytic)
analytic$residuals <- residuals(m_main)

# --- Get the map data for Africa ---
africa_map <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")


################################################################################
# MAP VISUAL 1: "Bubble Map" of Residuals
################################################################################

# --- Create the plot ---
plot_map_bubbles_revised <- ggplot() +
  # 1. Draw the map of Africa as the base layer
  geom_sf(data = africa_map, fill = "gray90", color = "gray70", linewidth = 0.3) +

  # 2. Add the residual points on top
  geom_point(
    data = analytic,
    aes(x = lon, y = lat, color = residuals),
    alpha = 0.7,
    size = 1.5
  ) +

  # 3. Add the color scale
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +

  # 4. Set map limits
  coord_sf(xlim = c(-20, 55), ylim = c(-35, 38)) +

  # 5. Set revised labs and theme
  labs(
    title = "Spatial Distribution of Model Residuals",
    caption = "Each point is a DHS cluster", # Moved subtitle to caption
    x = NULL, # Removes x-axis label
    y = NULL, # Removes y-axis label
    color = "Residual"
  ) +
  theme_void(base_size = 14) + # Use theme_void to remove grid/axes
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center the title
    plot.caption = element_text(hjust = 0.5, size = 10), # Center the caption
    legend.position = "bottom" # Move legend to the bottom
  )

# --- Display the plot ---
print(plot_map_bubbles_revised)

# --- Save the plot ---
ggsave(
  filename = file.path(save_path, "map_plot_residuals.png"),
  plot = plot_map_bubbles_revised,
  width = 8, # Specify width in inches
  height = 7, # Specify height in inches
  dpi = 300 # Set resolution
)

cat("Successfully saved map plot to:", file.path(save_path, "map_plot_residuals.png"), "\n")

################################################################################
# TWFE ROBUSTNESS COMPARISON: Final Thesis Version
################################################################################

# 1. Data Cleaning and Label Standardization
twfe_comp_clean <- twfe_comp_final %>%
  # Filter out the specific panel
  filter(panel_name != "World Bank: Trade and Tourism") %>%
  # Rename models to match your desired legend text exactly
  mutate(model = case_when(
    model == "Baseline (No Covars)" ~ "Baseline (No Covars)",
    str_detect(model, "Adjusted") ~ "Adjusted (6 Covars)",
    TRUE ~ model
  ))

# 2. Generate the Final Plot
p_twfe_robustness_final <- ggplot(twfe_comp_clean, aes(x = time, y = estimate, color = model, group = model)) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +

  # Confidence Intervals (Ribbons)
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high, fill = model),
    alpha = 0.15, color = NA
  ) +

  # Lines and Points
  geom_line(linewidth = 1) +
  geom_point(size = 2) +

  # Faceting: 4 columns for balanced thesis layout
  facet_wrap(~panel_name, scales = "free_y", ncol = 4) +

  # Thesis Standard Colors: Light Blue (Baseline) vs. Light Red (Adjusted)
  # Note: Ensure names here match the mutated names in Step 1
  scale_color_manual(values = c(
    "Baseline (No Covars)" = "#3498DB",
    "Adjusted (6 Covars)" = "#E74C3C"
  )) +
  scale_fill_manual(values = c(
    "Baseline (No Covars)" = "#3498DB",
    "Adjusted (6 Covars)" = "#E74C3C"
  )) +

  # Theme and Styling
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 7, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 10)
  ) +
  labs(
    title = "TWFE Models Comparison",
    subtitle = "Baseline (No Covariates) and with Covariates",
    x = "Event Time",
    y = "Estimated Effect (IWI)",
    color = "Model Specification:",
    fill = "Model Specification:"
  )

# 3. Save the Plot
# Creating directory if it doesn't exist just in case
if (!dir.exists(here::here("output", "plots"))) dir.create(here::here("output", "plots"), recursive = TRUE)

ggsave(
  filename = here::here("output", "plots", "TWFE_Comparison_Final_Thesis.png"),
  plot = p_twfe_robustness_final,
  width = 16,
  height = 12,
  dpi = 300
)

# Show the plot
print(p_twfe_robustness_final)

cat("\n[OK] SUCCESS: Final TWFE Comparison plot saved to output/plots/\n")

################################################################################
# dCdH ROBUSTNESS COMPARISON: Final Thesis Version
################################################################################

pacman::p_load(tidyverse, janitor, stringr)

# --- 1. LOAD AND PREP ADJUSTED DATA (WIDE -> LONG) ---
dcdh_adj_path <- here::here("output", "dCdH", "dCdH_Full_Stats_Wide.csv")
sector_lookup <- read_csv(here::here("Data", "Archive_enriched", "sector_group_names.csv"), show_col_types = FALSE) %>%
  clean_names()

dcdh_adj_long <- read_csv(dcdh_adj_path, show_col_types = FALSE) %>%
  clean_names() %>%
  # Pivot the estimates and CIs
  pivot_longer(
    cols = starts_with(c("estimate", "conf_low", "conf_high")),
    names_to = c(".value", "time_raw"),
    names_pattern = "(.*)_(t.*)"
  ) %>%
  # Map time strings to integers
  mutate(time = case_when(
    time_raw == "t_minus_1" ~ -1,
    time_raw == "t0" ~ 0,
    time_raw == "t_plus_1" ~ 1,
    time_raw == "t_plus_2" ~ 2,
    time_raw == "t_plus_3" ~ 3
  )) %>%
  mutate(model = "Adjusted (6 Covars)") %>%
  select(funder, ad_sector_names, time, estimate, conf_low, conf_high, model)

# --- 2. LOAD AND PREP BASELINE DATA (PARSING STRINGS) ---
baseline_dir <- here::here("Stata", "dcdh_results")
baseline_files <- list.files(baseline_dir, pattern = "\\.csv$", full.names = TRUE)

dcdh_base_list <- list()

for (f in baseline_files) {
  # Parse funder and sector code from filename (e.g., DiD_Table_ch_110.csv)
  fname <- basename(f)
  m <- str_match(fname, "Table_([a-z]+)_([0-9]+)")
  funder_abbr <- m[, 2]
  sector_code <- as.numeric(m[, 3])

  # Read and clean headers
  df <- read_csv(f, show_col_types = FALSE, skip = 1, col_names = c("label", "est_str", "se_str")) %>%
    filter(str_detect(label, "Effect|Placebo")) %>%
    mutate(
      # Extract numbers from strings like ="0.440***" or ="(0.119)"
      estimate = as.numeric(str_extract(est_str, "-?[0-9]+\\.[0-9]+")),
      se = as.numeric(str_extract(se_str, "[0-9]+\\.[0-9]+")),
      # Calculate CIs (95%)
      conf_low = estimate - (1.96 * se),
      conf_high = estimate + (1.96 * se),
      # Map labels to time
      time = case_when(
        str_detect(label, "Placebo_1") ~ -1,
        str_detect(label, "Effect_1") ~ 1,
        str_detect(label, "Effect_2") ~ 2,
        str_detect(label, "Effect_3") ~ 3
      ),
      model = "Baseline (No Covars)",
      funder = ifelse(funder_abbr == "ch", "China", "World Bank"),
      ad_sector_codes = sector_code
    ) %>%
    # Add reference period 0
    add_row(
      time = 0, estimate = 0, conf_low = 0, conf_high = 0,
      model = "Baseline (No Covars)", funder = ifelse(funder_abbr == "ch", "China", "World Bank"),
      ad_sector_codes = sector_code
    )

  dcdh_base_list[[fname]] <- df
}

# Merge all baseline and join names
dcdh_base_long <- bind_rows(dcdh_base_list) %>%
  left_join(sector_lookup, by = "ad_sector_codes") %>%
  select(funder, ad_sector_names, time, estimate, conf_low, conf_high, model)

# --- 3. MERGE AND CLEAN FINAL DATASET ---
dcdh_comp_final <- bind_rows(dcdh_adj_long, dcdh_base_long) %>%
  # Standardize panel names
  mutate(panel_name = paste0(funder, ": ", ad_sector_names)) %>%
  # EXCLUDE the requested panel
  filter(panel_name != "World Bank: Trade and Tourism") %>%
  # Filter for the common time window
  filter(!is.na(time) & !is.na(estimate))

# --- 4. GENERATE PLOT ---
p_dcdh_comparison_final <- ggplot(dcdh_comp_final, aes(x = time, y = estimate, color = model, group = model)) +
  # Reference Lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +

  # Error Ribbons
  geom_ribbon(aes(ymin = conf_low, ymax = conf_high, fill = model),
    alpha = 0.15, color = NA
  ) +

  # Lines and Points
  geom_line(linewidth = 1) +
  geom_point(size = 2) +

  # Faceting
  facet_wrap(~panel_name, scales = "free_y", ncol = 4) +

  # Colors: Blue (Old) vs Red (New)
  scale_color_manual(values = c(
    "Baseline (No Covars)" = "#3498DB",
    "Adjusted (6 Covars)" = "#E74C3C"
  )) +
  scale_fill_manual(values = c(
    "Baseline (No Covars)" = "#3498DB",
    "Adjusted (6 Covars)" = "#E74C3C"
  )) +

  # Theme and Labels
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 7, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 10)
  ) +
  labs(
    title = "dCdH Models Comparison",
    subtitle = "Baseline (No Covariates) and with Covariates (Residualized)",
    x = "Event Time",
    y = "Estimated Effect (IWI)",
    color = "Model Specification:",
    fill = "Model Specification:"
  )

# --- 5. SAVE ---
# Create directory if it doesn't exist
if (!dir.exists(here::here("output", "plots"))) dir.create(here::here("output", "plots"), recursive = TRUE)

ggsave(
  filename = here::here("output", "plots", "dCdH_Comparison_Final_Thesis.png"),
  plot = p_dcdh_comparison_final,
  width = 16,
  height = 12,
  dpi = 300
)

# Show the plot
print(p_dcdh_comparison_final)

cat("\n[OK] SUCCESS: Final dCdH Comparison plot saved to output/plots/\n")
