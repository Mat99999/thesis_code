################################################################################
# TWFE Event-Study Code
################################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(fixest, tidyverse, stringr, here)

# --- 1) DATA PATH (using here::here for portability) ---
data_path <- here::here("Data", "Archive_enriched")

# Optional sector lookup
sector_lookup_path <- here::here("Data", "Archive_enriched", "sector_group_names.csv")
sector_lookup <- tryCatch(
  readr::read_csv(sector_lookup_path, show_col_types = FALSE),
  error = function(e) tibble(ad_sector_codes = NA_real_, ad_sector_names = NA_character_)
)

funder_lookup <- tibble(
  funder_abbr = c("ch", "wb"),
  funder_full = c("China", "World Bank")
)

# --- 2) FILES ---
data_files <- list.files(path = data_path, pattern = "^InputData_.*\\.csv$", full.names = TRUE)
cat("Found", length(data_files), "files in Archive.\n")

# --- F) SETUP SAVE PATHS (using here::here for portability) ---
base_save_path <- here::here("output", "TWFE", "without_covariates")
tables_path <- here::here("output", "TWFE", "without_covariates", "tables")
plots_path <- here::here("output", "TWFE", "without_covariates", "plots")

# Ensure directories exist (Quality Control)
dir.create(tables_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)
cat(paste0("Saving tables to: ", tables_path, "\n"))
cat(paste0("Saving plots to: ", plots_path, "\n"))

# Initialize plot_data to accumulate results for parallel trends analysis
plot_data <- tibble()

# --- 3) MAIN LOOP ---
for (file in data_files) {
  # A) Strong regex parse: InputData_<funder>_<sector>*.csv
  filename <- basename(file)
  m <- str_match(filename, "^InputData_([A-Za-z]+)_([0-9]+).*\\.csv$")
  funder_abbr <- tolower(m[, 2])
  sector_code <- suppressWarnings(as.numeric(m[, 3]))

  # Funder name (1 value, never NA)
  funder_full <- funder_lookup %>%
    filter(funder_abbr == !!funder_abbr) %>%
    distinct(funder_full) %>%
    slice_head(n = 1) %>%
    pull(funder_full)
  if (length(funder_full) == 0 || is.na(funder_full)) {
    funder_full <- str_to_title(ifelse(is.na(funder_abbr), "Unknown", funder_abbr))
  }

  # Sector name (1 value, never NA)
  sector_full <- sector_lookup %>%
    filter(ad_sector_codes == !!sector_code) %>%
    distinct(ad_sector_names) %>%
    slice_head(n = 1) %>%
    pull(ad_sector_names)
  if (length(sector_full) == 0 || is.na(sector_full)) {
    sector_full <- paste0(sector_code) # fallback to code
  }

  cat("\n*****************************************************\n")
  cat(paste0("*** Processing: ", funder_full, " – ", sector_full, " ***\n"))
  cat("*****************************************************\n")

  # B) Load + prepare
  df <- suppressMessages(readr::read_csv(file, show_col_types = FALSE))
  needed <- c("dhs_id", "time", "treated", "iwi_est_post_oda")
  if (!all(needed %in% names(df))) {
    cat("[WARNING] Missing required columns - skipping.\n")
    next
  }

  df <- df %>%
    mutate(
      dhs_id = as.integer(dhs_id),
      time = as.integer(time),
      treated = as.numeric(treated),
      iwi_est_post_oda = as.numeric(iwi_est_post_oda)
    ) %>%
    group_by(dhs_id) %>%
    mutate(
      treatment_cohort = if_else(any(treated == 1, na.rm = TRUE),
        min(time[treated == 1], na.rm = TRUE),
        Inf
      )
    ) %>%
    ungroup() %>%
    mutate(relative_time = if_else(is.finite(treatment_cohort), time - treatment_cohort, NA_real_)) %>%
    # bin endpoints to match your original robustness choice
    mutate(relative_time = case_when(
      relative_time <= -5 ~ -5,
      relative_time >= 5 ~ 5,
      TRUE ~ relative_time
    ))

  # Reference period: last pre
  ref_period <- df %>%
    filter(relative_time < 0) %>%
    summarise(ref = suppressWarnings(max(relative_time, na.rm = TRUE))) %>%
    pull(ref)
  if (is.na(ref_period) || !is.finite(ref_period)) {
    cat("[WARNING] No pre-treatment - skipping.\n")
    next
  }
  cat("Reference period t =", ref_period, "\n")

  # C) TWFE
  twfe_model <- tryCatch(
    {
      feols(iwi_est_post_oda ~ i(relative_time, ref = ref_period) | dhs_id + time,
        data = df, cluster = ~dhs_id
      )
    },
    error = function(e) {
      cat("[ERROR] Model error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(twfe_model)) next

  # D) Tidy + build 5 event-time points
  est <- broom::tidy(twfe_model, conf.int = TRUE) %>%
    mutate(relative_time = suppressWarnings(as.numeric(str_extract(term, "-?[0-9]+")))) %>%
    filter(!is.na(relative_time)) %>%
    # include the reference (effect = 0)
    add_row(
      term = "Reference", estimate = 0, conf.low = 0, conf.high = 0,
      relative_time = ref_period, .before = 1
    ) %>%
    arrange(relative_time)

  # Available RT grid
  rt_vals <- sort(unique(est$relative_time))
  pre_vals <- rt_vals[rt_vals < ref_period]
  post_vals <- rt_vals[rt_vals > ref_period]

  # Map to exactly five event-time points (-1,0,1,2,3)
  # -1 := closest pre period to ref (handles -5 bin)
  #  0 := ref
  # +1,+2 := first/second post
  # +3 := third post if present, else the far-right bin (handles +5 bin)
  et_map <- tibble(
    event_time = c(-1, 0, 1, 2, 3),
    relative_time = c(
      if (length(pre_vals) >= 1) max(pre_vals) else NA_real_,
      ref_period,
      if (length(post_vals) >= 1) post_vals[1] else NA_real_,
      if (length(post_vals) >= 2) post_vals[2] else if (length(post_vals) >= 1) max(post_vals) else NA_real_,
      if (length(post_vals) >= 3) post_vals[3] else if (length(post_vals) >= 1) max(post_vals) else NA_real_
    )
  )

  results_for_plot <- et_map %>%
    left_join(est %>% select(relative_time, estimate, conf.low, conf.high),
      by = "relative_time"
    )

  # E) Plot (dark blue line + darker CI; x: -1..3)
  p <- ggplot(results_for_plot, aes(x = event_time, y = estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
      fill = "#08306B", alpha = 0.35, na.rm = TRUE
    ) +
    geom_line(color = "#084594", linewidth = 1.1, na.rm = TRUE) +
    geom_point(color = "#084594", size = 2.8, na.rm = TRUE) +
    geom_hline(yintercept = 0, color = "gray40", linetype = "dashed") +
    scale_x_continuous(breaks = -1:3, limits = c(-1, 3)) +
    labs(
      title = paste0("Effect of Aid: ", funder_full, ", ", sector_full, " Sector"),
      x = "Event Time",
      y = "Effect on IWI"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  print(p)

  # F) SAVE RESULTS (UPDATED CODE)

  # Create a clean filename slug
  file_slug <- paste(funder_abbr, sector_code, sep = "_")

  # Quality Control: Check if file_slug is valid
  if (is.na(file_slug) || file_slug == "NA_NA") {
    cat("[WARNING] Skipping save: Invalid file slug.\n")
    next
  }

  # Save the plot
  plot_filename <- file.path(plots_path, paste0("twfe_plot_", file_slug, ".pdf"))
  ggsave(plot_filename, plot = p, width = 8, height = 6)
  cat(paste0("Saved plot to: ", plot_filename, "\n"))

  # Save the table (tidy coefficients)
  table_filename <- file.path(tables_path, paste0("twfe_results_", file_slug, ".csv"))
  readr::write_csv(results_for_plot, table_filename)
  cat(paste0("Saved table to: ", table_filename, "\n"))

  # Optional quick check in console
  cat("\n— Coefficients at event times —\n")
  print(results_for_plot %>% select(event_time, relative_time, estimate, conf.low, conf.high))

  # Accumulate results for parallel trends analysis
  plot_data <- bind_rows(plot_data, results_for_plot %>%
    mutate(
      funder_abbr = funder_abbr,
      funder_full = funder_full,
      sector_code = sector_code,
      ad_sector_names = sector_full
    ))
}

############### PARALLELL TRENDS ASSUMPTION ANALYSIS ########
################################################################################
# 1. SETUP: Ensure data from the previous script is available
################################################################################

# This script assumes 'plot_data' from the script above has been created.
# If running this separately, you must run the "DATA PREPARATION" section
# of the first script before this one.


################################################################################
# 2. ANALYSIS: Isolate Pre-Treatment Coefficients
################################################################################

# --- Define save paths (using here::here for portability) ---
plots_save_path <- here::here("output", "TWFE", "plots")
tables_save_path <- here::here("output", "TWFE", "tables")
dir.create(plots_save_path, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_save_path, recursive = TRUE, showWarnings = FALSE)

# Filter for only the pre-treatment period (event_time = -1)
pre_trend_data <- plot_data %>%
  filter(event_time == -1) %>%
  # Create a clean panel name for the plot's y-axis
  mutate(
    panel_name = str_glue("{funder_full} - {ad_sector_names}"),
    # Create a variable to check if the estimate is statistically significant
    # (i.e., if the confidence interval does NOT contain zero)
    is_significant = !(conf.low < 0 & conf.high > 0)
  ) %>%
  # Order the panels by the size of the estimate for a cleaner look
  mutate(
    panel_name = fct_reorder(panel_name, estimate)
  )

# Quick check of the resulting data
glimpse(pre_trend_data)


################################################################################
# 3. VISUALIZATION: Create the Parallel Trends Dot-Whisker Plot
################################################################################

# Create the plot
parallel_trends_plot <- ggplot(pre_trend_data, aes(x = estimate, y = panel_name, color = is_significant)) +
  # Add a vertical line at x=0. This is the crucial reference line.
  geom_vline(xintercept = 0, color = "gray40", linetype = "dashed", linewidth = 1) +

  # Add the point estimate and its confidence interval as a horizontal range
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.4) +

  # Use color to distinguish between significant and non-significant pre-trends
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    labels = c("Not Significant", "Significant (p < 0.05)"),
    name = "Pre-Trend Violation"
  ) +

  # Set labels and title
  labs(
    title = "Parallel Trends Assumption",
    subtitle = "Pre-Treatment Coefficient (t-1) for Each Funder-Sector Panel",
    x = "Estimated Effect on IWI",
    y = "" # The panel names on the y-axis are self-explanatory
  ) +

  # Apply a clean theme
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines for clarity
    plot.title.position = "plot",
    plot.title = element_text(face = "bold")
  )

# --- Display the plot for inspection ---
print(parallel_trends_plot)

# --- Save the plot ---
ggsave(
  filename = file.path(plots_save_path, "Fig_TWFE_Parallel_Trends_Test.png"),
  plot = parallel_trends_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# --- Save the table ---
write_csv(
  pre_trend_data,
  file = file.path(tables_save_path, "Table_TWFE_Parallel_Trends_Test.csv")
)

cat("Successfully saved parallel trends plot to:", plots_save_path, "\n")
cat("Successfully saved parallel trends table to:", tables_save_path, "\n")
