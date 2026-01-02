#!/usr/bin/env Rscript
# =============================================================================
# dCdH Difference-in-Differences Analysis with 6 Covariates (Unified Version)
# =============================================================================
# This script runs the de Chaisemartin & D'Haultfoeuille (dCdH) estimator
# with covariates on all 23 panel datasets. Includes fallback strategies
# for panels that fail due to collinearity or insufficient variation.
#
# Usage:
#   Rscript dcdh_with_covariates.R              # Run all panels
#   Rscript dcdh_with_covariates.R --test       # Test on single panel
#
# Fallback Strategies (tried in order when covariates fail):
#   1. effects=3, no covariates
#   2. effects=3, reduced covariates
#   3. effects=2, no covariates
#   4. effects=1, no covariates
# =============================================================================

# --- Configuration -----------------------------------------------------------

# Handle Mac-specific issues with rgl package
Sys.setenv(RGL_USE_NULL = TRUE)

# Load here package for relative paths
if (!require("here", quietly = TRUE)) {
    install.packages("here", repos = "https://cloud.r-project.org")
    library(here)
}

# Set paths using here::here() for portability
DATA_DIR <- here::here("Data", "Archive_enriched")
OUTPUT_DIR <- here::here("output", "dCdH")

# Create output directory if needed
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Check if running in test mode
args <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args

# --- Package Installation/Loading --------------------------------------------

cat("=== dCdH Analysis with Covariates (Unified) ===\n\n")

# Try to load the package, install if needed
if (!require("DIDmultiplegtDYN", quietly = TRUE)) {
    cat("Installing DIDmultiplegtDYN package...\n")

    # First try CRAN
    tryCatch(
        {
            install.packages("DIDmultiplegtDYN", repos = "https://cloud.r-project.org")
        },
        error = function(e) {
            cat("CRAN install failed. Trying GitHub...\n")
            if (!require("devtools", quietly = TRUE)) {
                install.packages("devtools", repos = "https://cloud.r-project.org")
            }
            devtools::install_github("chaisemartinPackages/did_multiplegt_dyn/R")
        }
    )

    library(DIDmultiplegtDYN)
}

cat("DIDmultiplegtDYN package loaded successfully!\n\n")

# --- Define Covariates -------------------------------------------------------

# Full set of 6 covariates
COVARIATES_FULL <- c(
    "log_avg_pop_dens", # Population density
    "log_3yr_pre_conflict_deaths", # Conflict
    "log_disasters", # Natural disasters
    "election_year", # Election cycles
    "political_stability", # Institutional stability
    "leader_birthplace" # Leader birthplace
)

# Reduced set for collinearity issues
COVARIATES_REDUCED <- c(
    "log_avg_pop_dens",
    "log_disasters",
    "political_stability"
)

cat("Full covariates:\n")
for (cov in COVARIATES_FULL) cat("  -", cov, "\n")
cat("\n")

# --- Define Panels -----------------------------------------------------------

# List all panel CSV files
panel_files <- list.files(DATA_DIR, pattern = "InputData_.*_DiD_enriched\\.csv$", full.names = TRUE)

if (length(panel_files) == 0) {
    stop("No panel files found in: ", DATA_DIR)
}

cat("Found", length(panel_files), "panel files.\n\n")

# If test mode, only run on first panel
if (TEST_MODE) {
    cat("** TEST MODE: Running on first panel only **\n\n")
    panel_files <- panel_files[1]
}

# --- Helper Functions --------------------------------------------------------

extract_panel_info <- function(filepath) {
    # Extract funder and sector from filename
    basename_file <- basename(filepath)
    parts <- strsplit(gsub("InputData_|_DiD_enriched\\.csv", "", basename_file), "_")[[1]]
    list(funder = parts[1], sector = parts[2])
}

run_dcdh_with_fallback <- function(filepath) {
    # Run dCdH estimator with fallback strategies if covariates fail

    info <- extract_panel_info(filepath)
    panel_name <- paste(toupper(info$funder), info$sector, sep = "-")
    cat(sprintf("\n%s\n", strrep("=", 60)))
    cat(sprintf("Processing: %s\n", panel_name))
    cat(sprintf("%s\n", strrep("=", 60)))

    # Load data
    df <- read.csv(filepath, stringsAsFactors = FALSE)
    cat(sprintf("Loaded %d observations, %d groups\n", nrow(df), length(unique(df$dhs_id))))
    cat(sprintf("Treated: %d, Control: %d\n", sum(df$treated == 1), sum(df$treated == 0)))

    # Check for required columns
    required_cols <- c("dhs_id", "time", "treated", "iwi_est_post_oda")
    missing_cols <- setdiff(required_cols, names(df))

    if (length(missing_cols) > 0) {
        cat("  WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
        return(NULL)
    }

    # Check if there are any switchers
    if (sum(df$treated == 1) == 0) {
        cat("  WARNING: No treated units in this panel. Skipping.\n")
        return(NULL)
    }

    # Define fallback strategies - prioritize effects=3 with covariates
    strategies <- list(
        list(
            name = "Strategy 1: effects=3, full covariates",
            effects = 3, placebo = 1, controls = COVARIATES_FULL
        ),
        list(
            name = "Strategy 2: effects=3, no covariates",
            effects = 3, placebo = 1, controls = NULL
        ),
        list(
            name = "Strategy 3: effects=3, reduced covariates",
            effects = 3, placebo = 1, controls = COVARIATES_REDUCED
        ),
        list(
            name = "Strategy 4: effects=2, no covariates",
            effects = 2, placebo = 1, controls = NULL
        ),
        list(
            name = "Strategy 5: effects=1, no covariates",
            effects = 1, placebo = 1, controls = NULL
        )
    )

    for (strat in strategies) {
        cat(sprintf("\nTrying: %s\n", strat$name))

        tryCatch(
            {
                result <- did_multiplegt_dyn(
                    df = df,
                    outcome = "iwi_est_post_oda",
                    group = "dhs_id",
                    time = "time",
                    treatment = "treated",
                    effects = strat$effects,
                    placebo = strat$placebo,
                    controls = strat$controls,
                    cluster = "dhs_id",
                    graph_off = TRUE
                )

                cat("SUCCESS!\n\n")

                # Extract results
                r <- result$results
                row <- data.frame(
                    funder = info$funder,
                    sector = info$sector,
                    strategy = strat$name,
                    covariates_used = ifelse(is.null(strat$controls), "None",
                        ifelse(length(strat$controls) == 6, "Full (6)", "Reduced (3)")
                    ),
                    stringsAsFactors = FALSE
                )

                # Extract effects from Effects matrix
                effects_mat <- r$Effects
                if (!is.null(effects_mat) && nrow(effects_mat) >= 1) {
                    row$effect_1 <- effects_mat[1, "Estimate"]
                    row$se_1 <- effects_mat[1, "SE"]
                    row$ci_lower_1 <- effects_mat[1, "LB CI"]
                    row$ci_upper_1 <- effects_mat[1, "UB CI"]
                }
                if (!is.null(effects_mat) && nrow(effects_mat) >= 2) {
                    row$effect_2 <- effects_mat[2, "Estimate"]
                    row$se_2 <- effects_mat[2, "SE"]
                    row$ci_lower_2 <- effects_mat[2, "LB CI"]
                    row$ci_upper_2 <- effects_mat[2, "UB CI"]
                }
                if (!is.null(effects_mat) && nrow(effects_mat) >= 3) {
                    row$effect_3 <- effects_mat[3, "Estimate"]
                    row$se_3 <- effects_mat[3, "SE"]
                    row$ci_lower_3 <- effects_mat[3, "LB CI"]
                    row$ci_upper_3 <- effects_mat[3, "UB CI"]
                }

                # Placebo
                placebos_mat <- r$Placebos
                if (!is.null(placebos_mat) && nrow(placebos_mat) >= 1) {
                    row$placebo <- placebos_mat[1, "Estimate"]
                    row$placebo_se <- placebos_mat[1, "SE"]
                    row$placebo_ci_lower <- placebos_mat[1, "LB CI"]
                    row$placebo_ci_upper <- placebos_mat[1, "UB CI"]
                }

                return(row)
            },
            error = function(e) {
                cat(sprintf("  FAILED: %s\n", conditionMessage(e)))
            }
        )
    }

    # All strategies failed
    cat("\nAll strategies failed for this panel.\n")
    return(NULL)
}

# --- Main Execution ----------------------------------------------------------

cat("Starting estimation...\n")
cat(strrep("=", 60), "\n")

results_list <- list()

for (filepath in panel_files) {
    result <- run_dcdh_with_fallback(filepath)

    if (!is.null(result)) {
        results_list[[length(results_list) + 1]] <- result
    }
}

cat(sprintf("\n%s\n", strrep("=", 60)))
cat(sprintf(
    "Completed: %d/%d panels successfully estimated\n\n",
    length(results_list), length(panel_files)
))

# --- Compile Results ---------------------------------------------------------

if (length(results_list) > 0) {
    cat("Compiling results...\n")

    # Combine all result rows
    results_df <- do.call(rbind, results_list)

    # Save results
    output_file <- file.path(OUTPUT_DIR, "dcdh_results.csv")
    write.csv(results_df, output_file, row.names = FALSE)
    cat(sprintf("\nResults saved to: %s\n", output_file))

    # Print summary
    cat("\n=== Results Summary ===\n\n")

    # Print strategy usage
    cat("Strategy usage:\n")
    print(table(results_df$strategy))

    cat("\n")

    # Print covariate usage
    cat("Covariate usage:\n")
    print(table(results_df$covariates_used))

    cat("\n")

    # Print key results
    print_cols <- intersect(c("funder", "sector", "covariates_used", "effect_1", "effect_2", "effect_3", "placebo"), names(results_df))
    print(results_df[, print_cols, drop = FALSE])
} else {
    cat("No results to compile.\n")
}

cat("\nDone!\n")
