################################################################################
# Spatial Spillovers & Nearby Neighborhoods Robustness
# Mattias Antar, 2025
# (1) Distance-based dose–response
# (2) Exclusion buffers
# (3) Spatial autocorrelation diagnostics
################################################################################

# ---- Packages -----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, broom, geosphere, sandwich, lmtest, spdep, ggplot2, here)

# ---- Paths (using here::here for portability) ---------------------------------
data_path <- here::here("Data", "Archive_enriched")
save_path <- here::here("output", "spillover")
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

# ---- Locate & read inputs -----------------------------------------------------
files <- list.files(
  data_path,
  pattern = "InputData_.*_DiD_enriched\\.csv$",
  full.names = TRUE
)
stopifnot(length(files) > 0)

read_one <- function(f) {
  suppressMessages(readr::read_csv(f, show_col_types = FALSE)) %>%
    clean_names()
}
raw <- bind_rows(lapply(files, read_one))

# ---- Minimal sanity checks ----------------------------------------------------
need <- c(
  "lat", "lon", "treated", "time", "adm2", "iso3",
  "iwi_2014_2016", "iwi_2017_2019",
  "log_avg_min_to_city", "log_avg_pop_dens", "log_3yr_pre_conflict_deaths"
)
missing <- setdiff(need, names(raw))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

# ---- Collapse to unit-level (ADM2) -------------------------------------------
raw <- raw %>%
  mutate(
    time = as.integer(time),
    treated = as.integer(treated),
    cluster_id = paste0(iso3, "_", adm2)
  )

unit_level <- raw %>%
  group_by(cluster_id) %>%
  summarise(
    iso3 = first(iso3),
    adm2 = first(adm2),
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    treated = as.integer(any(treated == 1 & time <= 2012, na.rm = TRUE)),
    iwi_2014_2016 = mean(iwi_2014_2016, na.rm = TRUE),
    iwi_2017_2019 = mean(iwi_2017_2019, na.rm = TRUE),
    log_avg_min_to_city = mean(log_avg_min_to_city, na.rm = TRUE),
    log_avg_pop_dens = mean(log_avg_pop_dens, na.rm = TRUE),
    log_3yr_pre_conflict_deaths = mean(log_3yr_pre_conflict_deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    country = iso3,
    iwi_est_post_oda = coalesce(
      rowMeans(cbind(iwi_2014_2016, iwi_2017_2019), na.rm = TRUE),
      iwi_2017_2019, iwi_2014_2016
    )
  ) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(iwi_est_post_oda))

combined <- unit_level %>%
  mutate(time = 2018L) %>%
  relocate(cluster_id, iso3, adm2, country, time, treated, lat, lon)

# ---- Distances to nearest treated cluster -------------------------------------
centroids <- combined %>%
  distinct(cluster_id, .keep_all = TRUE) %>%
  select(cluster_id, lat, lon, treated)
treated_centroids <- centroids %>% filter(treated == 1)

if (nrow(treated_centroids) == 0) {
  warning("No treated clusters detected; distance = Inf")
  centroids$dist_to_treated_km <- Inf
} else {
  dist_matrix <- geosphere::distm(
    as.matrix(centroids[, c("lon", "lat")]),
    as.matrix(treated_centroids[, c("lon", "lat")]),
    fun = geosphere::distHaversine
  ) / 1000
  centroids$dist_to_treated_km <- apply(dist_matrix, 1, min)
}

combined <- combined %>%
  left_join(centroids %>% select(cluster_id, dist_to_treated_km), by = "cluster_id")

breaks <- c(0, 10, 20, 50, Inf)
labels <- c("[0,10)", "[10,20)", "[20,50)", "[50,∞)")
combined <- combined %>%
  mutate(
    dist_band = cut(dist_to_treated_km, breaks, right = FALSE, labels = labels),
    dist_band = forcats::fct_relevel(dist_band, "[50,∞)")
  )

# ---- Model-ready analytic sample ----------------------------------------------
analytic <- combined %>%
  select(
    iwi_est_post_oda, treated, dist_band, dist_to_treated_km,
    log_avg_min_to_city, log_avg_pop_dens, log_3yr_pre_conflict_deaths,
    time, country, lat, lon
  ) %>%
  drop_na() %>%
  mutate(time = as.factor(time), country = as.factor(country), .row_id = row_number())

# ==============================================================================
# (1) Distance-based dose–response
# ==============================================================================
m_band <- lm(iwi_est_post_oda ~ treated + dist_band +
  log_avg_min_to_city + log_avg_pop_dens + log_3yr_pre_conflict_deaths +
  country, data = analytic)
vcov_hc_band <- sandwich::vcovHC(m_band, type = "HC2")
band_coeftest <- lmtest::coeftest(m_band, vcov. = vcov_hc_band)
band_tidy <- broom::tidy(band_coeftest) %>%
  mutate(term = rownames(band_coeftest)) %>%
  select(term, estimate, std.error, statistic, p.value)

# ---- Plot 1: Predicted IWI gradient by distance band --------------------------
pred_data <- broom::augment(m_band, analytic) %>%
  group_by(dist_band) %>%
  summarise(
    mean_pred = mean(.fitted, na.rm = TRUE),
    mean_obs = mean(iwi_est_post_oda, na.rm = TRUE), .groups = "drop"
  )

fig1 <- ggplot(pred_data, aes(x = dist_band)) +
  geom_col(aes(y = mean_pred), fill = "#1b9e77", width = 0.6, alpha = 0.8) +
  geom_point(aes(y = mean_obs), color = "black", size = 3) +
  geom_line(aes(y = mean_obs, group = 1), color = "black", linewidth = 0.8) +
  labs(
    title = "Predicted IWI Gradient by Distance to Treated Cluster",
    x = "Distance Band (km)", y = "Mean Predicted IWI"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(save_path, "Fig1_IWI_gradient_by_distance.png"), fig1,
  width = 7, height = 5, dpi = 300
)

# ==============================================================================
# (2) Exclusion buffers
# ==============================================================================
Ks <- c(0, 10, 20, 30)
buffer_results <- purrr::map_dfr(Ks, function(K) {
  df_k <- analytic %>% filter(dist_to_treated_km >= K | treated == 1)
  m_k <- lm(iwi_est_post_oda ~ treated +
    log_avg_min_to_city + log_avg_pop_dens + log_3yr_pre_conflict_deaths +
    country, data = df_k)
  V_k <- sandwich::vcovHC(m_k, type = "HC2")
  ct <- lmtest::coeftest(m_k, vcov. = V_k)
  tibble(
    K_buffer_km = K,
    estimate = ct["treated", "Estimate"],
    std.error = ct["treated", "Std. Error"],
    statistic = ct["treated", "t value"],
    p.value = ct["treated", "Pr(>|t|)"],
    n_obs = nobs(m_k)
  )
})

# ---- Plot 2: Treated coefficient vs buffer radius -----------------------------
fig2 <- ggplot(buffer_results, aes(x = K_buffer_km, y = estimate)) +
  geom_line(color = "#d95f02", linewidth = 1) +
  geom_point(size = 3, color = "#d95f02") +
  geom_errorbar(
    aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ),
    width = 1.5, color = "gray40"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Effect of Excluding Controls Near Treated Clusters",
    x = "Exclusion Buffer (km)", y = "Estimated Treatment Effect (β̂)"
  )

ggsave(file.path(save_path, "Fig2_TreatedEffect_vs_Buffer.png"), fig2,
  width = 7, height = 5, dpi = 300
)

# ==============================================================================
# (3) Spatial autocorrelation diagnostics
# ==============================================================================
resid_vec <- residuals(m_band)
coords_used <- analytic %>% select(lon, lat)
k_nn <- 5
knn_nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords_used), k = k_nn))
listw <- spdep::nb2listw(knn_nb, style = "W")
moran_global <- spdep::moran.test(resid_vec, listw)
moran_mc <- spdep::moran.mc(resid_vec, listw, nsim = 999)

# ---- Plot 3: Spatial map of residuals -----------------------------------------
coords_with_resid <- analytic %>%
  mutate(
    residual = resid_vec,
    resid_std = scale(residual)[, 1]
  )
fig3 <- ggplot(coords_with_resid, aes(lon, lat, color = resid_std)) +
  geom_point(size = 1.4, alpha = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Residual (std)"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Spatial Distribution of Residuals (Moran’s I Map)",
    subtitle = sprintf(
      "Global Moran’s I = %.3f (p < %.3g)",
      moran_global$estimate[1], moran_global$p.value
    ),
    x = "Longitude", y = "Latitude"
  )

ggsave(file.path(save_path, "Fig3_Residuals_Moran_Map.png"), fig3,
  width = 7, height = 5, dpi = 300
)

# ==============================================================================
# Save Results
# ==============================================================================
write_csv(band_tidy, file.path(save_path, "Table1_band_tidy.csv"))
write_csv(buffer_results, file.path(save_path, "Table2_buffer_results.csv"))

# Save Moran’s I summary
moran_summary <- tibble(
  moran_i = moran_global$estimate[1],
  expected = moran_global$estimate[2],
  variance = moran_global$estimate[3],
  p_value = moran_global$p.value
)
write_csv(moran_summary, file.path(save_path, "Table3_Moran_summary.csv"))

# Save session metadata
sessionInfo() %>%
  capture.output(file = file.path(save_path, "R_sessionInfo.txt"))

# ---- Console confirmation -----------------------------------------------------
cat("\nAll spillover robustness outputs saved to:\n", save_path, "\n")
cat("\nSaved files:\n")
print(list.files(save_path, full.names = TRUE))
