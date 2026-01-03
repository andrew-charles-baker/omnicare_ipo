# =============================================================================
# 02_pre_post_summaries.R
# Pre vs Post Omnicare summary visualizations
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# --- Load Data ---
df <- load_ipo_data()

df_complete <- df[complete.cases(df[, c("first_day_return", "rf_fact_t07",
                                         "log_assets", "log_proceeds",
                                         "roll_vwretd", "ff12")]), ]

df_pre <- df[!df$post_omnicare, ]
df_post <- df[df$post_omnicare, ]

cat("Pre-Omnicare:", nrow(df_pre), "\n")
cat("Post-Omnicare:", nrow(df_post), "\n")

# =============================================================================
# 1. CDF COMPARISONS
# =============================================================================

# --- 1a. CDF of rf_fact: Pre vs Post ---
p_cdf_rf <- ggplot() +
  stat_ecdf(data = df_pre, aes(x = rf_fact_t07, color = "Pre-Omnicare"), linewidth = 1.2) +
  stat_ecdf(data = df_post, aes(x = rf_fact_t07, color = "Post-Omnicare"), linewidth = 1.2) +
  scale_color_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  labs(
    title = "CDF of Fact-Intensity: Pre vs Post Omnicare",
    subtitle = sprintf("Pre mean: %.3f | Post mean: %.3f",
                       mean(df_pre$rf_fact_t07, na.rm = TRUE),
                       mean(df_post$rf_fact_t07, na.rm = TRUE)),
    x = "Fact-Intensity (rf_fact_t07)",
    y = "Cumulative Probability",
    color = "Period"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

save_plot(p_cdf_rf, "02a_cdf_rf_fact.png")

# --- 1b. CDF of First-Day Returns: Pre vs Post ---
p_cdf_returns <- ggplot() +
  stat_ecdf(data = df_pre, aes(x = first_day_return, color = "Pre-Omnicare"), linewidth = 1.2) +
  stat_ecdf(data = df_post, aes(x = first_day_return, color = "Post-Omnicare"), linewidth = 1.2) +
  scale_color_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  coord_cartesian(xlim = c(-0.5, 1.5)) +
  labs(
    title = "CDF of First-Day Returns: Pre vs Post Omnicare",
    subtitle = sprintf("Pre mean: %.1f%% | Post mean: %.1f%%",
                       100 * mean(df_pre$first_day_return, na.rm = TRUE),
                       100 * mean(df_post$first_day_return, na.rm = TRUE)),
    x = "First-Day Return",
    y = "Cumulative Probability",
    color = "Period"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

save_plot(p_cdf_returns, "02b_cdf_returns.png")

# --- Combined CDF plot ---
p_cdf_combined <- p_cdf_rf + p_cdf_returns +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(p_cdf_combined, "02c_cdf_combined.png", width = 14, height = 6)


# =============================================================================
# 2. PDF/HISTOGRAM COMPARISONS
# =============================================================================

# --- 2a. Histogram of rf_fact ---
p_hist_rf <- ggplot() +
  geom_histogram(data = df_pre, aes(x = rf_fact_t07, y = after_stat(density), fill = "Pre-Omnicare"),
                 alpha = 0.5, bins = 30, color = "white") +
  geom_histogram(data = df_post, aes(x = rf_fact_t07, y = after_stat(density), fill = "Post-Omnicare"),
                 alpha = 0.5, bins = 30, color = "white") +
  geom_density(data = df_pre, aes(x = rf_fact_t07, color = "Pre-Omnicare"), linewidth = 1) +
  geom_density(data = df_post, aes(x = rf_fact_t07, color = "Post-Omnicare"), linewidth = 1) +
  scale_fill_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  scale_color_manual(values = c("Pre-Omnicare" = "darkred", "Post-Omnicare" = "darkblue")) +
  labs(
    title = "Distribution of Fact-Intensity: Pre vs Post Omnicare",
    x = "Fact-Intensity (rf_fact_t07)",
    y = "Density",
    fill = "Period",
    color = "Period"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5)))

save_plot(p_hist_rf, "02d_histogram_rf_fact.png")

# --- 2b. Histogram of First-Day Returns ---
p_hist_returns <- ggplot() +
  geom_histogram(data = df_pre, aes(x = first_day_return, y = after_stat(density), fill = "Pre-Omnicare"),
                 alpha = 0.5, bins = 50, color = "white") +
  geom_histogram(data = df_post, aes(x = first_day_return, y = after_stat(density), fill = "Post-Omnicare"),
                 alpha = 0.5, bins = 50, color = "white") +
  geom_density(data = df_pre, aes(x = first_day_return, color = "Pre-Omnicare"), linewidth = 1) +
  geom_density(data = df_post, aes(x = first_day_return, color = "Post-Omnicare"), linewidth = 1) +
  scale_fill_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  scale_color_manual(values = c("Pre-Omnicare" = "darkred", "Post-Omnicare" = "darkblue")) +
  coord_cartesian(xlim = c(-0.5, 2)) +
  labs(
    title = "Distribution of First-Day Returns: Pre vs Post Omnicare",
    x = "First-Day Return",
    y = "Density",
    fill = "Period",
    color = "Period"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5)))

save_plot(p_hist_returns, "02e_histogram_returns.png")

# --- Combined histogram plot ---
p_hist_combined <- p_hist_rf + p_hist_returns +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(p_hist_combined, "02f_histogram_combined.png", width = 14, height = 6)


# =============================================================================
# 3. DISTRIBUTION REGRESSION BY PERIOD
# =============================================================================

cat("Running distribution regression for pre-Omnicare...\n")
dr_pre <- distribution_regression(df_pre, "rf_fact_t07", "first_day_return",
                                   n_thresholds = 20, n_bootstrap = 200)
dr_pre$period <- "Pre-Omnicare"

cat("Running distribution regression for post-Omnicare...\n")
dr_post <- distribution_regression(df_post, "rf_fact_t07", "first_day_return",
                                    n_thresholds = 25, n_bootstrap = 200)
dr_post$period <- "Post-Omnicare"

dr_periods <- rbind(dr_pre, dr_post)

p_dr_periods <- ggplot(dr_periods, aes(x = threshold, y = marginal_effect,
                                        color = period, fill = period)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = me_lower, ymax = me_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  scale_fill_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  labs(
    title = "Distribution Regression: Pre vs Post Omnicare (No Controls)",
    subtitle = "Marginal effect of rf_fact on P(Return > threshold)",
    x = "First-Day Return Threshold",
    y = "Marginal Effect",
    color = "Period",
    fill = "Period"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

save_plot(p_dr_periods, "02g_dist_reg_pre_post.png")


# =============================================================================
# 4. DISTRIBUTION REGRESSION BY PERIOD (WITH CONTROLS)
# =============================================================================

df_complete_pre <- df_complete[!df_complete$post_omnicare, ]
df_complete_post <- df_complete[df_complete$post_omnicare, ]

controls <- c("log_assets", "log_proceeds", "roll_vwretd", "ff12")

cat("Running distribution regression with controls for pre-Omnicare...\n")
dr_pre_ctrl <- distribution_regression_controls(
  df_complete_pre, "rf_fact_t07", "first_day_return", controls,
  n_thresholds = 20, n_bootstrap = 200
)
dr_pre_ctrl$period <- "Pre-Omnicare"

cat("Running distribution regression with controls for post-Omnicare...\n")
dr_post_ctrl <- distribution_regression_controls(
  df_complete_post, "rf_fact_t07", "first_day_return", controls,
  n_thresholds = 25, n_bootstrap = 200
)
dr_post_ctrl$period <- "Post-Omnicare"

dr_periods_ctrl <- rbind(dr_pre_ctrl, dr_post_ctrl)

p_dr_periods_ctrl <- ggplot(dr_periods_ctrl, aes(x = threshold, y = marginal_effect,
                                                   color = period, fill = period)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = me_lower, ymax = me_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  scale_fill_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
  labs(
    title = "Distribution Regression: Pre vs Post Omnicare (With Controls)",
    subtitle = "Controls: log(assets), log(proceeds), roll_vwretd, FF12 industry",
    x = "First-Day Return Threshold",
    y = "Marginal Effect",
    color = "Period",
    fill = "Period"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

save_plot(p_dr_periods_ctrl, "02h_dist_reg_pre_post_controls.png")

# --- Combined: With and Without Controls ---
dr_periods$specification <- "No Controls"
dr_periods_ctrl$specification <- "With Controls"

p_dr_periods_comparison <- (p_dr_periods + labs(title = "No Controls")) +
  (p_dr_periods_ctrl + labs(title = "With Controls")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(p_dr_periods_comparison, "02i_dist_reg_pre_post_comparison.png", width = 14, height = 6)


# =============================================================================
# 5. SUMMARY STATISTICS TABLE
# =============================================================================

summary_stats <- df %>%
  group_by(post_omnicare) %>%
  summarise(
    n = n(),
    mean_rf_fact = mean(rf_fact_t07, na.rm = TRUE),
    sd_rf_fact = sd(rf_fact_t07, na.rm = TRUE),
    mean_return = mean(first_day_return, na.rm = TRUE),
    sd_return = sd(first_day_return, na.rm = TRUE),
    median_return = median(first_day_return, na.rm = TRUE),
    mean_assets = mean(assets_thous, na.rm = TRUE),
    mean_proceeds = mean(total_proceeds_thous, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(period = ifelse(post_omnicare, "Post-Omnicare", "Pre-Omnicare"))

cat("\n=== SUMMARY STATISTICS ===\n")
print(summary_stats)

# Save summary stats
write.csv(summary_stats, file.path(OUTPUT_DIR, "02_summary_statistics.csv"), row.names = FALSE)

cat("\n=== Script 02 Complete ===\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
