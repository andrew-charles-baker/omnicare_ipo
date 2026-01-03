# =============================================================================
# 01_relationship_plots.R
# Plots showing overall relationship between rf_fact and IPO underpricing
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# --- Load Data ---
df <- load_ipo_data()

# Create complete cases dataset for controlled analyses
df_complete <- df[complete.cases(df[, c("first_day_return", "rf_fact_t07",
                                         "log_assets", "log_proceeds",
                                         "roll_vwretd", "ff12")]), ]

cat("Full sample:", nrow(df), "\n")
cat("Complete cases:", nrow(df_complete), "\n")

# =============================================================================
# 1. BINSCATTER PLOTS
# =============================================================================

# --- 1a. Binscatter without controls (using binsreg or fallback) ---
create_binscatter <- function(data, x_var, y_var, n_bins = 20, title = NULL) {

  if (binsreg_available) {
    # Use binsreg for optimal binning
    bins_result <- binsreg(
      y = data[[y_var]],
      x = data[[x_var]],
      nbins = n_bins,
      polyreg = 1
    )

    bin_data <- data.frame(
      x = bins_result$data.plot$`Group Full Sample`$data.dots$x,
      y = bins_result$data.plot$`Group Full Sample`$data.dots$fit
    )
  } else {
    # Fallback: quantile-based binning
    data$bin <- cut(data[[x_var]],
                    breaks = quantile(data[[x_var]], probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
                    include.lowest = TRUE, labels = FALSE)

    bin_data <- data %>%
      group_by(bin) %>%
      summarise(
        x = mean(.data[[x_var]], na.rm = TRUE),
        y = mean(.data[[y_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(!is.na(x))
  }

  # Fit line
  fit <- lm(y ~ x, data = bin_data)

  p <- ggplot(bin_data, aes(x = x, y = y)) +
    geom_point(size = 3, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", fill = "darkred", alpha = 0.2) +
    labs(
      title = title %||% paste("Binscatter:", y_var, "vs", x_var),
      x = "Fact-Intensity (rf_fact_t07)",
      y = "First-Day Return",
      caption = sprintf("Slope: %.4f", coef(fit)[2])
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  return(p)
}

# Binscatter without controls
p_binscatter_baseline <- create_binscatter(
  df, "rf_fact_t07", "first_day_return",
  title = "Binscatter: IPO Underpricing vs Fact-Intensity (No Controls)"
)
save_plot(p_binscatter_baseline, "01a_binscatter_baseline.png")

# --- 1b. Binscatter with controls (Frisch-Waugh-Lovell residualization) ---
create_binscatter_controls <- function(data, x_var, y_var, controls, n_bins = 20, title = NULL) {

  # Residualize x
  formula_x <- as.formula(paste(x_var, "~", paste(controls, collapse = " + ")))
  model_x <- lm(formula_x, data = data)
  data$x_resid <- residuals(model_x)

  # Residualize y
  formula_y <- as.formula(paste(y_var, "~", paste(controls, collapse = " + ")))
  model_y <- lm(formula_y, data = data)
  data$y_resid <- residuals(model_y)

  # Create bins on residualized x
  data$bin <- cut(data$x_resid,
                  breaks = quantile(data$x_resid, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
                  include.lowest = TRUE, labels = FALSE)

  bin_data <- data %>%
    group_by(bin) %>%
    summarise(
      x = mean(x_resid, na.rm = TRUE),
      y = mean(y_resid, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(x))

  fit <- lm(y ~ x, data = bin_data)

  p <- ggplot(bin_data, aes(x = x, y = y)) +
    geom_point(size = 3, color = "darkgreen") +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", fill = "darkred", alpha = 0.2) +
    labs(
      title = title %||% "Binscatter with Controls (FWL Residualized)",
      subtitle = paste("Controls:", paste(controls, collapse = ", ")),
      x = "Residualized Fact-Intensity",
      y = "Residualized First-Day Return",
      caption = sprintf("Slope: %.4f", coef(fit)[2])
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  return(p)
}

controls <- c("log_assets", "log_proceeds", "roll_vwretd", "ff12")
p_binscatter_controls <- create_binscatter_controls(
  df_complete, "rf_fact_t07", "first_day_return", controls,
  title = "Binscatter: IPO Underpricing vs Fact-Intensity (With Controls)"
)
save_plot(p_binscatter_controls, "01b_binscatter_controls.png")

# --- Combined binscatter plot ---
p_binscatter_combined <- p_binscatter_baseline + p_binscatter_controls +
  plot_annotation(title = "Binscatter Analysis: rf_fact vs First-Day Returns",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_plot(p_binscatter_combined, "01c_binscatter_combined.png", width = 14, height = 6)


# =============================================================================
# 2. HEATMAP PLOTS
# =============================================================================

create_heatmap <- function(data, x_var, y_var, n_bins_x = 20, n_bins_y = 20, title = NULL) {

  # Create bins
  data$x_bin <- cut(data[[x_var]],
                    breaks = quantile(data[[x_var]], probs = seq(0, 1, length.out = n_bins_x + 1), na.rm = TRUE),
                    include.lowest = TRUE)

  data$y_bin <- cut(data[[y_var]],
                    breaks = quantile(data[[y_var]], probs = seq(0, 1, length.out = n_bins_y + 1), na.rm = TRUE),
                    include.lowest = TRUE)

  # Count observations in each cell
  heatmap_data <- data %>%
    filter(!is.na(x_bin) & !is.na(y_bin)) %>%
    group_by(x_bin, y_bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      x_mid = as.numeric(x_bin),
      y_mid = as.numeric(y_bin)
    )

  p <- ggplot(heatmap_data, aes(x = x_mid, y = y_mid, fill = count)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", name = "Count") +
    labs(
      title = title %||% "Heatmap: Joint Distribution",
      x = "Fact-Intensity Quantile",
      y = "First-Day Return Quantile"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

p_heatmap <- create_heatmap(
  df, "rf_fact_t07", "first_day_return",
  title = "Heatmap: rf_fact vs First-Day Returns (Joint Distribution)"
)
save_plot(p_heatmap, "01d_heatmap.png")

# --- Heatmap with residualized data ---
df_resid <- df_complete
model_x <- lm(rf_fact_t07 ~ log_assets + log_proceeds + roll_vwretd + ff12, data = df_resid)
model_y <- lm(first_day_return ~ log_assets + log_proceeds + roll_vwretd + ff12, data = df_resid)
df_resid$rf_fact_resid <- residuals(model_x)
df_resid$return_resid <- residuals(model_y)

p_heatmap_controls <- create_heatmap(
  df_resid, "rf_fact_resid", "return_resid",
  title = "Heatmap: rf_fact vs Returns (Residualized, With Controls)"
)
save_plot(p_heatmap_controls, "01e_heatmap_controls.png")


# =============================================================================
# 3. DISTRIBUTION REGRESSION PLOTS
# =============================================================================

cat("Running distribution regression (baseline)...\n")
dr_baseline <- distribution_regression(df, "rf_fact_t07", "first_day_return",
                                        n_thresholds = 25, n_bootstrap = 200)

p_dr_baseline <- ggplot(dr_baseline, aes(x = threshold, y = marginal_effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = me_lower, ymax = me_upper), alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  labs(
    title = "Distribution Regression: Baseline (No Controls)",
    subtitle = "Marginal effect of rf_fact on P(Return > threshold)",
    x = "First-Day Return Threshold",
    y = "Marginal Effect",
    caption = "Shaded area: 95% bootstrap CI"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

save_plot(p_dr_baseline, "01f_dist_reg_baseline.png")

cat("Running distribution regression (with controls)...\n")
dr_controls <- distribution_regression_controls(
  df_complete, "rf_fact_t07", "first_day_return",
  controls = c("log_assets", "log_proceeds", "roll_vwretd", "ff12"),
  n_thresholds = 25, n_bootstrap = 200
)

p_dr_controls <- ggplot(dr_controls, aes(x = threshold, y = marginal_effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = me_lower, ymax = me_upper), alpha = 0.2, fill = "darkgreen") +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 2) +
  labs(
    title = "Distribution Regression: With Controls",
    subtitle = "Controls: log(assets), log(proceeds), roll_vwretd, FF12 industry",
    x = "First-Day Return Threshold",
    y = "Marginal Effect",
    caption = "Shaded area: 95% bootstrap CI"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

save_plot(p_dr_controls, "01g_dist_reg_controls.png")

# --- Combined distribution regression plot ---
dr_baseline$specification <- "Baseline"
dr_controls$specification <- "With Controls"
dr_combined <- rbind(dr_baseline, dr_controls)

p_dr_combined <- ggplot(dr_combined, aes(x = threshold, y = marginal_effect,
                                          color = specification, fill = specification)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = me_lower, ymax = me_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Baseline" = "steelblue", "With Controls" = "darkgreen")) +
  scale_fill_manual(values = c("Baseline" = "steelblue", "With Controls" = "darkgreen")) +
  labs(
    title = "Distribution Regression: Baseline vs. With Controls",
    x = "First-Day Return Threshold",
    y = "Marginal Effect",
    color = "Specification",
    fill = "Specification"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

save_plot(p_dr_combined, "01h_dist_reg_combined.png")

cat("\n=== Script 01 Complete ===\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
