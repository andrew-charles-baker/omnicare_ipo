# =============================================================================
# 04_circuit_heterogeneity.R
# Circuit-level heterogeneity analysis
# 6th Circuit (reduced liability) vs 2nd/9th Circuit (expanded liability)
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# Source CFM functions from script 03
# (we'll redefine key functions here for standalone use)

# --- Load Data ---
df <- load_ipo_data()

df_complete <- df[complete.cases(df[, c("first_day_return", "rf_fact_t07",
                                         "log_assets", "log_proceeds",
                                         "roll_vwretd", "ff12")]), ]

# Add circuit classification
df_complete$circuit <- assign_circuit(df_complete$headquarter)

cat("=== CIRCUIT DISTRIBUTION ===\n")
print(table(df_complete$circuit, df_complete$post_omnicare))


# =============================================================================
# CFM FUNCTIONS (copied from 03 for standalone use)
# =============================================================================

estimate_cond_dist_multi <- function(data, outcome_var, covariates, y_grid = NULL) {
  if (is.null(y_grid)) {
    y_values <- data[[outcome_var]]
    y_grid <- quantile(y_values, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE)
    y_grid <- unique(y_grid)
  }

  distribution_results <- list()
  successful_fits <- c()

  for (j in seq_along(y_grid)) {
    y_threshold <- y_grid[j]
    temp_data <- data
    temp_data$binary_indicator <- as.numeric(data[[outcome_var]] <= y_threshold)

    prop_below <- mean(temp_data$binary_indicator)
    if (prop_below < 0.02 || prop_below > 0.98) next

    tryCatch({
      formula_str <- paste("binary_indicator ~", paste(covariates, collapse = " + "))
      logit_fit <- glm(as.formula(formula_str), family = binomial(link = "logit"), data = temp_data)
      if (!logit_fit$converged) next
      distribution_results[[length(distribution_results) + 1]] <- list(y_threshold = y_threshold, fit = logit_fit)
      successful_fits <- c(successful_fits, y_threshold)
    }, error = function(e) {})
  }
  return(list(distribution_results = distribution_results, y_grid = successful_fits))
}

predict_cdf_multi <- function(distribution_results, newdata, y_grid) {
  n_obs <- nrow(newdata)
  n_thresholds <- length(y_grid)
  cdf_matrix <- matrix(NA, nrow = n_obs, ncol = n_thresholds)
  actual_thresholds <- sapply(distribution_results, function(x) x$y_threshold)
  for (i in 1:n_thresholds) {
    match_idx <- which(abs(actual_thresholds - y_grid[i]) < 1e-10)
    if (length(match_idx) == 1) {
      fit <- distribution_results[[match_idx]]$fit
      cdf_matrix[, i] <- predict(fit, newdata = newdata, type = "response")
    }
  }
  return(cdf_matrix)
}

cfm_full_decomposition <- function(data, outcome_var, main_covariate, other_covariates,
                                    group_var, y_grid = NULL) {
  all_covariates <- c(main_covariate, other_covariates)
  before_data <- data[data[[group_var]] == FALSE, ]
  after_data <- data[data[[group_var]] == TRUE, ]
  before_data <- before_data[complete.cases(before_data[, c(outcome_var, all_covariates)]), ]
  after_data <- after_data[complete.cases(after_data[, c(outcome_var, all_covariates)]), ]

  if (is.null(y_grid)) {
    all_outcome <- c(before_data[[outcome_var]], after_data[[outcome_var]])
    y_grid <- quantile(all_outcome, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE)
    y_grid <- unique(y_grid)
  }

  cond_dist_before <- estimate_cond_dist_multi(before_data, outcome_var, all_covariates, y_grid)
  cond_dist_after <- estimate_cond_dist_multi(after_data, outcome_var, all_covariates, y_grid)
  common_y_grid <- intersect(cond_dist_before$y_grid, cond_dist_after$y_grid)

  if (length(common_y_grid) < 5) stop("Insufficient common thresholds.")

  before_thresh <- sapply(cond_dist_before$distribution_results, function(x) x$y_threshold)
  after_thresh <- sapply(cond_dist_after$distribution_results, function(x) x$y_threshold)
  before_idx <- match(common_y_grid, before_thresh)
  after_idx <- match(common_y_grid, after_thresh)
  before_results <- cond_dist_before$distribution_results[before_idx]
  after_results <- cond_dist_after$distribution_results[after_idx]

  X_before <- before_data[, all_covariates, drop = FALSE]
  X_after <- after_data[, all_covariates, drop = FALSE]

  set.seed(42)
  n_counterfactual <- min(nrow(X_before), nrow(X_after))
  before_sample_idx <- sample(1:nrow(X_before), n_counterfactual, replace = TRUE)
  after_sample_idx <- sample(1:nrow(X_after), n_counterfactual, replace = TRUE)
  X_hybrid <- X_before[before_sample_idx, , drop = FALSE]
  X_hybrid[[main_covariate]] <- X_after[after_sample_idx, main_covariate]

  F_00 <- colMeans(predict_cdf_multi(before_results, X_before, common_y_grid), na.rm = TRUE)
  F_11 <- colMeans(predict_cdf_multi(after_results, X_after, common_y_grid), na.rm = TRUE)
  F_01 <- colMeans(predict_cdf_multi(before_results, X_after, common_y_grid), na.rm = TRUE)
  F_0_hybrid <- colMeans(predict_cdf_multi(before_results, X_hybrid, common_y_grid), na.rm = TRUE)

  results <- data.frame(
    y_threshold = common_y_grid,
    total_change = F_11 - F_00,
    structure_effect = F_11 - F_01,
    other_composition = F_01 - F_0_hybrid,
    rf_composition = F_0_hybrid - F_00
  )

  return(list(decomposition = results, sample_sizes = list(before = nrow(before_data), after = nrow(after_data)), y_grid = common_y_grid))
}

bootstrap_cfm_full <- function(data, outcome_var, main_covariate, other_covariates,
                                group_var, y_grid = NULL, n_bootstrap = 200) {
  original <- cfm_full_decomposition(data, outcome_var, main_covariate, other_covariates, group_var, y_grid)
  y_grid_common <- original$y_grid
  boot_results <- list()

  for (b in 1:n_bootstrap) {
    tryCatch({
      before_data <- data[data[[group_var]] == FALSE, ]
      after_data <- data[data[[group_var]] == TRUE, ]
      before_boot <- before_data[sample(nrow(before_data), replace = TRUE), ]
      after_boot <- after_data[sample(nrow(after_data), replace = TRUE), ]
      boot_data <- rbind(before_boot, after_boot)
      boot_decomp <- cfm_full_decomposition(boot_data, outcome_var, main_covariate, other_covariates, group_var, y_grid_common)
      boot_results[[b]] <- boot_decomp$decomposition
    }, error = function(e) {})
  }

  boot_results <- boot_results[!sapply(boot_results, is.null)]

  if (length(boot_results) > 10) {
    n_success <- length(boot_results)
    n_thresh <- length(y_grid_common)
    total_boot <- structure_boot <- other_boot <- rf_boot <- matrix(NA, n_success, n_thresh)
    for (i in 1:n_success) {
      if (nrow(boot_results[[i]]) == n_thresh) {
        total_boot[i, ] <- boot_results[[i]]$total_change
        structure_boot[i, ] <- boot_results[[i]]$structure_effect
        other_boot[i, ] <- boot_results[[i]]$other_composition
        rf_boot[i, ] <- boot_results[[i]]$rf_composition
      }
    }
    original$decomposition$total_ci_lower <- apply(total_boot, 2, quantile, 0.025, na.rm = TRUE)
    original$decomposition$total_ci_upper <- apply(total_boot, 2, quantile, 0.975, na.rm = TRUE)
    original$decomposition$structure_ci_lower <- apply(structure_boot, 2, quantile, 0.025, na.rm = TRUE)
    original$decomposition$structure_ci_upper <- apply(structure_boot, 2, quantile, 0.975, na.rm = TRUE)
    original$decomposition$rf_ci_lower <- apply(rf_boot, 2, quantile, 0.025, na.rm = TRUE)
    original$decomposition$rf_ci_upper <- apply(rf_boot, 2, quantile, 0.975, na.rm = TRUE)
    original$decomposition$other_ci_lower <- apply(other_boot, 2, quantile, 0.025, na.rm = TRUE)
    original$decomposition$other_ci_upper <- apply(other_boot, 2, quantile, 0.975, na.rm = TRUE)
  }

  return(original)
}


# =============================================================================
# 1. DISTRIBUTION REGRESSION BY CIRCUIT AND PERIOD
# =============================================================================

# Subset by circuit
df_6th <- df_complete[df_complete$circuit == "6th Circuit", ]
df_2nd9th <- df_complete[df_complete$circuit == "2nd/9th Circuit", ]

controls <- c("log_assets", "log_proceeds", "roll_vwretd", "ff12")

# --- 6th Circuit ---
cat("\n=== 6TH CIRCUIT ANALYSIS ===\n")
cat("Pre-Omnicare:", sum(!df_6th$post_omnicare), "\n")
cat("Post-Omnicare:", sum(df_6th$post_omnicare), "\n")

if (sum(!df_6th$post_omnicare) >= 20 && sum(df_6th$post_omnicare) >= 20) {
  dr_6th_pre <- distribution_regression_controls(
    df_6th[!df_6th$post_omnicare, ], "rf_fact_t07", "first_day_return", controls,
    n_thresholds = 15, n_bootstrap = 150
  )
  dr_6th_pre$period <- "Pre-Omnicare"
  dr_6th_pre$circuit <- "6th Circuit"

  dr_6th_post <- distribution_regression_controls(
    df_6th[df_6th$post_omnicare, ], "rf_fact_t07", "first_day_return", controls,
    n_thresholds = 15, n_bootstrap = 150
  )
  dr_6th_post$period <- "Post-Omnicare"
  dr_6th_post$circuit <- "6th Circuit"

  dr_6th <- rbind(dr_6th_pre, dr_6th_post)
} else {
  cat("Insufficient data for 6th Circuit analysis\n")
  dr_6th <- NULL
}

# --- 2nd/9th Circuit ---
cat("\n=== 2ND/9TH CIRCUIT ANALYSIS ===\n")
cat("Pre-Omnicare:", sum(!df_2nd9th$post_omnicare), "\n")
cat("Post-Omnicare:", sum(df_2nd9th$post_omnicare), "\n")

if (sum(!df_2nd9th$post_omnicare) >= 20 && sum(df_2nd9th$post_omnicare) >= 20) {
  dr_2nd9th_pre <- distribution_regression_controls(
    df_2nd9th[!df_2nd9th$post_omnicare, ], "rf_fact_t07", "first_day_return", controls,
    n_thresholds = 20, n_bootstrap = 150
  )
  dr_2nd9th_pre$period <- "Pre-Omnicare"
  dr_2nd9th_pre$circuit <- "2nd/9th Circuit"

  dr_2nd9th_post <- distribution_regression_controls(
    df_2nd9th[df_2nd9th$post_omnicare, ], "rf_fact_t07", "first_day_return", controls,
    n_thresholds = 20, n_bootstrap = 150
  )
  dr_2nd9th_post$period <- "Post-Omnicare"
  dr_2nd9th_post$circuit <- "2nd/9th Circuit"

  dr_2nd9th <- rbind(dr_2nd9th_pre, dr_2nd9th_post)
} else {
  cat("Insufficient data for 2nd/9th Circuit analysis\n")
  dr_2nd9th <- NULL
}


# =============================================================================
# 2. PLOT DISTRIBUTION REGRESSION BY CIRCUIT
# =============================================================================

plot_dr_circuit <- function(dr_data, circuit_name, subtitle_text) {
  if (is.null(dr_data)) return(NULL)

  ggplot(dr_data, aes(x = threshold, y = marginal_effect, color = period, fill = period)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = me_lower, ymax = me_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
    scale_fill_manual(values = c("Pre-Omnicare" = "coral", "Post-Omnicare" = "steelblue")) +
    labs(
      title = circuit_name,
      subtitle = subtitle_text,
      x = "First-Day Return Threshold",
      y = "Marginal Effect",
      color = "Period",
      fill = "Period"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
}

p_dr_6th <- plot_dr_circuit(dr_6th, "6th Circuit (KY, MI, OH, TN)",
                             "Omnicare REDUCED liability")
p_dr_2nd9th <- plot_dr_circuit(dr_2nd9th, "2nd/9th Circuit (NY, CA, etc.)",
                                "Omnicare EXPANDED liability")

if (!is.null(p_dr_6th) && !is.null(p_dr_2nd9th)) {
  p_dr_circuits <- p_dr_6th + p_dr_2nd9th +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  save_plot(p_dr_circuits, "04a_dist_reg_by_circuit.png", width = 14, height = 6)
}


# =============================================================================
# 3. CFM DECOMPOSITION BY CIRCUIT
# =============================================================================

main_cov <- "rf_fact_t07"
other_covs <- c("log_assets", "log_proceeds", "roll_vwretd")

# --- 6th Circuit Decomposition ---
cat("\n=== Running CFM Decomposition for 6th Circuit ===\n")
if (sum(!df_6th$post_omnicare) >= 20 && sum(df_6th$post_omnicare) >= 20) {
  decomp_6th <- tryCatch({
    bootstrap_cfm_full(df_6th, "first_day_return", main_cov, other_covs, "post_omnicare", n_bootstrap = 200)
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    NULL
  })
} else {
  decomp_6th <- NULL
}

# --- 2nd/9th Circuit Decomposition ---
cat("\n=== Running CFM Decomposition for 2nd/9th Circuit ===\n")
if (sum(!df_2nd9th$post_omnicare) >= 20 && sum(df_2nd9th$post_omnicare) >= 20) {
  decomp_2nd9th <- tryCatch({
    bootstrap_cfm_full(df_2nd9th, "first_day_return", main_cov, other_covs, "post_omnicare", n_bootstrap = 200)
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    NULL
  })
} else {
  decomp_2nd9th <- NULL
}


# =============================================================================
# 4. PLOT DECOMPOSITION BY CIRCUIT
# =============================================================================

plot_decomp_circuit <- function(decomp_result, circuit_name, subtitle_text) {
  if (is.null(decomp_result)) return(NULL)

  decomp_df <- decomp_result$decomposition

  effect_long <- decomp_df %>%
    dplyr::select(y_threshold, total_change, structure_effect, rf_composition, other_composition) %>%
    tidyr::pivot_longer(cols = -y_threshold, names_to = "effect_type", values_to = "effect_value") %>%
    dplyr::mutate(
      effect_label = case_when(
        effect_type == "total_change" ~ "Total Change",
        effect_type == "structure_effect" ~ "Structure Effect",
        effect_type == "rf_composition" ~ "rf_fact Composition",
        effect_type == "other_composition" ~ "Other Composition"
      )
    )

  has_ci <- "total_ci_lower" %in% names(decomp_df)

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

  if (has_ci) {
    p <- p +
      geom_ribbon(data = decomp_df, aes(x = y_threshold, ymin = total_ci_lower, ymax = total_ci_upper),
                  fill = "black", alpha = 0.15) +
      geom_ribbon(data = decomp_df, aes(x = y_threshold, ymin = structure_ci_lower, ymax = structure_ci_upper),
                  fill = "red", alpha = 0.15) +
      geom_ribbon(data = decomp_df, aes(x = y_threshold, ymin = rf_ci_lower, ymax = rf_ci_upper),
                  fill = "blue", alpha = 0.15) +
      geom_ribbon(data = decomp_df, aes(x = y_threshold, ymin = other_ci_lower, ymax = other_ci_upper),
                  fill = "darkgreen", alpha = 0.15)
  }

  p <- p +
    geom_line(data = effect_long, aes(x = y_threshold, y = effect_value, color = effect_label), linewidth = 1.1) +
    scale_color_manual(values = c(
      "Total Change" = "black", "Structure Effect" = "red",
      "rf_fact Composition" = "blue", "Other Composition" = "darkgreen"
    )) +
    labs(
      title = circuit_name,
      subtitle = subtitle_text,
      x = "First-Day Return Threshold",
      y = "Change in Cumulative Probability",
      color = "Effect"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  return(p)
}

p_decomp_6th <- plot_decomp_circuit(decomp_6th, "6th Circuit (KY, MI, OH, TN)",
                                     "Omnicare REDUCED liability")
p_decomp_2nd9th <- plot_decomp_circuit(decomp_2nd9th, "2nd/9th Circuit (NY, CA, etc.)",
                                        "Omnicare EXPANDED liability")

if (!is.null(p_decomp_6th) && !is.null(p_decomp_2nd9th)) {
  p_decomp_circuits <- p_decomp_6th / p_decomp_2nd9th +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  save_plot(p_decomp_circuits, "04b_cfm_decomposition_by_circuit.png", width = 10, height = 12)
}


# =============================================================================
# 5. SUMMARY COMPARISON TABLE
# =============================================================================

cat("\n=== CIRCUIT COMPARISON SUMMARY ===\n\n")

if (!is.null(decomp_6th)) {
  cat("6TH CIRCUIT (Omnicare reduced liability):\n")
  cat("  Sample: Pre =", decomp_6th$sample_sizes$before, ", Post =", decomp_6th$sample_sizes$after, "\n")
  cat(sprintf("  Avg Total Change:        %7.4f\n", mean(decomp_6th$decomposition$total_change)))
  cat(sprintf("  Avg Structure Effect:    %7.4f\n", mean(decomp_6th$decomposition$structure_effect)))
  cat(sprintf("  Avg rf_fact Composition: %7.4f\n", mean(decomp_6th$decomposition$rf_composition)))
  cat(sprintf("  Avg Other Composition:   %7.4f\n\n", mean(decomp_6th$decomposition$other_composition)))
}

if (!is.null(decomp_2nd9th)) {
  cat("2ND/9TH CIRCUIT (Omnicare expanded liability):\n")
  cat("  Sample: Pre =", decomp_2nd9th$sample_sizes$before, ", Post =", decomp_2nd9th$sample_sizes$after, "\n")
  cat(sprintf("  Avg Total Change:        %7.4f\n", mean(decomp_2nd9th$decomposition$total_change)))
  cat(sprintf("  Avg Structure Effect:    %7.4f\n", mean(decomp_2nd9th$decomposition$structure_effect)))
  cat(sprintf("  Avg rf_fact Composition: %7.4f\n", mean(decomp_2nd9th$decomposition$rf_composition)))
  cat(sprintf("  Avg Other Composition:   %7.4f\n", mean(decomp_2nd9th$decomposition$other_composition)))
}


# =============================================================================
# 6. RF_FACT SHIFT BY CIRCUIT
# =============================================================================

rf_shift <- df_complete %>%
  group_by(circuit, post_omnicare) %>%
  summarise(
    n = n(),
    mean_rf_fact = mean(rf_fact_t07, na.rm = TRUE),
    mean_return = mean(first_day_return, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(period = ifelse(post_omnicare, "Post", "Pre"))

cat("\n=== RF_FACT SHIFT BY CIRCUIT ===\n")
print(rf_shift %>% filter(circuit != "Other") %>% arrange(circuit, period))

# Save summary data
write.csv(rf_shift, file.path(OUTPUT_DIR, "04_circuit_summary_stats.csv"), row.names = FALSE)

cat("\n=== Script 04 Complete ===\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
