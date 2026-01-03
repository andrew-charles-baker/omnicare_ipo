# =============================================================================
# 03_cfm_decomposition.R
# CFM (Chernozhukov, Fernandez-Val, Melly 2013) Distributional Decomposition
# Three-way decomposition: Structure, rf_fact Composition, Other Composition
# =============================================================================

source(here::here("scripts", "00_setup.R"))

# --- Load Data ---
df <- load_ipo_data()

df_complete <- df[complete.cases(df[, c("first_day_return", "rf_fact_t07",
                                         "log_assets", "log_proceeds",
                                         "roll_vwretd", "ff12")]), ]

cat("Complete cases for decomposition:", nrow(df_complete), "\n")
cat("Pre-Omnicare:", sum(!df_complete$post_omnicare), "\n")
cat("Post-Omnicare:", sum(df_complete$post_omnicare), "\n")

# =============================================================================
# CFM DECOMPOSITION FUNCTIONS
# =============================================================================

# Estimate conditional distribution with multiple covariates
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
      logit_formula <- as.formula(formula_str)
      logit_fit <- glm(logit_formula, family = binomial(link = "logit"), data = temp_data)

      if (!logit_fit$converged) next

      distribution_results[[length(distribution_results) + 1]] <- list(
        y_threshold = y_threshold,
        fit = logit_fit
      )

      successful_fits <- c(successful_fits, y_threshold)

    }, error = function(e) {})
  }

  return(list(distribution_results = distribution_results, y_grid = successful_fits))
}

# Predict CDF given model and covariate matrix
predict_cdf_multi <- function(distribution_results, newdata, y_grid) {

  n_obs <- nrow(newdata)
  n_thresholds <- length(y_grid)
  cdf_matrix <- matrix(NA, nrow = n_obs, ncol = n_thresholds)

  actual_thresholds <- sapply(distribution_results, function(x) x$y_threshold)

  for (i in 1:n_thresholds) {
    match_idx <- which(abs(actual_thresholds - y_grid[i]) < 1e-10)

    if (length(match_idx) == 1) {
      fit <- distribution_results[[match_idx]]$fit
      pred <- predict(fit, newdata = newdata, type = "response")
      cdf_matrix[, i] <- pred
    }
  }

  return(cdf_matrix)
}

# Main CFM decomposition function
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

  if (length(common_y_grid) < 5) {
    stop("Insufficient common thresholds for decomposition.")
  }

  before_thresh <- sapply(cond_dist_before$distribution_results, function(x) x$y_threshold)
  after_thresh <- sapply(cond_dist_after$distribution_results, function(x) x$y_threshold)

  before_idx <- match(common_y_grid, before_thresh)
  after_idx <- match(common_y_grid, after_thresh)

  before_results <- cond_dist_before$distribution_results[before_idx]
  after_results <- cond_dist_after$distribution_results[after_idx]

  X_before <- before_data[, all_covariates, drop = FALSE]
  X_after <- after_data[, all_covariates, drop = FALSE]

  # Create hybrid: rf_fact from post, other covariates from pre
  set.seed(42)
  n_counterfactual <- min(nrow(X_before), nrow(X_after))
  before_sample_idx <- sample(1:nrow(X_before), n_counterfactual, replace = TRUE)
  after_sample_idx <- sample(1:nrow(X_after), n_counterfactual, replace = TRUE)

  X_hybrid <- X_before[before_sample_idx, , drop = FALSE]
  X_hybrid[[main_covariate]] <- X_after[after_sample_idx, main_covariate]

  # Four distributions
  cdf_00 <- predict_cdf_multi(before_results, X_before, common_y_grid)
  F_00 <- colMeans(cdf_00, na.rm = TRUE)

  cdf_11 <- predict_cdf_multi(after_results, X_after, common_y_grid)
  F_11 <- colMeans(cdf_11, na.rm = TRUE)

  cdf_01 <- predict_cdf_multi(before_results, X_after, common_y_grid)
  F_01 <- colMeans(cdf_01, na.rm = TRUE)

  cdf_0_hybrid <- predict_cdf_multi(before_results, X_hybrid, common_y_grid)
  F_0_hybrid <- colMeans(cdf_0_hybrid, na.rm = TRUE)

  # Decomposition
  total_change <- F_11 - F_00
  structure_effect <- F_11 - F_01
  other_composition <- F_01 - F_0_hybrid
  rf_composition <- F_0_hybrid - F_00

  results <- data.frame(
    y_threshold = common_y_grid,
    F_00 = F_00,
    F_11 = F_11,
    F_01 = F_01,
    F_0_hybrid = F_0_hybrid,
    total_change = total_change,
    structure_effect = structure_effect,
    other_composition = other_composition,
    rf_composition = rf_composition
  )

  return(list(
    decomposition = results,
    sample_sizes = list(before = nrow(before_data), after = nrow(after_data)),
    y_grid = common_y_grid
  ))
}

# Bootstrap version
bootstrap_cfm_full <- function(data, outcome_var, main_covariate, other_covariates,
                                group_var, y_grid = NULL, n_bootstrap = 200) {

  original <- cfm_full_decomposition(data, outcome_var, main_covariate,
                                      other_covariates, group_var, y_grid)
  y_grid_common <- original$y_grid

  boot_results <- list()

  cat("Running bootstrap (", n_bootstrap, "iterations)...\n")
  pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)

  for (b in 1:n_bootstrap) {
    setTxtProgressBar(pb, b)
    tryCatch({
      before_data <- data[data[[group_var]] == FALSE, ]
      after_data <- data[data[[group_var]] == TRUE, ]

      before_boot <- before_data[sample(nrow(before_data), replace = TRUE), ]
      after_boot <- after_data[sample(nrow(after_data), replace = TRUE), ]

      boot_data <- rbind(before_boot, after_boot)

      boot_decomp <- cfm_full_decomposition(boot_data, outcome_var, main_covariate,
                                             other_covariates, group_var, y_grid_common)
      boot_results[[b]] <- boot_decomp$decomposition

    }, error = function(e) {})
  }
  close(pb)

  boot_results <- boot_results[!sapply(boot_results, is.null)]
  cat("Successful bootstrap iterations:", length(boot_results), "\n")

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
    original$decomposition$other_ci_lower <- apply(other_boot, 2, quantile, 0.025, na.rm = TRUE)
    original$decomposition$other_ci_upper <- apply(other_boot, 2, quantile, 0.975, na.rm = TRUE)
    original$decomposition$rf_ci_lower <- apply(rf_boot, 2, quantile, 0.025, na.rm = TRUE)
    original$decomposition$rf_ci_upper <- apply(rf_boot, 2, quantile, 0.975, na.rm = TRUE)
  }

  return(original)
}


# =============================================================================
# RUN DECOMPOSITION
# =============================================================================

cat("\n=== Running Full CFM Decomposition ===\n")

main_cov <- "rf_fact_t07"
other_covs <- c("log_assets", "log_proceeds", "roll_vwretd")

full_decomp <- bootstrap_cfm_full(
  data = df_complete,
  outcome_var = "first_day_return",
  main_covariate = main_cov,
  other_covariates = other_covs,
  group_var = "post_omnicare",
  n_bootstrap = 300
)

decomp_df <- full_decomp$decomposition

# Summary
cat("\n=== DECOMPOSITION SUMMARY ===\n")
cat("Sample sizes: Pre =", full_decomp$sample_sizes$before,
    ", Post =", full_decomp$sample_sizes$after, "\n\n")
cat(sprintf("Average Total Change:        %7.4f\n", mean(decomp_df$total_change)))
cat(sprintf("Average Structure Effect:    %7.4f\n", mean(decomp_df$structure_effect)))
cat(sprintf("Average rf_fact Composition: %7.4f\n", mean(decomp_df$rf_composition)))
cat(sprintf("Average Other Composition:   %7.4f\n", mean(decomp_df$other_composition)))


# =============================================================================
# PLOT DECOMPOSITION
# =============================================================================

# Reshape for plotting
effect_long <- decomp_df %>%
  dplyr::select(y_threshold, total_change, structure_effect, rf_composition, other_composition) %>%
  tidyr::pivot_longer(cols = -y_threshold, names_to = "effect_type", values_to = "effect_value") %>%
  dplyr::mutate(
    effect_label = case_when(
      effect_type == "total_change" ~ "Total Change",
      effect_type == "structure_effect" ~ "Structure Effect (Pricing)",
      effect_type == "rf_composition" ~ "rf_fact Composition",
      effect_type == "other_composition" ~ "Other Covariates"
    ),
    effect_label = factor(effect_label, levels = c(
      "Total Change", "Structure Effect (Pricing)",
      "rf_fact Composition", "Other Covariates"
    ))
  )

has_ci <- "total_ci_lower" %in% names(decomp_df)

p_decomp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

if (has_ci) {
  p_decomp <- p_decomp +
    geom_ribbon(data = decomp_df,
                aes(x = y_threshold, ymin = total_ci_lower, ymax = total_ci_upper),
                fill = "black", alpha = 0.15) +
    geom_ribbon(data = decomp_df,
                aes(x = y_threshold, ymin = structure_ci_lower, ymax = structure_ci_upper),
                fill = "red", alpha = 0.15) +
    geom_ribbon(data = decomp_df,
                aes(x = y_threshold, ymin = rf_ci_lower, ymax = rf_ci_upper),
                fill = "blue", alpha = 0.15) +
    geom_ribbon(data = decomp_df,
                aes(x = y_threshold, ymin = other_ci_lower, ymax = other_ci_upper),
                fill = "darkgreen", alpha = 0.15)
}

p_decomp <- p_decomp +
  geom_line(data = effect_long, aes(x = y_threshold, y = effect_value, color = effect_label),
            linewidth = 1.2) +
  scale_color_manual(values = c(
    "Total Change" = "black",
    "Structure Effect (Pricing)" = "red",
    "rf_fact Composition" = "blue",
    "Other Covariates" = "darkgreen"
  )) +
  labs(
    title = "CFM Decomposition: Three-Way Attribution",
    subtitle = "Total = Structure + rf_fact Composition + Other Composition",
    x = "First-Day Return Threshold",
    y = "Change in Cumulative Probability",
    color = "Effect",
    caption = "Shaded areas: 95% bootstrap CI"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

save_plot(p_decomp, "03a_cfm_decomposition.png", width = 10, height = 7)


# =============================================================================
# OBSERVED VS COUNTERFACTUAL CDFs
# =============================================================================

cdf_long <- decomp_df %>%
  dplyr::select(y_threshold, F_00, F_11, F_01, F_0_hybrid) %>%
  tidyr::pivot_longer(cols = -y_threshold, names_to = "distribution", values_to = "cdf_value") %>%
  dplyr::mutate(
    dist_label = case_when(
      distribution == "F_00" ~ "Observed Pre-Omnicare",
      distribution == "F_11" ~ "Observed Post-Omnicare",
      distribution == "F_01" ~ "Pre Pricing + Post Covariates",
      distribution == "F_0_hybrid" ~ "Pre Pricing + Post rf_fact + Pre Other"
    ),
    dist_type = case_when(
      distribution %in% c("F_00", "F_11") ~ "Observed",
      TRUE ~ "Counterfactual"
    )
  )

p_cdfs <- ggplot(cdf_long, aes(x = y_threshold, y = cdf_value,
                                color = dist_label, linetype = dist_type)) +
  geom_line(linewidth = 1.1) +
  scale_linetype_manual(values = c("Observed" = "solid", "Counterfactual" = "dashed")) +
  scale_color_manual(values = c(
    "Observed Pre-Omnicare" = "coral",
    "Observed Post-Omnicare" = "steelblue",
    "Pre Pricing + Post Covariates" = "darkgreen",
    "Pre Pricing + Post rf_fact + Pre Other" = "purple"
  )) +
  labs(
    title = "CFM Decomposition: Observed and Counterfactual CDFs",
    x = "First-Day Return",
    y = "Cumulative Probability F(return)",
    color = "Distribution",
    linetype = "Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

save_plot(p_cdfs, "03b_cfm_counterfactual_cdfs.png", width = 10, height = 7)


# =============================================================================
# RELATIVE CONTRIBUTIONS BAR CHART
# =============================================================================

avg_effects <- data.frame(
  effect = c("Structure Effect", "rf_fact Composition", "Other Composition"),
  value = c(
    mean(decomp_df$structure_effect),
    mean(decomp_df$rf_composition),
    mean(decomp_df$other_composition)
  )
)
avg_effects$pct <- 100 * abs(avg_effects$value) / sum(abs(avg_effects$value))
avg_effects$effect <- factor(avg_effects$effect, levels = avg_effects$effect)

p_bar <- ggplot(avg_effects, aes(x = effect, y = pct, fill = effect)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Structure Effect" = "red",
                                "rf_fact Composition" = "blue",
                                "Other Composition" = "darkgreen")) +
  labs(
    title = "Relative Contributions to Total Change",
    subtitle = "Based on average absolute effect across distribution",
    x = "",
    y = "Percentage of Total"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  ylim(0, max(avg_effects$pct) * 1.15)

save_plot(p_bar, "03c_cfm_relative_contributions.png", width = 8, height = 6)

# Save decomposition data
write.csv(decomp_df, file.path(OUTPUT_DIR, "03_cfm_decomposition_data.csv"), row.names = FALSE)

cat("\n=== Script 03 Complete ===\n")
cat("Output files saved to:", OUTPUT_DIR, "\n")
