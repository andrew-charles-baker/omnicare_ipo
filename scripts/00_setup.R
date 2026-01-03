# =============================================================================
# 00_setup.R
# Setup script: Load packages, data, and define common functions
# =============================================================================

# --- Load Packages ---
library(tidyverse)
library(ggplot2)
library(patchwork)

# Try to load binsreg if available
binsreg_available <- requireNamespace("binsreg", quietly = TRUE)
if (binsreg_available) {
  library(binsreg)
  cat("binsreg package loaded successfully\n")
} else {
  cat("binsreg not available - using fallback binning\n")
}

# --- Set Theme ---
theme_set(theme_minimal(base_size = 12))

# Output directory
OUTPUT_DIR <- here::here("output", "figures")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# --- Load and Prepare Data ---
load_ipo_data <- function() {
  df <- read.csv(here::here("data", "omnicare_rf.csv"))

  # Filter to offer_price >= $5

df <- df[df$offer_price >= 5, ]

  # Create derived variables
  df$offer_date <- as.Date(df$offer_date)
  df$post_omnicare <- df$offer_date >= as.Date("2015-03-24")
  df$log_assets <- log(df$assets_thous + 1)
  df$log_proceeds <- log(df$total_proceeds_thous + 1)

  # Fama-French 12 Industry Classification
  df$ff12 <- factor(assign_ff12(df$sic_code))

  return(df)
}

# --- Fama-French 12 Industry Classification ---
assign_ff12 <- function(sic) {
  ff12 <- rep(NA_character_, length(sic))

  # 1. Consumer NonDurables
  ff12[sic >= 100 & sic <= 999] <- "NoDur"
  ff12[sic >= 2000 & sic <= 2399] <- "NoDur"
  ff12[sic >= 2700 & sic <= 2749] <- "NoDur"
  ff12[sic >= 2770 & sic <= 2799] <- "NoDur"
  ff12[sic >= 3100 & sic <= 3199] <- "NoDur"
  ff12[sic >= 3940 & sic <= 3989] <- "NoDur"

 # 2. Consumer Durables
  ff12[sic >= 2500 & sic <= 2519] <- "Durbl"
  ff12[sic >= 2590 & sic <= 2599] <- "Durbl"
  ff12[sic >= 3630 & sic <= 3659] <- "Durbl"
  ff12[sic >= 3710 & sic <= 3711] <- "Durbl"
  ff12[sic >= 3714 & sic <= 3714] <- "Durbl"
  ff12[sic >= 3716 & sic <= 3716] <- "Durbl"
  ff12[sic >= 3750 & sic <= 3751] <- "Durbl"
  ff12[sic >= 3792 & sic <= 3792] <- "Durbl"
  ff12[sic >= 3900 & sic <= 3939] <- "Durbl"
  ff12[sic >= 3990 & sic <= 3999] <- "Durbl"

  # 3. Manufacturing
  ff12[sic >= 2520 & sic <= 2589] <- "Manuf"
  ff12[sic >= 2600 & sic <= 2699] <- "Manuf"
  ff12[sic >= 2750 & sic <= 2769] <- "Manuf"
  ff12[sic >= 2800 & sic <= 2829] <- "Manuf"
  ff12[sic >= 2840 & sic <= 2899] <- "Manuf"
  ff12[sic >= 3000 & sic <= 3099] <- "Manuf"
  ff12[sic >= 3200 & sic <= 3569] <- "Manuf"
  ff12[sic >= 3580 & sic <= 3621] <- "Manuf"
  ff12[sic >= 3623 & sic <= 3629] <- "Manuf"
  ff12[sic >= 3700 & sic <= 3709] <- "Manuf"
  ff12[sic >= 3712 & sic <= 3713] <- "Manuf"
  ff12[sic >= 3715 & sic <= 3715] <- "Manuf"
  ff12[sic >= 3717 & sic <= 3749] <- "Manuf"
  ff12[sic >= 3752 & sic <= 3791] <- "Manuf"
  ff12[sic >= 3793 & sic <= 3799] <- "Manuf"
  ff12[sic >= 3830 & sic <= 3839] <- "Manuf"
  ff12[sic >= 3860 & sic <= 3899] <- "Manuf"

  # 4. Energy
  ff12[sic >= 1200 & sic <= 1399] <- "Enrgy"
  ff12[sic >= 2900 & sic <= 2999] <- "Enrgy"

  # 5. Chemicals
  ff12[sic >= 2830 & sic <= 2839] <- "Chems"

  # 6. Business Equipment
  ff12[sic >= 3570 & sic <= 3579] <- "BusEq"
  ff12[sic >= 3622 & sic <= 3622] <- "BusEq"
  ff12[sic >= 3660 & sic <= 3692] <- "BusEq"
  ff12[sic >= 3694 & sic <= 3699] <- "BusEq"
  ff12[sic >= 3810 & sic <= 3829] <- "BusEq"
  ff12[sic >= 7370 & sic <= 7379] <- "BusEq"

  # 7. Telecom
  ff12[sic >= 4800 & sic <= 4899] <- "Telcm"

  # 8. Utilities
  ff12[sic >= 4900 & sic <= 4949] <- "Utils"

  # 9. Shops (Wholesale, Retail)
  ff12[sic >= 5000 & sic <= 5999] <- "Shops"
  ff12[sic >= 7200 & sic <= 7299] <- "Shops"
  ff12[sic >= 7600 & sic <= 7699] <- "Shops"

  # 10. Healthcare
  ff12[sic >= 2830 & sic <= 2839] <- "Hlth"
  ff12[sic >= 3693 & sic <= 3693] <- "Hlth"
  ff12[sic >= 3840 & sic <= 3859] <- "Hlth"
  ff12[sic >= 8000 & sic <= 8099] <- "Hlth"

  # 11. Finance
  ff12[sic >= 6000 & sic <= 6999] <- "Money"

  # 12. Other
  ff12[is.na(ff12)] <- "Other"

  return(ff12)
}

# --- Circuit Classification ---
assign_circuit <- function(state) {
  sixth_circuit <- c("KY", "MI", "OH", "TN")
  second_circuit <- c("CT", "NY", "VT")
  ninth_circuit <- c("AK", "AZ", "CA", "HI", "ID", "MT", "NV", "OR", "WA")

  case_when(
    state %in% sixth_circuit ~ "6th Circuit",
    state %in% c(second_circuit, ninth_circuit) ~ "2nd/9th Circuit",
    TRUE ~ "Other"
  )
}

# --- Distribution Regression Function ---
distribution_regression <- function(data, x_var, y_var, n_thresholds = 25, n_bootstrap = 200) {

  y_values <- data[[y_var]]
  thresholds <- quantile(y_values, probs = seq(0.05, 0.95, length.out = n_thresholds), na.rm = TRUE)
  thresholds <- unique(thresholds)

  results <- data.frame()

  for (tau in thresholds) {
    data$y_binary <- as.numeric(data[[y_var]] > tau)

    prop_above <- mean(data$y_binary, na.rm = TRUE)
    if (prop_above < 0.02 || prop_above > 0.98) next

    tryCatch({
      formula <- as.formula(paste("y_binary ~", x_var))
      model <- glm(formula, family = binomial(link = "logit"), data = data)

      if (!model$converged) next

      beta <- coef(model)[2]
      pred <- predict(model, type = "response")
      avg_deriv <- mean(pred * (1 - pred), na.rm = TRUE)
      marginal_effect <- beta * avg_deriv

      # Bootstrap CI
      boot_me <- numeric(n_bootstrap)
      for (b in 1:n_bootstrap) {
        tryCatch({
          boot_idx <- sample(nrow(data), replace = TRUE)
          boot_data <- data[boot_idx, ]
          boot_model <- glm(formula, family = binomial(link = "logit"), data = boot_data)
          if (boot_model$converged) {
            boot_beta <- coef(boot_model)[2]
            boot_pred <- predict(boot_model, type = "response")
            boot_me[b] <- boot_beta * mean(boot_pred * (1 - boot_pred), na.rm = TRUE)
          }
        }, error = function(e) {})
      }

      boot_me <- boot_me[boot_me != 0 & !is.na(boot_me)]
      me_lower <- if (length(boot_me) > 10) quantile(boot_me, 0.025) else NA
      me_upper <- if (length(boot_me) > 10) quantile(boot_me, 0.975) else NA

      results <- rbind(results, data.frame(
        threshold = tau,
        marginal_effect = marginal_effect,
        me_lower = me_lower,
        me_upper = me_upper
      ))

    }, error = function(e) {})
  }

  return(results)
}

# --- Distribution Regression with Controls ---
distribution_regression_controls <- function(data, x_var, y_var, controls,
                                              n_thresholds = 25, n_bootstrap = 200) {

  y_values <- data[[y_var]]
  thresholds <- quantile(y_values, probs = seq(0.05, 0.95, length.out = n_thresholds), na.rm = TRUE)
  thresholds <- unique(thresholds)

  results <- data.frame()

  for (tau in thresholds) {
    data$y_binary <- as.numeric(data[[y_var]] > tau)

    prop_above <- mean(data$y_binary, na.rm = TRUE)
    if (prop_above < 0.02 || prop_above > 0.98) next

    tryCatch({
      formula_str <- paste("y_binary ~", x_var, "+", paste(controls, collapse = " + "))
      formula <- as.formula(formula_str)
      model <- glm(formula, family = binomial(link = "logit"), data = data)

      if (!model$converged) next

      beta <- coef(model)[x_var]
      pred <- predict(model, type = "response")
      avg_deriv <- mean(pred * (1 - pred), na.rm = TRUE)
      marginal_effect <- beta * avg_deriv

      # Bootstrap CI
      boot_me <- numeric(n_bootstrap)
      for (b in 1:n_bootstrap) {
        tryCatch({
          boot_idx <- sample(nrow(data), replace = TRUE)
          boot_data <- data[boot_idx, ]
          boot_model <- glm(formula, family = binomial(link = "logit"), data = boot_data)
          if (boot_model$converged) {
            boot_beta <- coef(boot_model)[x_var]
            boot_pred <- predict(boot_model, type = "response")
            boot_me[b] <- boot_beta * mean(boot_pred * (1 - boot_pred), na.rm = TRUE)
          }
        }, error = function(e) {})
      }

      boot_me <- boot_me[boot_me != 0 & !is.na(boot_me)]
      me_lower <- if (length(boot_me) > 10) quantile(boot_me, 0.025) else NA
      me_upper <- if (length(boot_me) > 10) quantile(boot_me, 0.975) else NA

      results <- rbind(results, data.frame(
        threshold = tau,
        marginal_effect = marginal_effect,
        me_lower = me_lower,
        me_upper = me_upper
      ))

    }, error = function(e) {})
  }

  return(results)
}

# --- Save Plot Function ---
save_plot <- function(plot, filename, width = 8, height = 6) {
  filepath <- file.path(OUTPUT_DIR, filename)
  ggsave(filepath, plot, width = width, height = height, dpi = 300)
  cat("Saved:", filepath, "\n")
}

cat("Setup complete. Use load_ipo_data() to load data.\n")
