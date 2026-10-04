source("config.R")

# ============================================================
# 23_standard_product_integral_point_estimates.R
#
# Point-estimate check using the standard product-integral /
# Aalen-Johansen update for each origin-specific competing-risk
# process, keeping fixed:
#   - the Cox models fitted in file 05
#   - Efron baseline cumulative-hazard increments
#   - time-varying implant-order effects
#   - common origin-specific standardization distributions
#
# No bootstrap is run here.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(openxlsx)
})

out_dir <- output_dir

input_file <- file.path(
  out_dir,
  "p1_multistate_state_probabilities.RData"
)

if (!file.exists(input_file)) {
  stop(
    "Cannot find: ", input_file,
    "\nRun 05_p1_multistate_state_probabilities.R first."
  )
}

# Load into a dedicated environment so the script is independent
# of the current Global Environment.
e <- new.env(parent = emptyenv())
loaded_objects <- load(input_file, envir = e)

required_objects <- c(
  "transition_map",
  "prob_model_objects",
  "baseline_increments",
  "standardization_patterns",
  "adjusted_probabilities_report"
)

missing_objects <- setdiff(required_objects, loaded_objects)

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects in ", input_file, ": ",
    paste(missing_objects, collapse = ", ")
  )
}

transition_map <- e$transition_map
prob_model_objects <- e$prob_model_objects
baseline_increments <- e$baseline_increments
standardization_patterns <- e$standardization_patterns
pe_report <- e$adjusted_probabilities_report

times_report <- c(1, 3, 5, 10, 15)
tol <- 1e-10

# ------------------------------------------------------------
# Linear predictor, matching file 05
# ------------------------------------------------------------

calc_lp_times_pi <- function(pattern, beta, times) {
  if (!is.data.frame(pattern)) {
    pattern <- as.data.frame(pattern)
  }
  
  n_pattern <- nrow(pattern)
  n_times <- length(times)
  n <- max(n_pattern, n_times)
  
  if (n_pattern == 1 && n > 1) {
    pattern <- pattern[rep(1, n), , drop = FALSE]
  } else if (n_pattern != n) {
    stop("Incompatible dimensions between pattern and times.")
  }
  
  if (n_times == 1 && n > 1) {
    times <- rep(times, n)
  } else if (n_times != n) {
    stop("Incompatible dimensions between pattern and times.")
  }
  
  hip_value <- as.numeric(pattern$hip_binary)
  sex_value <- as.character(pattern$sex_ms)
  bmi_value <- as.character(pattern$bmi_p1_f_ms)
  
  lp <- rep(0, n)
  
  if ("hipfirst_hip" %in% names(beta)) {
    lp <- lp + beta["hipfirst_hip"] * hip_value
  }
  
  if ("hip_binary" %in% names(beta)) {
    lp <- lp + beta["hip_binary"] * hip_value
  }
  
  if ("tt(hip_binary)" %in% names(beta)) {
    lp <- lp +
      beta["tt(hip_binary)"] *
      hip_value *
      log(times + 1)
  }
  
  if ("sex_msfemale" %in% names(beta)) {
    lp <- lp +
      beta["sex_msfemale"] *
      as.numeric(sex_value == "female")
  }
  
  if ("bmi_p1_f_msoverweight" %in% names(beta)) {
    lp <- lp +
      beta["bmi_p1_f_msoverweight"] *
      as.numeric(bmi_value == "overweight")
  }
  
  if ("bmi_p1_f_msobese" %in% names(beta)) {
    lp <- lp +
      beta["bmi_p1_f_msobese"] *
      as.numeric(bmi_value == "obese")
  }
  
  if ("bmi_p1_f_msmissing" %in% names(beta)) {
    lp <- lp +
      beta["bmi_p1_f_msmissing"] *
      as.numeric(bmi_value == "missing")
  }
  
  as.numeric(lp)
}

make_pattern_pi <- function(age_stratum, sex_value, bmi_value, hip_value) {
  tibble(
    age_cat = factor(
      age_stratum,
      levels = c("<60", "60-69", "70-79", "80+")
    ),
    sex_ms = factor(
      sex_value,
      levels = c("male", "female")
    ),
    bmi_p1_f_ms = factor(
      bmi_value,
      levels = c(
        "normal_underweight",
        "overweight",
        "obese",
        "missing"
      )
    ),
    hip_binary = as.numeric(hip_value)
  )
}

# ------------------------------------------------------------
# Predicted cause-specific cumulative-hazard increments
# for one covariate pattern
# ------------------------------------------------------------

make_hazard_increments_pi <- function(
    origin_name,
    age_stratum,
    sex_value,
    bmi_value,
    hip_value
) {
  meta_origin <- transition_map %>%
    filter(as.character(origin) == origin_name)
  
  pattern <- make_pattern_pi(
    age_stratum = age_stratum,
    sex_value = sex_value,
    bmi_value = bmi_value,
    hip_value = hip_value
  )
  
  hazard_parts <- map_dfr(
    as.character(meta_origin$trans),
    function(tr) {
      inc <- baseline_increments %>%
        filter(
          .data$trans == .env$tr,
          .data$age_stratum == .env$age_stratum,
          !is.na(.data$time),
          !is.na(.data$dLambda0),
          .data$time <= max(times_report) + tol
        )
      
      if (nrow(inc) == 0) {
        return(
          tibble(
            time = numeric(),
            destination = character(),
            dh = numeric()
          )
        )
      }
      
      obj <- prob_model_objects[[tr]]
      
      if (is.null(obj) || inherits(obj$fit, "error")) {
        stop("Invalid model object for transition: ", tr)
      }
      
      beta <- coef(obj$fit)
      
      lp <- calc_lp_times_pi(
        pattern = pattern,
        beta = beta,
        times = inc$time
      )
      
      destination <- meta_origin %>%
        filter(as.character(trans) == tr) %>%
        pull(destination) %>%
        as.character()
      
      tibble(
        time = inc$time,
        destination = destination,
        dh = inc$dLambda0 * exp(lp)
      )
    }
  )
  
  hazard_parts %>%
    group_by(time, destination) %>%
    summarise(
      dh = sum(dh),
      .groups = "drop"
    ) %>%
    arrange(time, destination)
}

# ------------------------------------------------------------
# Standard product-integral / Aalen-Johansen update
#
# For this origin-specific competing-risk process:
#
#   S(t_j) = S(t_j-) * {1 - sum_k dH_k(t_j)}
#   F_k(t_j) = F_k(t_j-) + S(t_j-) dH_k(t_j)
#
# This is the product-integral update for one transient origin
# state with absorbing competing destinations.
#
# IMPORTANT:
# If sum_k dH_k(t_j) > 1, the estimated increment matrix
# I + dA(t_j) is not a valid stochastic transition matrix for
# that predicted covariate pattern. We do NOT truncate, rescale,
# exponentiate, or silently drop that pattern. It is flagged.
# ------------------------------------------------------------

predict_pattern_product_integral <- function(
    origin_name,
    age_stratum,
    sex_value,
    bmi_value,
    hip_value
) {
  meta_origin <- transition_map %>%
    filter(as.character(origin) == origin_name)
  
  dest_states <- unique(as.character(meta_origin$destination))
  
  hz <- make_hazard_increments_pi(
    origin_name = origin_name,
    age_stratum = age_stratum,
    sex_value = sex_value,
    bmi_value = bmi_value,
    hip_value = hip_value
  )
  
  event_times <- sort(unique(hz$time))
  
  S <- 1
  F <- setNames(rep(0, length(dest_states)), dest_states)
  
  event_index <- 1
  invalid <- FALSE
  invalid_time <- NA_real_
  invalid_dtotal <- NA_real_
  max_dtotal <- 0
  
  out <- vector("list", length(times_report))
  
  for (ii in seq_along(times_report)) {
    target <- times_report[ii]
    
    while (
      event_index <= length(event_times) &&
      event_times[event_index] <= target + tol
    ) {
      tt <- event_times[event_index]
      
      d_tt <- hz %>%
        filter(abs(time - tt) < tol)
      
      d_total <- sum(d_tt$dh)
      max_dtotal <- max(max_dtotal, d_total, na.rm = TRUE)
      
      if (
        !is.finite(d_total) ||
        d_total < -tol ||
        d_total > 1 + tol
      ) {
        invalid <- TRUE
        invalid_time <- tt
        invalid_dtotal <- d_total
        break
      }
      
      S_before <- S
      
      if (d_total > tol) {
        for (jj in seq_len(nrow(d_tt))) {
          dest <- d_tt$destination[jj]
          
          F[dest] <- F[dest] +
            S_before * d_tt$dh[jj]
        }
        
        S <- S_before * (1 - d_total)
      }
      
      event_index <- event_index + 1
    }
    
    if (invalid) {
      out[[ii]] <- bind_rows(
        tibble(
          origin = origin_name,
          time_years = target,
          state = origin_name,
          probability = NA_real_
        ),
        tibble(
          origin = origin_name,
          time_years = target,
          state = names(F),
          probability = NA_real_
        )
      ) %>%
        mutate(
          valid_pi = FALSE,
          invalid_time = invalid_time,
          invalid_dtotal = invalid_dtotal,
          max_dtotal_to_target = max_dtotal
        )
      
      # Once invalid, all later requested horizons are invalid too.
      if (ii < length(times_report)) {
        for (kk in (ii + 1):length(times_report)) {
          out[[kk]] <- bind_rows(
            tibble(
              origin = origin_name,
              time_years = times_report[kk],
              state = origin_name,
              probability = NA_real_
            ),
            tibble(
              origin = origin_name,
              time_years = times_report[kk],
              state = names(F),
              probability = NA_real_
            )
          ) %>%
            mutate(
              valid_pi = FALSE,
              invalid_time = invalid_time,
              invalid_dtotal = invalid_dtotal,
              max_dtotal_to_target = max_dtotal
            )
        }
      }
      
      break
    }
    
    out[[ii]] <- bind_rows(
      tibble(
        origin = origin_name,
        time_years = target,
        state = origin_name,
        probability = S
      ),
      tibble(
        origin = origin_name,
        time_years = target,
        state = names(F),
        probability = as.numeric(F)
      )
    ) %>%
      mutate(
        valid_pi = TRUE,
        invalid_time = NA_real_,
        invalid_dtotal = NA_real_,
        max_dtotal_to_target = max_dtotal
      )
  }
  
  bind_rows(out)
}

# ------------------------------------------------------------
# Compute every standardization pattern
# ------------------------------------------------------------

origins <- unique(as.character(transition_map$origin))

pattern_predictions <- map_dfr(
  origins,
  function(origin_name) {
    pats <- standardization_patterns %>%
      filter(as.character(origin) == origin_name)
    
    map_dfr(
      seq_len(nrow(pats)),
      function(i) {
        pat <- pats[i, , drop = FALSE]
        
        map_dfr(
          c(0, 1),
          function(hip_value) {
            predict_pattern_product_integral(
              origin_name = origin_name,
              age_stratum = as.character(pat$age_cat),
              sex_value = as.character(pat$sex_ms),
              bmi_value = as.character(pat$bmi_p1_f_ms),
              hip_value = hip_value
            ) %>%
              mutate(
                group = if_else(
                  hip_value == 1,
                  "first_hip",
                  "first_knee"
                ),
                age_cat = as.character(pat$age_cat),
                sex_ms = as.character(pat$sex_ms),
                bmi_p1_f_ms = as.character(pat$bmi_p1_f_ms),
                pattern_weight = pat$weight,
                n_pattern = pat$n_pattern,
                n_standardization = pat$n_standardization
              )
          }
        )
      }
    )
  }
)

# ------------------------------------------------------------
# Validity summary
# ------------------------------------------------------------

pattern_validity <- pattern_predictions %>%
  group_by(
    origin,
    group,
    age_cat,
    sex_ms,
    bmi_p1_f_ms,
    pattern_weight,
    n_pattern,
    time_years
  ) %>%
  summarise(
    valid_pi = all(valid_pi),
    invalid_time = {
      x <- invalid_time[!is.na(invalid_time)]
      if (length(x) == 0) NA_real_ else min(x)
    },
    invalid_dtotal = {
      x <- invalid_dtotal[!is.na(invalid_dtotal)]
      if (length(x) == 0) NA_real_ else max(x)
    },
    max_dtotal_to_target = max(max_dtotal_to_target, na.rm = TRUE),
    .groups = "drop"
  )

validity_summary <- pattern_validity %>%
  group_by(origin, group, time_years) %>%
  summarise(
    n_patterns = n(),
    n_invalid_patterns = sum(!valid_pi),
    weight_invalid = sum(pattern_weight[!valid_pi]),
    all_patterns_valid = all(valid_pi),
    max_dtotal = max(max_dtotal_to_target, na.rm = TRUE),
    first_invalid_time = {
      x <- invalid_time[!is.na(invalid_time)]
      if (length(x) == 0) NA_real_ else min(x)
    },
    max_invalid_dtotal = {
      x <- invalid_dtotal[!is.na(invalid_dtotal)]
      if (length(x) == 0) NA_real_ else max(x)
    },
    .groups = "drop"
  )

# ------------------------------------------------------------
# Standardize ONLY when all positive-weight patterns are valid.
# We do not redefine the target population by excluding invalid
# patterns.
# ------------------------------------------------------------

pi_standardized <- pattern_predictions %>%
  group_by(
    origin,
    group,
    time_years,
    state
  ) %>%
  summarise(
    all_patterns_valid = all(valid_pi),
    probability = if (
      all(valid_pi)
    ) {
      sum(probability * pattern_weight)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# Probability sums where PI is valid
pi_sum_check <- pi_standardized %>%
  group_by(origin, group, time_years) %>%
  summarise(
    all_patterns_valid = all(all_patterns_valid),
    total_probability = if (
      all(all_patterns_valid)
    ) {
      sum(probability)
    } else {
      NA_real_
    },
    min_probability = if (
      all(all_patterns_valid)
    ) {
      min(probability)
    } else {
      NA_real_
    },
    max_probability = if (
      all(all_patterns_valid)
    ) {
      max(probability)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    abs_error = abs(total_probability - 1),
    valid_probability_vector =
      all_patterns_valid &
      abs_error < 1e-8 &
      min_probability >= -1e-10 &
      max_probability <= 1 + 1e-10
  )

# ------------------------------------------------------------
# Compare with current piecewise-exponential estimates
# ------------------------------------------------------------

pe_for_compare <- pe_report %>%
  filter(time_years %in% times_report) %>%
  transmute(
    origin = as.character(origin),
    group = as.character(group),
    time_years,
    state = as.character(state),
    probability_piecewise = probability
  )

comparison <- pi_standardized %>%
  mutate(
    origin = as.character(origin),
    group = as.character(group),
    state = as.character(state)
  ) %>%
  left_join(
    pe_for_compare,
    by = c(
      "origin",
      "group",
      "time_years",
      "state"
    )
  ) %>%
  mutate(
    difference_pi_minus_piecewise =
      probability - probability_piecewise,
    difference_pp =
      100 * difference_pi_minus_piecewise
  ) %>%
  arrange(
    origin,
    state,
    group,
    time_years
  )

# Contrasts for non-origin destination states
pi_contrasts <- pi_standardized %>%
  mutate(
    origin = as.character(origin),
    group = as.character(group),
    state = as.character(state)
  ) %>%
  filter(state != origin) %>%
  select(
    origin,
    time_years,
    state,
    group,
    probability
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = probability
  ) %>%
  mutate(
    risk_difference_pi =
      first_hip - first_knee,
    risk_difference_pi_pp =
      100 * risk_difference_pi
  ) %>%
  arrange(origin, state, time_years)

pe_contrasts <- pe_for_compare %>%
  filter(state != origin) %>%
  select(
    origin,
    time_years,
    state,
    group,
    probability_piecewise
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = probability_piecewise
  ) %>%
  mutate(
    risk_difference_piecewise =
      first_hip - first_knee,
    risk_difference_piecewise_pp =
      100 * risk_difference_piecewise
  )

contrast_comparison <- pi_contrasts %>%
  left_join(
    pe_contrasts %>%
      select(
        origin,
        time_years,
        state,
        risk_difference_piecewise,
        risk_difference_piecewise_pp
      ),
    by = c(
      "origin",
      "time_years",
      "state"
    )
  ) %>%
  mutate(
    rd_difference_pi_minus_piecewise =
      risk_difference_pi -
      risk_difference_piecewise,
    rd_difference_pp =
      100 * rd_difference_pi_minus_piecewise
  )

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

output_xlsx <- file.path(
  out_dir,
  "p1_standard_product_integral_point_estimates.xlsx"
)

wb <- createWorkbook()

tabs <- list(
  validity_summary = validity_summary,
  standardized_PI = pi_standardized,
  comparison_PI_vs_piecewise = comparison,
  contrast_comparison = contrast_comparison,
  probability_sum_check = pi_sum_check,
  pattern_validity = pattern_validity
)

for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tabs[[nm]])
  setColWidths(wb, nm, cols = 1:50, widths = "auto")
}

saveWorkbook(
  wb,
  output_xlsx,
  overwrite = TRUE
)

save(
  pattern_predictions,
  pattern_validity,
  validity_summary,
  pi_standardized,
  pi_sum_check,
  comparison,
  pi_contrasts,
  contrast_comparison,
  file = file.path(
    out_dir,
    "p1_standard_product_integral_point_estimates.RData"
  )
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STANDARD PRODUCT-INTEGRAL: VALIDITY SUMMARY\n")
cat("============================================\n")
print(validity_summary)

cat("\n============================================\n")
cat("PI vs CURRENT PIECEWISE: DESTINATION STATES\n")
cat("============================================\n")
print(
  comparison %>%
    filter(state != origin) %>%
    select(
      origin,
      state,
      group,
      time_years,
      all_patterns_valid,
      probability,
      probability_piecewise,
      difference_pp
    )
)

cat("\n============================================\n")
cat("RISK-DIFFERENCE COMPARISON\n")
cat("============================================\n")
print(
  contrast_comparison %>%
    select(
      origin,
      state,
      time_years,
      risk_difference_pi_pp,
      risk_difference_piecewise_pp,
      rd_difference_pp
    )
)

cat("\n============================================\n")
cat("PROBABILITY-SUM CHECK\n")
cat("============================================\n")
print(pi_sum_check)

cat("\nSaved:\n")
cat(output_xlsx, "\n")
cat(
  file.path(
    out_dir,
    "p1_standard_product_integral_point_estimates.RData"
  ),
  "\n"
)