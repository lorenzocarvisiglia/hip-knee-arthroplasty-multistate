source("config.R")

# ============================================================
# 22_validate_piecewise_vs_product_integral.R
#
# Validation of the model-based origin-specific probabilities:
# current piecewise-exponential (PE) update versus a canonical
# product-integral / Aalen-Johansen-type (PI) update applied to
# exactly the same fitted cause-specific cumulative-hazard jumps.
#
# This script DOES NOT refit the Cox models. It loads the output
# created by 05_p1_multistate_state_probabilities.R.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(survival)
  library(openxlsx)
})

out_dir <- output_dir
input_file <- file.path(
  out_dir,
  "p1_multistate_state_probabilities.RData"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "File not found: ", input_file, "\n",
      "Run 05_p1_multistate_state_probabilities.R first."
    )
  )
}

load(input_file, envir = .GlobalEnv)

needed_objects <- c(
  "transition_map",
  "prob_model_objects",
  "baseline_increments",
  "standardization_patterns",
  "adjusted_probabilities_report",
  "adjusted_probability_contrasts"
)

missing_objects <- needed_objects[
  !vapply(
    needed_objects,
    function(x) exists(x, envir = .GlobalEnv, inherits = FALSE),
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "Missing required objects in ", input_file, ": ",
      paste(missing_objects, collapse = ", ")
    )
  )
}

# Main reporting horizons used for the validation.
times_validate <- c(5, 10, 15)
max_time <- max(times_validate)
tol <- 1e-10

# Diagnostic threshold only: this is NOT an inferential threshold.
practical_threshold_pp <- 0.10

origin_order <- c("P1", "P2", "R_pre")

# ============================================================
# Helpers reproducing the linear predictor used in file 05
# ============================================================

calc_lp_times <- function(pattern, beta, times) {
  if (!is.data.frame(pattern)) {
    pattern <- as.data.frame(pattern)
  }
  
  n_pattern <- nrow(pattern)
  n_times <- length(times)
  n <- max(n_pattern, n_times)
  
  if (n_pattern == 1 && n > 1) {
    pattern <- pattern[rep(1, n), , drop = FALSE]
  } else if (n_pattern != n) {
    stop("Incompatible dimensions between pattern and times")
  }
  
  if (n_times == 1 && n > 1) {
    times <- rep(times, n)
  } else if (n_times != n) {
    stop("Incompatible dimensions between pattern and times")
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

make_pattern <- function(age_stratum, sex_value, bmi_value, hip_value) {
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
    hip = factor(
      if_else(hip_value == 1, "first_hip", "first_knee"),
      levels = c("first_knee", "first_hip")
    ),
    hip_binary = as.numeric(hip_value)
  )
}

# Build the model-based outgoing cumulative-hazard jumps for one
# covariate pattern and one implant-order scenario.
build_hazard_table <- function(
    origin_name,
    age_stratum,
    sex_value,
    bmi_value,
    hip_value,
    max_time = max_time
) {
  meta_origin <- transition_map %>%
    filter(as.character(origin) == origin_name)
  
  pattern <- make_pattern(
    age_stratum = age_stratum,
    sex_value = sex_value,
    bmi_value = bmi_value,
    hip_value = hip_value
  )
  
  map_dfr(
    meta_origin$trans,
    function(tr) {
      inc <- baseline_increments %>%
        filter(
          .data$trans == .env$tr,
          .data$age_stratum == .env$age_stratum,
          !is.na(.data$time),
          !is.na(.data$dLambda0),
          .data$time <= .env$max_time + tol
        )
      
      if (nrow(inc) == 0) {
        return(
          tibble(
            time = numeric(),
            trans = character(),
            destination = character(),
            dh = numeric()
          )
        )
      }
      
      obj <- prob_model_objects[[tr]]
      
      if (is.null(obj)) {
        stop(paste0("Model object not found for transition: ", tr))
      }
      
      if (inherits(obj$fit, "error")) {
        stop(paste0(tr, ": ", obj$fit$message))
      }
      
      beta <- coef(obj$fit)
      
      if (length(beta) == 0 || any(!is.finite(beta))) {
        stop(paste0("Non-finite coefficients for transition: ", tr))
      }
      
      lp <- calc_lp_times(
        pattern = pattern,
        beta = beta,
        times = inc$time
      )
      
      destination <- meta_origin %>%
        filter(.data$trans == .env$tr) %>%
        pull(destination)
      
      tibble(
        time = inc$time,
        trans = tr,
        destination = as.character(destination),
        dh = inc$dLambda0 * exp(lp)
      )
    }
  ) %>%
    group_by(time, destination) %>%
    summarise(
      dh = sum(dh),
      .groups = "drop"
    ) %>%
    arrange(time, destination)
}

# ============================================================
# Canonical product-integral update
#
# At each event time u_j:
#   dF_k = S(u_j-) dH_k(u_j)
#   S(u_j) = S(u_j-) [1 - sum_k dH_k(u_j)]
#
# This is the direct PI/AJ-type update for the same model-based
# cumulative-hazard jumps used by the PE estimator.
# ============================================================

predict_pattern_pi <- function(
    origin_name,
    age_stratum,
    sex_value,
    bmi_value,
    hip_value,
    times = times_validate
) {
  meta_origin <- transition_map %>%
    filter(as.character(origin) == origin_name)
  
  dest_states <- unique(as.character(meta_origin$destination))
  
  hazard_by_time <- build_hazard_table(
    origin_name = origin_name,
    age_stratum = age_stratum,
    sex_value = sex_value,
    bmi_value = bmi_value,
    hip_value = hip_value,
    max_time = max(times)
  )
  
  all_event_times <- sort(unique(hazard_by_time$time))
  
  survival_origin <- 1
  cif <- setNames(rep(0, length(dest_states)), dest_states)
  event_index <- 1
  pi_valid <- TRUE
  first_invalid_time <- NA_real_
  max_d_total_seen <- 0
  
  out <- vector("list", length(times))
  
  for (ii in seq_along(times)) {
    target <- times[ii]
    
    while (
      pi_valid &&
      event_index <= length(all_event_times) &&
      all_event_times[event_index] <= target + tol
    ) {
      tt <- all_event_times[event_index]
      
      d_tt <- hazard_by_time %>%
        filter(abs(time - tt) < tol)
      
      d_total <- sum(d_tt$dh)
      max_d_total_seen <- max(max_d_total_seen, d_total, na.rm = TRUE)
      
      if (!is.finite(d_total) || d_total < -tol) {
        stop(
          paste0(
            "Invalid cumulative-hazard increment: origin=", origin_name,
            ", age=", age_stratum,
            ", sex=", sex_value,
            ", bmi=", bmi_value,
            ", hip=", hip_value,
            ", t=", signif(tt, 8),
            ", dH=", signif(d_total, 8)
          )
        )
      }
      
      # A canonical product-integral step requires the total jump
      # to be <= 1. If not, the discrete PI update would produce
      # negative survival and is not numerically admissible.
      if (d_total > 1 + tol) {
        pi_valid <- FALSE
        first_invalid_time <- tt
        break
      }
      
      survival_before <- survival_origin
      
      if (d_total > tol) {
        for (jj in seq_len(nrow(d_tt))) {
          dest <- d_tt$destination[jj]
          cif[dest] <-
            cif[dest] +
            survival_before * d_tt$dh[jj]
        }
        
        survival_origin <-
          survival_before * (1 - d_total)
      }
      
      event_index <- event_index + 1
    }
    
    if (pi_valid) {
      out[[ii]] <- bind_rows(
        tibble(
          origin = origin_name,
          time_years = target,
          state = origin_name,
          probability = survival_origin,
          pi_valid = TRUE,
          first_invalid_time = NA_real_,
          max_d_total_seen = max_d_total_seen
        ),
        tibble(
          origin = origin_name,
          time_years = target,
          state = names(cif),
          probability = as.numeric(cif),
          pi_valid = TRUE,
          first_invalid_time = NA_real_,
          max_d_total_seen = max_d_total_seen
        )
      )
    } else {
      out[[ii]] <- bind_rows(
        tibble(
          origin = origin_name,
          time_years = target,
          state = origin_name,
          probability = NA_real_,
          pi_valid = FALSE,
          first_invalid_time = first_invalid_time,
          max_d_total_seen = max_d_total_seen
        ),
        tibble(
          origin = origin_name,
          time_years = target,
          state = dest_states,
          probability = NA_real_,
          pi_valid = FALSE,
          first_invalid_time = first_invalid_time,
          max_d_total_seen = max_d_total_seen
        )
      )
    }
  }
  
  bind_rows(out)
}

# Standardize to exactly the same empirical covariate distribution
# saved by file 05.
predict_standardized_pi <- function(origin_name, hip_value) {
  patterns <- standardization_patterns %>%
    filter(
      as.character(origin) == origin_name,
      weight > 0
    ) %>%
    mutate(
      age_cat = as.character(age_cat),
      sex_ms = as.character(sex_ms),
      bmi_p1_f_ms = as.character(bmi_p1_f_ms)
    )
  
  if (nrow(patterns) == 0) {
    stop(paste0(origin_name, ": empty standardization population"))
  }
  
  group_name <- if_else(
    hip_value == 1,
    "first_hip",
    "first_knee"
  )
  
  pred_patterns <- map_dfr(
    seq_len(nrow(patterns)),
    function(i) {
      pat <- patterns[i, , drop = FALSE]
      
      predict_pattern_pi(
        origin_name = origin_name,
        age_stratum = pat$age_cat,
        sex_value = pat$sex_ms,
        bmi_value = pat$bmi_p1_f_ms,
        hip_value = hip_value,
        times = times_validate
      ) %>%
        mutate(
          pattern_id = i,
          pattern_weight = pat$weight,
          n_pattern = pat$n_pattern,
          age_cat = pat$age_cat,
          sex_ms = pat$sex_ms,
          bmi_p1_f_ms = pat$bmi_p1_f_ms
        )
    }
  )
  
  standardized <- pred_patterns %>%
    group_by(origin, time_years, state) %>%
    summarise(
      pi_valid = all(pi_valid),
      invalid_weight = sum(
        unique(pattern_weight[!pi_valid]),
        na.rm = TRUE
      ),
      probability_pi = if (all(pi_valid)) {
        sum(probability * pattern_weight)
      } else {
        NA_real_
      },
      max_d_total_seen = max(max_d_total_seen, na.rm = TRUE),
      first_invalid_time = if (all(pi_valid)) {
        NA_real_
      } else {
        min(first_invalid_time[!pi_valid], na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    mutate(group = group_name)
  
  list(
    standardized = standardized,
    pattern_predictions = pred_patterns %>%
      mutate(group = group_name)
  )
}

# ============================================================
# Run PI predictions
# ============================================================

message("Computing product-integral validation probabilities...")

pi_results <- list()
pi_pattern_results <- list()

for (origin_name in origin_order) {
  for (hip_value in c(0, 1)) {
    key <- paste(origin_name, hip_value, sep = "_")
    tmp <- predict_standardized_pi(origin_name, hip_value)
    pi_results[[key]] <- tmp$standardized
    pi_pattern_results[[key]] <- tmp$pattern_predictions
  }
}

pi_standardized <- bind_rows(pi_results) %>%
  mutate(
    origin = as.character(origin),
    state = as.character(state),
    group = as.character(group)
  ) %>%
  arrange(origin, group, time_years, state)

pi_pattern_predictions <- bind_rows(pi_pattern_results)

# ============================================================
# Compare against the current PE estimator from file 05
# ============================================================

pe_report <- adjusted_probabilities_report %>%
  mutate(
    origin = as.character(origin),
    state = as.character(state),
    group = as.character(group)
  ) %>%
  filter(time_years %in% times_validate) %>%
  select(
    origin,
    group,
    time_years,
    state,
    probability_pe = probability
  )

comparison <- pe_report %>%
  left_join(
    pi_standardized,
    by = c("origin", "group", "time_years", "state")
  ) %>%
  mutate(
    difference_pi_minus_pe = probability_pi - probability_pe,
    difference_pp = 100 * difference_pi_minus_pe,
    abs_difference_pp = abs(difference_pp)
  ) %>%
  arrange(origin, state, group, time_years)

# Destination states only, because these are the reported clinical estimands.
comparison_destinations <- comparison %>%
  filter(state != origin)

# ============================================================
# Compare first-hip minus first-knee risk differences
# ============================================================

pe_rd <- comparison_destinations %>%
  select(
    origin,
    time_years,
    state,
    group,
    probability_pe
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = probability_pe
  ) %>%
  mutate(
    rd_pe = first_hip - first_knee
  ) %>%
  select(origin, time_years, state, rd_pe)

pi_rd <- comparison_destinations %>%
  select(
    origin,
    time_years,
    state,
    group,
    probability_pi,
    pi_valid
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(probability_pi, pi_valid)
  ) %>%
  mutate(
    pi_valid_both = pi_valid_first_hip & pi_valid_first_knee,
    rd_pi = if_else(
      pi_valid_both,
      probability_pi_first_hip - probability_pi_first_knee,
      NA_real_
    )
  ) %>%
  select(origin, time_years, state, pi_valid_both, rd_pi)

rd_comparison <- pe_rd %>%
  left_join(
    pi_rd,
    by = c("origin", "time_years", "state")
  ) %>%
  mutate(
    rd_pe_pp = 100 * rd_pe,
    rd_pi_pp = 100 * rd_pi,
    rd_difference_pp = rd_pi_pp - rd_pe_pp,
    abs_rd_difference_pp = abs(rd_difference_pp)
  ) %>%
  arrange(origin, state, time_years)

# ============================================================
# Increment diagnostics
# ============================================================

increment_diagnostics <- map_dfr(
  origin_order,
  function(origin_name) {
    patterns <- standardization_patterns %>%
      filter(
        as.character(origin) == origin_name,
        weight > 0
      ) %>%
      mutate(
        age_cat = as.character(age_cat),
        sex_ms = as.character(sex_ms),
        bmi_p1_f_ms = as.character(bmi_p1_f_ms)
      )
    
    map_dfr(
      c(0, 1),
      function(hip_value) {
        group_name <- if_else(
          hip_value == 1,
          "first_hip",
          "first_knee"
        )
        
        map_dfr(
          seq_len(nrow(patterns)),
          function(i) {
            pat <- patterns[i, , drop = FALSE]
            
            h <- build_hazard_table(
              origin_name = origin_name,
              age_stratum = pat$age_cat,
              sex_value = pat$sex_ms,
              bmi_value = pat$bmi_p1_f_ms,
              hip_value = hip_value,
              max_time = max_time
            )
            
            if (nrow(h) == 0) {
              return(
                tibble(
                  origin = origin_name,
                  group = group_name,
                  pattern_id = i,
                  age_cat = pat$age_cat,
                  sex_ms = pat$sex_ms,
                  bmi_p1_f_ms = pat$bmi_p1_f_ms,
                  pattern_weight = pat$weight,
                  max_total_dH = 0,
                  n_steps = 0,
                  n_dH_gt_0_01 = 0,
                  n_dH_gt_0_05 = 0,
                  n_dH_gt_0_10 = 0,
                  n_dH_gt_0_50 = 0,
                  n_dH_gt_1 = 0
                )
              )
            }
            
            totals <- h %>%
              group_by(time) %>%
              summarise(
                total_dH = sum(dh),
                .groups = "drop"
              )
            
            tibble(
              origin = origin_name,
              group = group_name,
              pattern_id = i,
              age_cat = pat$age_cat,
              sex_ms = pat$sex_ms,
              bmi_p1_f_ms = pat$bmi_p1_f_ms,
              pattern_weight = pat$weight,
              max_total_dH = max(totals$total_dH, na.rm = TRUE),
              n_steps = nrow(totals),
              n_dH_gt_0_01 = sum(totals$total_dH > 0.01),
              n_dH_gt_0_05 = sum(totals$total_dH > 0.05),
              n_dH_gt_0_10 = sum(totals$total_dH > 0.10),
              n_dH_gt_0_50 = sum(totals$total_dH > 0.50),
              n_dH_gt_1 = sum(totals$total_dH > 1 + tol)
            )
          }
        )
      }
    )
  }
)

increment_summary <- increment_diagnostics %>%
  group_by(origin, group) %>%
  summarise(
    max_total_dH = max(max_total_dH, na.rm = TRUE),
    patterns_with_dH_gt_0_10 = sum(n_dH_gt_0_10 > 0),
    weighted_mass_with_dH_gt_0_10 = sum(
      pattern_weight[n_dH_gt_0_10 > 0],
      na.rm = TRUE
    ),
    patterns_with_dH_gt_1 = sum(n_dH_gt_1 > 0),
    weighted_mass_with_dH_gt_1 = sum(
      pattern_weight[n_dH_gt_1 > 0],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ============================================================
# Overall diagnostic conclusion
# ============================================================

all_pi_valid <- all(
  comparison_destinations$pi_valid %in% TRUE
)

max_abs_probability_diff_pp <- if (all_pi_valid) {
  max(comparison_destinations$abs_difference_pp, na.rm = TRUE)
} else {
  NA_real_
}

all_rd_valid <- all(rd_comparison$pi_valid_both %in% TRUE)

max_abs_rd_diff_pp <- if (all_rd_valid) {
  max(rd_comparison$abs_rd_difference_pp, na.rm = TRUE)
} else {
  NA_real_
}

max_total_dH <- max(increment_summary$max_total_dH, na.rm = TRUE)

conclusion <- if (!all_pi_valid) {
  paste(
    "The canonical PI update was not numerically admissible for at least one",
    "standardization pattern before 15 years because a model-based total",
    "cumulative-hazard jump exceeded 1. Inspect increment_diagnostics."
  )
} else if (
  max_abs_probability_diff_pp <= practical_threshold_pp &&
  max_abs_rd_diff_pp <= practical_threshold_pp
) {
  paste0(
    "PE and PI estimates are practically indistinguishable at 5, 10 and 15 years: ",
    "maximum absolute difference in a transition probability = ",
    sprintf("%.4f", max_abs_probability_diff_pp),
    " percentage points; maximum absolute difference in a first-hip minus ",
    "first-knee risk difference = ",
    sprintf("%.4f", max_abs_rd_diff_pp),
    " percentage points."
  )
} else {
  paste0(
    "PE and PI estimates differ by more than the diagnostic threshold of ",
    practical_threshold_pp,
    " percentage points for at least one reported quantity. Inspect comparison."
  )
}

validation_summary <- tibble(
  item = c(
    "all_PI_updates_valid_through_15y",
    "maximum_total_model_based_dH",
    "max_abs_probability_difference_pp",
    "max_abs_risk_difference_difference_pp",
    "diagnostic_threshold_pp",
    "conclusion"
  ),
  value = c(
    as.character(all_pi_valid),
    format(max_total_dH, digits = 10),
    if_else(
      is.na(max_abs_probability_diff_pp),
      NA_character_,
      sprintf("%.6f", max_abs_probability_diff_pp)
    ),
    if_else(
      is.na(max_abs_rd_diff_pp),
      NA_character_,
      sprintf("%.6f", max_abs_rd_diff_pp)
    ),
    sprintf("%.2f", practical_threshold_pp),
    conclusion
  )
)

# ============================================================
# Save outputs
# ============================================================

output_xlsx <- file.path(
  out_dir,
  "p1_piecewise_vs_product_integral_validation.xlsx"
)

wb <- createWorkbook()

sheets <- list(
  validation_summary = validation_summary,
  probability_comparison = comparison_destinations,
  risk_difference_comparison = rd_comparison,
  PI_standardized = pi_standardized,
  increment_summary = increment_summary,
  increment_diagnostics = increment_diagnostics,
  PI_pattern_predictions = pi_pattern_predictions
)

for (nm in names(sheets)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, sheets[[nm]])
  setColWidths(wb, nm, cols = 1:100, widths = "auto")
}

saveWorkbook(
  wb,
  output_xlsx,
  overwrite = TRUE
)

save(
  validation_summary,
  comparison,
  comparison_destinations,
  rd_comparison,
  pi_standardized,
  pi_pattern_predictions,
  increment_summary,
  increment_diagnostics,
  file = file.path(
    out_dir,
    "p1_piecewise_vs_product_integral_validation.RData"
  )
)

# ============================================================
# Console output: paste this back into ChatGPT
# ============================================================

cat("\n=============================================\n")
cat("PE vs PRODUCT-INTEGRAL VALIDATION\n")
cat("=============================================\n")
print(validation_summary, n = Inf, width = Inf)

cat("\n=============================================\n")
cat("PROBABILITY COMPARISON (destination states)\n")
cat("=============================================\n")
print(
  comparison_destinations %>%
    select(
      origin,
      state,
      group,
      time_years,
      probability_pe,
      probability_pi,
      difference_pp,
      pi_valid,
      invalid_weight
    ),
  n = Inf,
  width = Inf
)

cat("\n=============================================\n")
cat("RISK-DIFFERENCE COMPARISON\n")
cat("=============================================\n")
print(
  rd_comparison,
  n = Inf,
  width = Inf
)

cat("\n=============================================\n")
cat("INCREMENT SUMMARY\n")
cat("=============================================\n")
print(increment_summary, n = Inf, width = Inf)

cat("\nOutput saved to:\n")
cat(output_xlsx, "\n")