source("config.R")

# ============================================================
# 21_p1_age_at_P2_and_history_sensitivity.R
#
# Sensitivity analysis after direct P2:
# disentangle previous-history duration from attained age at P2.
#
# Motivation:
# In script 20, time P1->P2 was strongly associated with
# post-P2 mortality but not with post-P2 revision. Since
#
#   age at P2 ~= age at P1 + time P1->P2,
#
# the mortality association may reflect attained age at entry
# into P2 rather than residual dependence on the previous
# transition history.
#
# This script:
#   1) reproduces the frozen P2 models from script 04;
#   2) replaces age-at-P1 strata with age-at-P2 strata;
#   3) adds time P1->P2 after adjustment for age at P2;
#   4) checks a nonlinear P1->P2 history effect with a spline;
#   5) for mortality only, uses attained age as the Cox time
#      scale with delayed entry at age P2, with and without the
#      P1->P2 history term.
#
# This is a sensitivity analysis for model adequacy/history
# dependence. It is not a causal adjustment for implant order.
#
# Outputs:
#   analisi_P1/p1_age_at_P2_and_history_sensitivity.xlsx
#   analisi_P1/p1_age_at_P2_and_history_sensitivity.RData
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(survival)
  library(openxlsx)
  library(splines)
})

# ============================================================
# settings
# ============================================================

out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

analysis_file <- file.path(
  out_dir,
  "p1_multistate_analysis_objects.RData"
)

final_models_file <- file.path(
  out_dir,
  "p1_multistate_final_models.RData"
)

out_xlsx <- file.path(
  out_dir,
  "p1_age_at_P2_and_history_sensitivity.xlsx"
)

out_rdata <- file.path(
  out_dir,
  "p1_age_at_P2_and_history_sensitivity.RData"
)

times_report <- c(1, 3, 5, 10, 15)
tol <- 1e-8

if (!file.exists(analysis_file)) {
  stop("Missing: ", analysis_file, ". Run script 03 first.")
}

if (!file.exists(final_models_file)) {
  stop("Missing: ", final_models_file, ". Run script 04 first.")
}

# ============================================================
# load frozen objects
# ============================================================

env03 <- new.env(parent = globalenv())
env04 <- new.env(parent = globalenv())

load(analysis_file, envir = env03)
load(final_models_file, envir = env04)

required_03 <- c("p1_ms_long", "p1_ms_wide")
required_04 <- c("final_main_hip_hr", "final_main_size")

missing_03 <- required_03[
  !vapply(required_03, exists, logical(1), envir = env03)
]

missing_04 <- required_04[
  !vapply(required_04, exists, logical(1), envir = env04)
]

if (length(missing_03) > 0) {
  stop(
    "Objects missing from script 03 RData: ",
    paste(missing_03, collapse = ", ")
  )
}

if (length(missing_04) > 0) {
  stop(
    "Objects missing from script 04 RData: ",
    paste(missing_04, collapse = ", ")
  )
}

p1_ms_long <- as_tibble(env03$p1_ms_long)
p1_ms_wide <- as_tibble(env03$p1_ms_wide)
final_main_hip_hr <- as_tibble(env04$final_main_hip_hr)
final_main_size <- as_tibble(env04$final_main_size)

# ============================================================
# helpers
# ============================================================

tt_log <- function(x, t, ...) {
  x * log(t + 1)
}

safe_cox <- function(formula, data, use_tt = FALSE) {
  tryCatch(
    {
      if (use_tt) {
        coxph(
          formula,
          data = data,
          ties = "efron",
          x = TRUE,
          y = TRUE,
          model = FALSE,
          tt = tt_log
        )
      } else {
        coxph(
          formula,
          data = data,
          ties = "efron",
          x = TRUE,
          y = TRUE,
          model = TRUE
        )
      }
    },
    error = function(e) e
  )
}

assert_fit <- function(fit, label) {
  if (inherits(fit, "error")) {
    stop(label, " failed: ", fit$message)
  }
  invisible(TRUE)
}

extract_tv_hip_hr <- function(
    fit,
    model_label,
    transition = "P2_to_Rpost",
    times = times_report
) {
  b <- coef(fit)
  V <- vcov(fit)
  
  term_main <- "hip_binary"
  term_tv <- "tt(hip_binary)"
  
  if (!all(c(term_main, term_tv) %in% names(b))) {
    stop(
      "Cannot find hip time-varying coefficients in ",
      model_label
    )
  }
  
  map_dfr(
    times,
    function(tt) {
      l <- log(tt + 1)
      
      log_hr <-
        b[term_main] +
        b[term_tv] * l
      
      var_log_hr <-
        V[term_main, term_main] +
        l^2 * V[term_tv, term_tv] +
        2 * l * V[term_main, term_tv]
      
      se <- sqrt(var_log_hr)
      
      tibble(
        transition = transition,
        model = model_label,
        time_years = tt,
        hr = exp(log_hr),
        conf_low = exp(log_hr - 1.96 * se),
        conf_high = exp(log_hr + 1.96 * se)
      )
    }
  )
}

extract_ph_hip_hr <- function(
    fit,
    model_label,
    transition = "P2_to_D"
) {
  b <- coef(fit)
  V <- vcov(fit)
  
  hip_terms <- names(b)[
    str_detect(names(b), "^hip") &
      !str_detect(names(b), "^tt")
  ]
  
  if (length(hip_terms) != 1) {
    stop(
      "Could not identify a unique proportional hip term in ",
      model_label,
      ". Terms: ",
      paste(hip_terms, collapse = ", ")
    )
  }
  
  term <- hip_terms[1]
  beta <- b[term]
  se <- sqrt(V[term, term])
  
  tibble(
    transition = transition,
    model = model_label,
    time_years = NA_real_,
    hr = exp(beta),
    conf_low = exp(beta - 1.96 * se),
    conf_high = exp(beta + 1.96 * se)
  )
}

extract_history <- function(
    fit,
    transition,
    model_label,
    term = "time_p1_p2_5y_c"
) {
  b <- coef(fit)
  V <- vcov(fit)
  
  if (!(term %in% names(b))) {
    stop(
      "History term not found in ",
      model_label,
      ": ", term
    )
  }
  
  beta <- b[term]
  se <- sqrt(V[term, term])
  z <- beta / se
  
  tibble(
    transition = transition,
    model = model_label,
    history_scale = "per 5 additional years from P1 to P2",
    hr = exp(beta),
    conf_low = exp(beta - 1.96 * se),
    conf_high = exp(beta + 1.96 * se),
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}

spline_joint_wald <- function(
    fit,
    transition,
    model_label
) {
  b <- coef(fit)
  V <- vcov(fit)
  
  spline_terms <- grep(
    "ns\\(time_p1_p2_years",
    names(b),
    value = TRUE
  )
  
  if (length(spline_terms) == 0) {
    spline_terms <- grep(
      "splines::ns\\(time_p1_p2_years",
      names(b),
      value = TRUE
    )
  }
  
  if (length(spline_terms) == 0) {
    stop(
      "Could not identify spline terms in ",
      model_label
    )
  }
  
  bs <- b[spline_terms]
  Vs <- V[
    spline_terms,
    spline_terms,
    drop = FALSE
  ]
  
  stat <- as.numeric(
    t(bs) %*% solve(Vs, bs)
  )
  
  df <- length(bs)
  
  tibble(
    transition = transition,
    model = model_label,
    test = "robust joint Wald test for 3-df P1->P2 history spline",
    df = df,
    chi_square = stat,
    p_value = pchisq(
      stat,
      df = df,
      lower.tail = FALSE
    )
  )
}

coef_table <- function(fit, model_label) {
  b <- coef(fit)
  V <- vcov(fit)
  se <- sqrt(diag(V))
  
  tibble(
    model = model_label,
    term = names(b),
    estimate = as.numeric(b),
    robust_se = se,
    hr = exp(as.numeric(b)),
    conf_low = exp(as.numeric(b) - 1.96 * se),
    conf_high = exp(as.numeric(b) + 1.96 * se)
  )
}

# ============================================================
# previous-history and age-at-P2 lookup
# ============================================================

needed_wide <- c(
  "patient_side_id",
  "first_event_from_p1_ms",
  "time_p2_from_p1",
  "age_p1"
)

if (!all(needed_wide %in% names(p1_ms_wide))) {
  stop(
    "p1_ms_wide does not contain: ",
    paste(
      setdiff(needed_wide, names(p1_ms_wide)),
      collapse = ", "
    )
  )
}

has_recorded_age_p2 <- "age_p2" %in% names(p1_ms_wide)

history_lookup <- p1_ms_wide %>%
  transmute(
    patient_side_id = as.character(patient_side_id),
    reached_direct_p2 =
      first_event_from_p1_ms == "p2",
    age_p1_wide = as.numeric(age_p1),
    time_p1_p2_years =
      as.numeric(time_p2_from_p1),
    age_p2_recorded = if (has_recorded_age_p2) {
      as.numeric(age_p2)
    } else {
      NA_real_
    },
    age_p2_reconstructed =
      age_p1_wide + time_p1_p2_years
  ) %>%
  filter(reached_direct_p2) %>%
  mutate(
    age_p2_use = if_else(
      is.finite(age_p2_recorded),
      age_p2_recorded,
      age_p2_reconstructed
    ),
    age_p2_source = if_else(
      is.finite(age_p2_recorded),
      "recorded_at_P2",
      "reconstructed_ageP1_plus_time"
    ),
    age_p2_difference_recorded_minus_reconstructed =
      age_p2_recorded - age_p2_reconstructed
  ) %>%
  distinct(
    patient_side_id,
    .keep_all = TRUE
  )

if (any(!is.finite(history_lookup$time_p1_p2_years))) {
  stop("Non-finite P1->P2 times among direct-P2 units.")
}

if (any(history_lookup$time_p1_p2_years <= 0)) {
  stop("Non-positive P1->P2 times among direct-P2 units.")
}

if (any(!is.finite(history_lookup$age_p2_use))) {
  stop("Could not define age at P2 for all direct-P2 units.")
}

history_median <- median(
  history_lookup$time_p1_p2_years,
  na.rm = TRUE
)

age_p2_median <- median(
  history_lookup$age_p2_use,
  na.rm = TRUE
)

age_p2_qc <- history_lookup %>%
  summarise(
    n_direct_p2 = n(),
    n_recorded_age_p2 =
      sum(age_p2_source == "recorded_at_P2"),
    n_reconstructed_age_p2 =
      sum(age_p2_source == "reconstructed_ageP1_plus_time"),
    median_age_p2 =
      median(age_p2_use),
    p25_age_p2 =
      quantile(age_p2_use, 0.25, names = FALSE),
    p75_age_p2 =
      quantile(age_p2_use, 0.75, names = FALSE),
    median_p1_p2_years =
      median(time_p1_p2_years),
    max_abs_recorded_reconstructed_difference =
      ifelse(
        any(is.finite(age_p2_difference_recorded_minus_reconstructed)),
        max(
          abs(age_p2_difference_recorded_minus_reconstructed),
          na.rm = TRUE
        ),
        NA_real_
      ),
    median_abs_recorded_reconstructed_difference =
      ifelse(
        any(is.finite(age_p2_difference_recorded_minus_reconstructed)),
        median(
          abs(age_p2_difference_recorded_minus_reconstructed),
          na.rm = TRUE
        ),
        NA_real_
      )
  )

# ============================================================
# prepare P2-origin data
# ============================================================

prepare_p2_transition <- function(tr) {
  p1_ms_long %>%
    filter(
      as.character(trans) == tr
    ) %>%
    mutate(
      patient_side_id =
        as.character(patient_side_id),
      
      hip = factor(
        as.character(hip),
        levels = c(
          "first_knee",
          "first_hip"
        )
      ),
      
      hip_binary = if_else(
        hip == "first_hip",
        1,
        0
      ),
      
      age_p1_cat = cut(
        age_p1,
        breaks = c(
          -Inf,
          60,
          70,
          80,
          Inf
        ),
        right = FALSE,
        labels = c(
          "<60",
          "60-69",
          "70-79",
          "80+"
        )
      ),
      
      sex_ms = factor(
        as.character(sex_ms),
        levels = c(
          "male",
          "female"
        )
      ),
      
      bmi_p1_f_ms = factor(
        as.character(bmi_p1_f_ms),
        levels = c(
          "normal_underweight",
          "overweight",
          "obese",
          "missing"
        )
      )
    ) %>%
    left_join(
      history_lookup %>%
        select(
          patient_side_id,
          time_p1_p2_years,
          age_p2_use,
          age_p2_source
        ),
      by = "patient_side_id"
    ) %>%
    mutate(
      time_p1_p2_5y_c =
        (time_p1_p2_years - history_median) / 5,
      
      age_p2_10_c =
        (age_p2_use - age_p2_median) / 10,
      
      age_p2_cat = cut(
        age_p2_use,
        breaks = c(
          -Inf,
          60,
          70,
          80,
          Inf
        ),
        right = FALSE,
        labels = c(
          "<60",
          "60-69",
          "70-79",
          "80+"
        )
      ),
      
      age_entry =
        age_p2_use,
      
      age_exit =
        age_p2_use + Tstop
    ) %>%
    filter(
      !is.na(Tstop),
      Tstop > 0,
      !is.na(status),
      !is.na(hip),
      !is.na(age_p1_cat),
      !is.na(age_p2_cat),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms),
      !is.na(CODPAT),
      is.finite(time_p1_p2_years),
      time_p1_p2_years > 0,
      is.finite(age_entry),
      is.finite(age_exit),
      age_exit > age_entry
    )
}

d_rpost <- prepare_p2_transition(
  "P2_to_Rpost"
)

d_death <- prepare_p2_transition(
  "P2_to_D"
)

# ============================================================
# risk-set QC against frozen script 04
# ============================================================

qc_observed <- tibble(
  transition = c(
    "P2_to_Rpost",
    "P2_to_D"
  ),
  n = c(
    nrow(d_rpost),
    nrow(d_death)
  ),
  events = c(
    sum(d_rpost$status == 1, na.rm = TRUE),
    sum(d_death$status == 1, na.rm = TRUE)
  )
)

qc_expected <- final_main_size %>%
  filter(
    transition %in% c(
      "P2_to_Rpost",
      "P2_to_D"
    ),
    model == "final_main_age_stratified"
  ) %>%
  select(
    transition,
    expected_n = n,
    expected_events = events
  )

qc_riskset <- qc_observed %>%
  left_join(
    qc_expected,
    by = "transition"
  ) %>%
  mutate(
    n_match = n == expected_n,
    events_match = events == expected_events,
    status = if_else(
      n_match & events_match,
      "PASS",
      "CHECK"
    )
  )

if (
  nrow(qc_riskset) != 2 ||
  any(qc_riskset$status != "PASS")
) {
  print(qc_riskset)
  stop(
    "P2 risk sets no longer match frozen script 04."
  )
}

# ============================================================
# descriptive age/history distributions
# ============================================================

age_history_distribution <- d_rpost %>%
  distinct(
    patient_side_id,
    hip,
    age_p2_use,
    time_p1_p2_years
  ) %>%
  group_by(hip) %>%
  summarise(
    n = n(),
    median_age_p2 = median(age_p2_use),
    p25_age_p2 =
      quantile(age_p2_use, 0.25, names = FALSE),
    p75_age_p2 =
      quantile(age_p2_use, 0.75, names = FALSE),
    median_p1_p2_years =
      median(time_p1_p2_years),
    p25_p1_p2_years =
      quantile(time_p1_p2_years, 0.25, names = FALSE),
    p75_p1_p2_years =
      quantile(time_p1_p2_years, 0.75, names = FALSE),
    correlation_ageP2_history =
      cor(age_p2_use, time_p1_p2_years),
    .groups = "drop"
  )

age_p2_strata_counts <- d_rpost %>%
  distinct(
    patient_side_id,
    hip,
    age_p2_cat
  ) %>%
  count(
    age_p2_cat,
    hip,
    name = "n"
  )

# ============================================================
# 1) frozen-specification refits: age at P1 strata
# ============================================================

fit_rpost_ageP1 <- safe_cox(
  Surv(Tstop, status) ~
    hip_binary +
    tt(hip_binary) +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p1_cat) +
    cluster(CODPAT),
  data = d_rpost,
  use_tt = TRUE
)

fit_death_ageP1 <- safe_cox(
  Surv(Tstop, status) ~
    hip +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p1_cat) +
    cluster(CODPAT),
  data = d_death,
  use_tt = FALSE
)

assert_fit(fit_rpost_ageP1, "P2->Rpost age-P1 refit")
assert_fit(fit_death_ageP1, "P2->death age-P1 refit")

# ============================================================
# 2) replace age-at-P1 strata with age-at-P2 strata
# ============================================================

fit_rpost_ageP2 <- safe_cox(
  Surv(Tstop, status) ~
    hip_binary +
    tt(hip_binary) +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p2_cat) +
    cluster(CODPAT),
  data = d_rpost,
  use_tt = TRUE
)

fit_death_ageP2 <- safe_cox(
  Surv(Tstop, status) ~
    hip +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p2_cat) +
    cluster(CODPAT),
  data = d_death,
  use_tt = FALSE
)

assert_fit(fit_rpost_ageP2, "P2->Rpost age-P2 model")
assert_fit(fit_death_ageP2, "P2->death age-P2 model")

# ============================================================
# 3) age at P2 + linear previous-history duration
# ============================================================

fit_rpost_ageP2_history <- safe_cox(
  Surv(Tstop, status) ~
    hip_binary +
    tt(hip_binary) +
    time_p1_p2_5y_c +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p2_cat) +
    cluster(CODPAT),
  data = d_rpost,
  use_tt = TRUE
)

fit_death_ageP2_history <- safe_cox(
  Surv(Tstop, status) ~
    hip +
    time_p1_p2_5y_c +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p2_cat) +
    cluster(CODPAT),
  data = d_death,
  use_tt = FALSE
)

assert_fit(
  fit_rpost_ageP2_history,
  "P2->Rpost age-P2 + history"
)

assert_fit(
  fit_death_ageP2_history,
  "P2->death age-P2 + history"
)

history_after_ageP2 <- bind_rows(
  extract_history(
    fit_rpost_ageP2_history,
    "P2_to_Rpost",
    "age_P2_strata_plus_linear_history"
  ),
  extract_history(
    fit_death_ageP2_history,
    "P2_to_D",
    "age_P2_strata_plus_linear_history"
  )
)

# ============================================================
# 4) age at P2 + nonlinear previous-history duration
# ============================================================

fit_rpost_ageP2_spline_history <- safe_cox(
  Surv(Tstop, status) ~
    hip_binary +
    tt(hip_binary) +
    ns(time_p1_p2_years, df = 3) +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p2_cat) +
    cluster(CODPAT),
  data = d_rpost,
  use_tt = TRUE
)

fit_death_ageP2_spline_history <- safe_cox(
  Surv(Tstop, status) ~
    hip +
    ns(time_p1_p2_years, df = 3) +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_p2_cat) +
    cluster(CODPAT),
  data = d_death,
  use_tt = FALSE
)

assert_fit(
  fit_rpost_ageP2_spline_history,
  "P2->Rpost age-P2 + spline history"
)

assert_fit(
  fit_death_ageP2_spline_history,
  "P2->death age-P2 + spline history"
)

spline_tests_after_ageP2 <- bind_rows(
  spline_joint_wald(
    fit_rpost_ageP2_spline_history,
    "P2_to_Rpost",
    "age_P2_strata_plus_spline_history"
  ),
  spline_joint_wald(
    fit_death_ageP2_spline_history,
    "P2_to_D",
    "age_P2_strata_plus_spline_history"
  )
)

# ============================================================
# 5) mortality with attained age as time scale
#
# This uses delayed entry:
#   entry = age at P2
#   exit  = age at P2 + time since P2
#
# It controls age nonparametrically through the baseline hazard.
# ============================================================

fit_death_attained_age <- safe_cox(
  Surv(age_entry, age_exit, status) ~
    hip +
    sex_ms +
    bmi_p1_f_ms +
    cluster(CODPAT),
  data = d_death,
  use_tt = FALSE
)

fit_death_attained_age_history <- safe_cox(
  Surv(age_entry, age_exit, status) ~
    hip +
    time_p1_p2_5y_c +
    sex_ms +
    bmi_p1_f_ms +
    cluster(CODPAT),
  data = d_death,
  use_tt = FALSE
)

assert_fit(
  fit_death_attained_age,
  "P2->death attained-age model"
)

assert_fit(
  fit_death_attained_age_history,
  "P2->death attained-age + history model"
)

attained_age_history_effect <- extract_history(
  fit_death_attained_age_history,
  "P2_to_D",
  "attained_age_timescale_plus_linear_history"
)

# ============================================================
# implant-order HR comparisons
# ============================================================

rpost_hr_comparison <- bind_rows(
  extract_tv_hip_hr(
    fit_rpost_ageP1,
    "age_P1_strata"
  ),
  extract_tv_hip_hr(
    fit_rpost_ageP2,
    "age_P2_strata"
  ),
  extract_tv_hip_hr(
    fit_rpost_ageP2_history,
    "age_P2_strata_plus_history"
  ),
  extract_tv_hip_hr(
    fit_rpost_ageP2_spline_history,
    "age_P2_strata_plus_spline_history"
  )
) %>%
  arrange(
    time_years,
    model
  )

death_hr_comparison <- bind_rows(
  extract_ph_hip_hr(
    fit_death_ageP1,
    "age_P1_strata"
  ),
  extract_ph_hip_hr(
    fit_death_ageP2,
    "age_P2_strata"
  ),
  extract_ph_hip_hr(
    fit_death_ageP2_history,
    "age_P2_strata_plus_history"
  ),
  extract_ph_hip_hr(
    fit_death_ageP2_spline_history,
    "age_P2_strata_plus_spline_history"
  ),
  extract_ph_hip_hr(
    fit_death_attained_age,
    "attained_age_timescale"
  ),
  extract_ph_hip_hr(
    fit_death_attained_age_history,
    "attained_age_timescale_plus_history"
  )
)

# ============================================================
# reproduce frozen P2->Rpost HRs
# ============================================================

frozen_rpost_hr <- final_main_hip_hr %>%
  filter(
    transition == "P2_to_Rpost",
    time_years %in% times_report
  ) %>%
  transmute(
    time_years,
    frozen_hr = hr
  )

refit_rpost_hr <- extract_tv_hip_hr(
  fit_rpost_ageP1,
  "age_P1_strata"
) %>%
  select(
    time_years,
    refit_hr = hr
  )

main_reproduction <- frozen_rpost_hr %>%
  left_join(
    refit_rpost_hr,
    by = "time_years"
  ) %>%
  mutate(
    abs_difference =
      abs(frozen_hr - refit_hr),
    status = if_else(
      abs_difference < tol,
      "PASS",
      "CHECK"
    )
  )
