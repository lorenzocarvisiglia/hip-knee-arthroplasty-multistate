source("config.R")

# ============================================================
# 08_p1_sensitivity_analyses_updated.R
# sensitivity analyses per analisi P1 multi-state finale
#
# include:
# 1) BMI complete-case sensitivity aligned with file 04
# 2) joint-specific hip/knee revision sensitivity after direct P2
#
# main alignment:
# - age stratified: <60, 60-69, 70-79, 80+
# - hip(t) only where selected in file 04
# - Efron ties
# - patient-cluster robust SE
# - direct P2 only for post-P2 analyses
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(survival)
  library(broom)
  library(openxlsx)
  library(ggplot2)
})

#impostazioni
out_dir <- output_dir
fig_dir <- file.path(out_dir, "figures_p1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
final_model_file <- file.path(out_dir, "p1_multistate_final_models.RData")

if (!file.exists(analysis_file)) {
  stop("non trovo p1_multistate_analysis_objects.RData. lancia prima lo script 03.")
}

if (!file.exists(final_model_file)) {
  stop("non trovo p1_multistate_final_models.RData. lancia prima lo script 04.")
}

load(analysis_file)
load(final_model_file)

times_out <- c(1, 3, 5, 10, 15)
max_time <- max(times_out)
ph_alpha <- 0.05
min_events_for_tv <- 50
tol <- 1e-10

adjusted_transitions <- c(
  "P1_to_P2",
  "P1_to_Rpre",
  "P1_to_D",
  "P2_to_Rpost",
  "P2_to_D",
  "Rpre_to_P2"
)

#funzioni generali
safe_cox <- function(formula, data, keep_model = TRUE, ...) {
  tryCatch(
    coxph(
      formula,
      data = data,
      ties = "efron",
      x = TRUE,
      y = TRUE,
      model = keep_model,
      ...
    ),
    error = function(e) e
  )
}

tt_log <- function(x, t, ...) {
  x * log(t + 1)
}

prepare_transition_data <- function(tr, complete_case_bmi = FALSE) {
  d <- p1_ms_long %>%
    filter(as.character(trans) == tr) %>%
    mutate(
      hip = factor(
        as.character(hip),
        levels = c("first_knee", "first_hip")
      ),
      hip_binary = if_else(hip == "first_hip", 1, 0),
      age_cat = cut(
        age_p1,
        breaks = c(-Inf, 60, 70, 80, Inf),
        right = FALSE,
        labels = c("<60", "60-69", "70-79", "80+")
      ),
      sex_ms = factor(
        as.character(sex_ms),
        levels = c("male", "female")
      ),
      bmi_p1_f_ms = factor(
        as.character(bmi_p1_f_ms),
        levels = c("normal_underweight", "overweight", "obese", "missing")
      )
    ) %>%
    filter(
      !is.na(Tstop),
      Tstop > 0,
      !is.na(status),
      !is.na(hip),
      !is.na(age_cat),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms),
      !is.na(CODPAT)
    )
  
  if (complete_case_bmi) {
    d <- d %>%
      filter(as.character(bmi_p1_f_ms) != "missing") %>%
      mutate(
        bmi_cc = factor(
          as.character(bmi_p1_f_ms),
          levels = c("normal_underweight", "overweight", "obese")
        )
      )
  }
  
  d
}

fit_final_spec <- function(tr, complete_case_bmi = FALSE) {
  d <- prepare_transition_data(tr, complete_case_bmi)
  hip_tv <- tr %in% tv_hip_transitions
  
  if (complete_case_bmi) {
    if (hip_tv) {
      fit <- safe_cox(
        Surv(Tstop, status) ~
          hip_binary + tt(hip_binary) + sex_ms + bmi_cc +
          strata(age_cat) + cluster(CODPAT),
        data = d,
        keep_model = FALSE,
        tt = tt_log
      )
    } else {
      fit <- safe_cox(
        Surv(Tstop, status) ~
          hip + sex_ms + bmi_cc +
          strata(age_cat) + cluster(CODPAT),
        data = d
      )
    }
  } else {
    if (hip_tv) {
      fit <- safe_cox(
        Surv(Tstop, status) ~
          hip_binary + tt(hip_binary) + sex_ms + bmi_p1_f_ms +
          strata(age_cat) + cluster(CODPAT),
        data = d,
        keep_model = FALSE,
        tt = tt_log
      )
    } else {
      fit <- safe_cox(
        Surv(Tstop, status) ~
          hip + sex_ms + bmi_p1_f_ms +
          strata(age_cat) + cluster(CODPAT),
        data = d
      )
    }
  }
  
  list(
    trans = tr,
    fit = fit,
    data = d,
    hip_tv = hip_tv,
    n = nrow(d),
    events = sum(d$status == 1, na.rm = TRUE)
  )
}

extract_hip_effect <- function(obj, model_name, times = times_out) {
  tr <- obj$trans
  fit <- obj$fit
  
  if (inherits(fit, "error")) {
    return(tibble(
      transition = tr,
      model = model_name,
      hip_time_varying = obj$hip_tv,
      time_years = if (obj$hip_tv) times else NA_real_,
      hr = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p = NA_real_,
      error = fit$message
    ))
  }
  
  b <- coef(fit)
  v <- vcov(fit)
  
  if (obj$hip_tv) {
    if (!all(c("hip_binary", "tt(hip_binary)") %in% names(b))) {
      return(tibble(
        transition = tr,
        model = model_name,
        hip_time_varying = TRUE,
        time_years = times,
        hr = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p = NA_real_,
        error = "hip time-varying terms not found"
      ))
    }
    
    return(map_dfr(times, function(tt) {
      l <- log(tt + 1)
      est <- b["hip_binary"] + b["tt(hip_binary)"] * l
      var_est <-
        v["hip_binary", "hip_binary"] +
        l^2 * v["tt(hip_binary)", "tt(hip_binary)"] +
        2 * l * v["hip_binary", "tt(hip_binary)"]
      se <- sqrt(var_est)
      z <- est / se
      
      tibble(
        transition = tr,
        model = model_name,
        hip_time_varying = TRUE,
        time_years = tt,
        hr = exp(est),
        conf_low = exp(est - 1.96 * se),
        conf_high = exp(est + 1.96 * se),
        p = 2 * pnorm(abs(z), lower.tail = FALSE),
        error = NA_character_
      )
    }))
  }
  
  hip_term <- names(b)[str_detect(names(b), "^hip") & !str_detect(names(b), "^tt")]
  
  if (length(hip_term) != 1) {
    return(tibble(
      transition = tr,
      model = model_name,
      hip_time_varying = FALSE,
      time_years = NA_real_,
      hr = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p = NA_real_,
      error = "proportional hip coefficient not uniquely identified"
    ))
  }
  
  term <- hip_term[1]
  est <- b[term]
  se <- sqrt(v[term, term])
  z <- est / se
  
  tibble(
    transition = tr,
    model = model_name,
    hip_time_varying = FALSE,
    time_years = NA_real_,
    hr = exp(est),
    conf_low = exp(est - 1.96 * se),
    conf_high = exp(est + 1.96 * se),
    p = 2 * pnorm(abs(z), lower.tail = FALSE),
    error = NA_character_
  )
}

# ============================================================
# 1) BMI complete-case sensitivity
# ============================================================

message("fit BMI sensitivity con specifica finale del file 04...")

bmi_main_fits <- setNames(
  map(adjusted_transitions, ~ fit_final_spec(.x, complete_case_bmi = FALSE)),
  adjusted_transitions
)

bmi_cc_fits <- setNames(
  map(adjusted_transitions, ~ fit_final_spec(.x, complete_case_bmi = TRUE)),
  adjusted_transitions
)

bmi_main_hip <- map_dfr(
  bmi_main_fits,
  ~ extract_hip_effect(.x, "main_missing_bmi_category")
)

bmi_cc_hip <- map_dfr(
  bmi_cc_fits,
  ~ extract_hip_effect(.x, "bmi_complete_case")
)

bmi_hip_comparison <- bmi_main_hip %>%
  select(
    transition,
    time_years,
    hip_time_varying,
    hr_main = hr,
    low_main = conf_low,
    high_main = conf_high,
    p_main = p,
    error_main = error
  ) %>%
  full_join(
    bmi_cc_hip %>%
      select(
        transition,
        time_years,
        hr_cc = hr,
        low_cc = conf_low,
        high_cc = conf_high,
        p_cc = p,
        error_cc = error
      ),
    by = c("transition", "time_years")
  ) %>%
  mutate(
    hr_ratio_cc_vs_main = hr_cc / hr_main,
    percent_difference = 100 * (hr_cc / hr_main - 1),
    abs_log_hr_difference = abs(log(hr_cc) - log(hr_main))
  ) %>%
  arrange(match(transition, adjusted_transitions), time_years)

bmi_counts <- map_dfr(adjusted_transitions, function(tr) {
  d_main <- bmi_main_fits[[tr]]$data
  d_cc <- bmi_cc_fits[[tr]]$data
  
  tibble(
    transition = tr,
    n_main = nrow(d_main),
    events_main = sum(d_main$status == 1, na.rm = TRUE),
    patients_main = n_distinct(d_main$CODPAT),
    n_cc = nrow(d_cc),
    events_cc = sum(d_cc$status == 1, na.rm = TRUE),
    patients_cc = n_distinct(d_cc$CODPAT),
    excluded_rows = nrow(d_main) - nrow(d_cc),
    excluded_percent = 100 * (nrow(d_main) - nrow(d_cc)) / nrow(d_main)
  )
})

#check against file 04 point estimates
#for PH transitions file 04 repeats the same HR at each reporting time;
#collapse those rows to one NA-time row before joining the refit table
file04_hip <- final_main_hip_hr %>%
  transmute(
    transition,
    hip_time_varying = as.logical(hip_time_varying),
    time_years = if_else(
      hip_time_varying,
      as.numeric(time_years),
      NA_real_
    ),
    hr_file04 = hr,
    low_file04 = conf.low,
    high_file04 = conf.high
  ) %>%
  distinct()

bmi_main_validation <- bmi_main_hip %>%
  select(
    transition,
    hip_time_varying,
    time_years,
    hr_refit = hr,
    low_refit = conf_low,
    high_refit = conf_high
  ) %>%
  left_join(
    file04_hip,
    by = c("transition", "hip_time_varying", "time_years")
  ) %>%
  mutate(
    abs_diff_hr = abs(hr_refit - hr_file04),
    match = !is.na(abs_diff_hr) & abs_diff_hr < 1e-8
  )

if (!all(bmi_main_validation$match)) {
  stop(
    "il refit main del file 08 non coincide con il file 04; controllare bmi_main_validation"
  )
}

# ============================================================
# 2) joint-specific revision sensitivity after direct P2
# ============================================================

message("preparo joint-specific sensitivity dopo P2 diretto...")

p2_risk_ids <- p1_ms_long %>%
  filter(as.character(trans) == "P2_to_Rpost") %>%
  distinct(patient_side_id) %>%
  pull(patient_side_id)

needed_cols <- c(
  "patient_side_id",
  "CODPAT",
  "hip",
  "age_p1",
  "sex_ms",
  "bmi_p1_f_ms",
  "p2_date",
  "hip_revision_date",
  "knee_revision_date",
  "death_date",
  "censor_date"
)

missing_cols <- setdiff(needed_cols, names(p1_ms_wide))
if (length(missing_cols) > 0) {
  stop(
    "mancano colonne in p1_ms_wide: ",
    paste(missing_cols, collapse = ", ")
  )
}

p2_joint_data <- p1_ms_wide %>%
  filter(patient_side_id %in% p2_risk_ids) %>%
  transmute(
    patient_side_id,
    CODPAT,
    hip = factor(as.character(hip), levels = c("first_knee", "first_hip")),
    hip_binary = if_else(hip == "first_hip", 1, 0),
    age_p1,
    age_cat = cut(
      age_p1,
      breaks = c(-Inf, 60, 70, 80, Inf),
      right = FALSE,
      labels = c("<60", "60-69", "70-79", "80+")
    ),
    sex_ms = factor(as.character(sex_ms), levels = c("male", "female")),
    bmi_p1_f_ms = factor(
      as.character(bmi_p1_f_ms),
      levels = c("normal_underweight", "overweight", "obese", "missing")
    ),
    p2_date = as.Date(p2_date),
    hip_revision_date = as.Date(hip_revision_date),
    knee_revision_date = as.Date(knee_revision_date),
    death_date = as.Date(death_date),
    censor_date = as.Date(censor_date)
  ) %>%
  mutate(
    hip_revision_after_p2 = case_when(
      !is.na(hip_revision_date) & hip_revision_date > p2_date ~ hip_revision_date,
      TRUE ~ as.Date(NA)
    ),
    knee_revision_after_p2 = case_when(
      !is.na(knee_revision_date) & knee_revision_date > p2_date ~ knee_revision_date,
      TRUE ~ as.Date(NA)
    ),
    death_after_p2 = case_when(
      !is.na(death_date) & death_date > p2_date ~ death_date,
      TRUE ~ as.Date(NA)
    )
  ) %>%
  rowwise() %>%
  mutate(
    first_date = min(
      c(
        hip_revision_after_p2,
        knee_revision_after_p2,
        death_after_p2,
        censor_date
      ),
      na.rm = TRUE
    ),
    hip_first = !is.na(hip_revision_after_p2) & hip_revision_after_p2 == first_date,
    knee_first = !is.na(knee_revision_after_p2) & knee_revision_after_p2 == first_date,
    death_first = !is.na(death_after_p2) & death_after_p2 == first_date,
    first_event = case_when(
      hip_first & knee_first ~ "both_revisions_same_day",
      hip_first ~ "hip_revision",
      knee_first ~ "knee_revision",
      death_first ~ "death",
      TRUE ~ "censor"
    ),
    time_from_p2 = as.numeric(first_date - p2_date) / 365.25
  ) %>%
  ungroup() %>%
  filter(
    !is.na(time_from_p2),
    time_from_p2 > 0,
    !is.na(hip),
    !is.na(age_cat),
    !is.na(sex_ms),
    !is.na(bmi_p1_f_ms),
    !is.na(CODPAT)
  )

joint_ties <- p2_joint_data %>%
  filter(first_event == "both_revisions_same_day")

if (nrow(joint_ties) > 0) {
  stop(
    paste0(
      "trovati ", nrow(joint_ties),
      " casi con revisione hip e knee nello stesso giorno dopo P2. ",
      "serve una decisione esplicita prima di costruire competing risks joint-specific."
    )
  )
}

joint_event_counts <- p2_joint_data %>%
  count(hip, first_event, name = "n") %>%
  group_by(hip) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

joint_causes <- c("hip_revision", "knee_revision", "death")

make_joint_long <- function(cause) {
  p2_joint_data %>%
    mutate(
      cause = cause,
      Tstop = time_from_p2,
      status = as.integer(first_event == cause)
    )
}

p2_joint_long <- map_dfr(joint_causes, make_joint_long)

joint_transition_counts <- p2_joint_long %>%
  group_by(cause) %>%
  summarise(
    n_rows = n(),
    events = sum(status == 1, na.rm = TRUE),
    patients = n_distinct(CODPAT),
    .groups = "drop"
  )

#diagnostic PH models with age stratification
fit_joint_ph <- function(cause) {
  d <- p2_joint_long %>%
    filter(.data$cause == .env$cause)
  
  fit <- safe_cox(
    Surv(Tstop, status) ~
      hip + sex_ms + bmi_p1_f_ms +
      strata(age_cat) + cluster(CODPAT),
    data = d
  )
  
  list(
    cause = cause,
    data = d,
    fit = fit,
    events = sum(d$status == 1, na.rm = TRUE)
  )
}

joint_ph_fits <- setNames(map(joint_causes, fit_joint_ph), joint_causes)

joint_ph_diagnostics <- map_dfr(joint_ph_fits, function(obj) {
  if (inherits(obj$fit, "error")) {
    return(tibble(
      cause = obj$cause,
      term = NA_character_,
      chisq = NA_real_,
      df = NA_real_,
      p = NA_real_,
      error = obj$fit$message
    ))
  }
  
  z <- tryCatch(cox.zph(obj$fit), error = function(e) e)
  
  if (inherits(z, "error")) {
    return(tibble(
      cause = obj$cause,
      term = NA_character_,
      chisq = NA_real_,
      df = NA_real_,
      p = NA_real_,
      error = z$message
    ))
  }
  
  as.data.frame(z$table) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    transmute(
      cause = obj$cause,
      term,
      chisq,
      df,
      p,
      error = NA_character_
    )
})

joint_model_decision <- joint_transition_counts %>%
  left_join(
    joint_ph_diagnostics %>%
      filter(str_detect(term, "^hip")) %>%
      group_by(cause) %>%
      summarise(
        hip_ph_p = min(p, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "cause"
  ) %>%
  mutate(
    use_hip_tv =
      !is.na(hip_ph_p) &
      hip_ph_p < ph_alpha &
      events >= min_events_for_tv,
    specification = if_else(
      use_hip_tv,
      "hip time-varying",
      "hip proportional"
    )
  )

fit_joint_final <- function(cause) {
  d <- p2_joint_long %>%
    filter(.data$cause == .env$cause)
  
  use_tv <- joint_model_decision %>%
    filter(.data$cause == .env$cause) %>%
    pull(use_hip_tv)
  
  if (use_tv) {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip_binary + tt(hip_binary) + sex_ms + bmi_p1_f_ms +
        strata(age_cat) + cluster(CODPAT),
      data = d,
      keep_model = FALSE,
      tt = tt_log
    )
  } else {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip + sex_ms + bmi_p1_f_ms +
        strata(age_cat) + cluster(CODPAT),
      data = d
    )
  }
  
  list(
    cause = cause,
    fit = fit,
    data = d,
    hip_tv = use_tv,
    events = sum(d$status == 1, na.rm = TRUE)
  )
}

joint_final_fits <- setNames(map(joint_causes, fit_joint_final), joint_causes)

joint_hip_effects <- map_dfr(joint_final_fits, function(obj) {
  wrapper <- list(
    trans = obj$cause,
    fit = obj$fit,
    hip_tv = obj$hip_tv
  )
  extract_hip_effect(wrapper, "joint_specific_final") %>%
    rename(cause = transition)
})

#linear predictor for selected joint-specific models
calc_joint_lp <- function(pattern, beta, times) {
  if (!is.data.frame(pattern)) {
    pattern <- as.data.frame(pattern)
  }
  
  n_pattern <- nrow(pattern)
  n_times <- length(times)
  n <- max(n_pattern, n_times)
  
  if (n_pattern == 1 && n > 1) {
    pattern <- pattern[rep(1, n), , drop = FALSE]
  } else if (n_pattern != n) {
    stop("dimensioni incompatibili tra pattern e times")
  }
  
  if (n_times == 1 && n > 1) {
    times <- rep(times, n)
  } else if (n_times != n) {
    stop("dimensioni incompatibili tra pattern e times")
  }
  
  lp <- rep(0, n)
  hip_value <- as.numeric(pattern$hip_binary)
  sex_value <- as.character(pattern$sex_ms)
  bmi_value <- as.character(pattern$bmi_p1_f_ms)
  
  if ("hipfirst_hip" %in% names(beta)) {
    lp <- lp + beta["hipfirst_hip"] * hip_value
  }
  if ("hip_binary" %in% names(beta)) {
    lp <- lp + beta["hip_binary"] * hip_value
  }
  if ("tt(hip_binary)" %in% names(beta)) {
    lp <- lp + beta["tt(hip_binary)"] * hip_value * log(times + 1)
  }
  if ("sex_msfemale" %in% names(beta)) {
    lp <- lp + beta["sex_msfemale"] * as.numeric(sex_value == "female")
  }
  if ("bmi_p1_f_msoverweight" %in% names(beta)) {
    lp <- lp + beta["bmi_p1_f_msoverweight"] * as.numeric(bmi_value == "overweight")
  }
  if ("bmi_p1_f_msobese" %in% names(beta)) {
    lp <- lp + beta["bmi_p1_f_msobese"] * as.numeric(bmi_value == "obese")
  }
  if ("bmi_p1_f_msmissing" %in% names(beta)) {
    lp <- lp + beta["bmi_p1_f_msmissing"] * as.numeric(bmi_value == "missing")
  }
  
  as.numeric(lp)
}

make_joint_efron_baseline <- function(obj) {
  cause <- obj$cause
  d <- obj$data
  fit <- obj$fit
  
  if (inherits(fit, "error")) {
    return(tibble(
      cause = cause,
      age_stratum = NA_character_,
      time = NA_real_,
      dLambda0 = NA_real_,
      n_events = NA_integer_,
      n_risk = NA_integer_,
      note = fit$message
    ))
  }
  
  beta <- coef(fit)
  
  map_dfr(levels(d$age_cat), function(age_level) {
    d_s <- d %>%
      filter(as.character(age_cat) == age_level)
    
    event_times <- sort(unique(
      d_s$Tstop[d_s$status == 1 & d_s$Tstop <= max_time]
    ))
    
    if (length(event_times) == 0) {
      return(tibble(
        cause = character(),
        age_stratum = character(),
        time = numeric(),
        dLambda0 = numeric(),
        n_events = integer(),
        n_risk = integer(),
        note = character()
      ))
    }
    
    map_dfr(event_times, function(tt) {
      risk_idx <- d_s$Tstop >= tt - tol
      event_idx <- d_s$status == 1 & abs(d_s$Tstop - tt) < tol
      n_events <- sum(event_idx, na.rm = TRUE)
      n_risk <- sum(risk_idx, na.rm = TRUE)
      
      risk_pattern <- d_s[risk_idx, , drop = FALSE]
      event_pattern <- d_s[event_idx, , drop = FALSE]
      
      risk_weight <- sum(exp(calc_joint_lp(
        risk_pattern,
        beta,
        rep(tt, nrow(risk_pattern))
      )))
      
      event_weight <- sum(exp(calc_joint_lp(
        event_pattern,
        beta,
        rep(tt, nrow(event_pattern))
      )))
      
      den <- risk_weight -
        (seq(0, n_events - 1) / n_events) * event_weight
      
      if (any(den <= 0)) {
        return(tibble(
          cause = cause,
          age_stratum = age_level,
          time = tt,
          dLambda0 = NA_real_,
          n_events = n_events,
          n_risk = n_risk,
          note = "non-positive Efron denominator"
        ))
      }
      
      tibble(
        cause = cause,
        age_stratum = age_level,
        time = tt,
        dLambda0 = sum(1 / den),
        n_events = n_events,
        n_risk = n_risk,
        note = NA_character_
      )
    })
  })
}

joint_baseline <- map_dfr(joint_final_fits, make_joint_efron_baseline)

joint_baseline_errors <- joint_baseline %>%
  filter(!is.na(note) | is.na(dLambda0))

if (nrow(joint_baseline_errors) > 0) {
  print(joint_baseline_errors)
  stop("errore nella baseline Efron joint-specific")
}

joint_patterns <- p2_joint_data %>%
  distinct(patient_side_id, .keep_all = TRUE) %>%
  count(age_cat, sex_ms, bmi_p1_f_ms, name = "n_pattern", .drop = FALSE) %>%
  filter(n_pattern > 0) %>%
  mutate(weight = n_pattern / sum(n_pattern))

make_joint_pattern <- function(age_stratum, sex_value, bmi_value, hip_value) {
  tibble(
    age_cat = factor(
      age_stratum,
      levels = c("<60", "60-69", "70-79", "80+")
    ),
    sex_ms = factor(sex_value, levels = c("male", "female")),
    bmi_p1_f_ms = factor(
      bmi_value,
      levels = c("normal_underweight", "overweight", "obese", "missing")
    ),
    hip = factor(
      if_else(hip_value == 1, "first_hip", "first_knee"),
      levels = c("first_knee", "first_hip")
    ),
    hip_binary = as.numeric(hip_value)
  )
}

predict_joint_pattern <- function(age_stratum, sex_value, bmi_value, hip_value) {
  pattern <- make_joint_pattern(age_stratum, sex_value, bmi_value, hip_value)
  
  hazard_parts <- map_dfr(joint_causes, function(cause) {
    inc <- joint_baseline %>%
      filter(
        .data$cause == .env$cause,
        .data$age_stratum == .env$age_stratum,
        !is.na(time),
        !is.na(dLambda0),
        time <= max_time
      )
    
    if (nrow(inc) == 0) {
      return(tibble(time = numeric(), cause = character(), dh = numeric()))
    }
    
    beta <- coef(joint_final_fits[[cause]]$fit)
    lp <- calc_joint_lp(pattern, beta, inc$time)
    
    tibble(
      time = inc$time,
      cause = cause,
      dh = inc$dLambda0 * exp(lp)
    )
  })
  
  event_times <- sort(unique(hazard_parts$time))
  S <- 1
  cif <- setNames(rep(0, length(joint_causes)), joint_causes)
  event_index <- 1
  out <- vector("list", length(times_out))
  
  for (ii in seq_along(times_out)) {
    target <- times_out[ii]
    
    while (
      event_index <= length(event_times) &&
      event_times[event_index] <= target + tol
    ) {
      tt <- event_times[event_index]
      z <- hazard_parts %>%
        filter(abs(time - tt) < tol)
      
      dH <- sum(z$dh)
      
      if (!is.finite(dH) || dH < 0) {
        stop("incremento di hazard non valido nel joint-specific")
      }
      
      if (dH > 0) {
        event_prob <- 1 - exp(-dH)
        S_before <- S
        
        for (jj in seq_len(nrow(z))) {
          cause <- z$cause[jj]
          cif[cause] <- cif[cause] +
            S_before * event_prob * z$dh[jj] / dH
        }
        
        S <- S_before * exp(-dH)
      }
      
      event_index <- event_index + 1
    }
    
    out[[ii]] <- bind_rows(
      tibble(
        time_years = target,
        state = "remain_P2",
        probability = S
      ),
      tibble(
        time_years = target,
        state = names(cif),
        probability = as.numeric(cif)
      )
    )
  }
  
  bind_rows(out)
}

predict_joint_standardized <- function(hip_value) {
  group_name <- if_else(hip_value == 1, "first_hip", "first_knee")
  
  pred <- map_dfr(seq_len(nrow(joint_patterns)), function(i) {
    pat <- joint_patterns[i, , drop = FALSE]
    
    predict_joint_pattern(
      age_stratum = as.character(pat$age_cat),
      sex_value = as.character(pat$sex_ms),
      bmi_value = as.character(pat$bmi_p1_f_ms),
      hip_value = hip_value
    ) %>%
      mutate(weight = pat$weight)
  })
  
  pred %>%
    group_by(time_years, state) %>%
    summarise(
      probability = sum(probability * weight),
      .groups = "drop"
    ) %>%
    mutate(group = group_name)
}

joint_probabilities <- bind_rows(
  predict_joint_standardized(0),
  predict_joint_standardized(1)
) %>%
  mutate(
    group = factor(group, levels = c("first_knee", "first_hip")),
    state_label = case_when(
      state == "hip_revision" ~ "P2 to hip revision",
      state == "knee_revision" ~ "P2 to knee revision",
      state == "death" ~ "P2 to death",
      state == "remain_P2" ~ "Remain in direct P2",
      TRUE ~ state
    )
  )

joint_probability_sum_check <- joint_probabilities %>%
  group_by(group, time_years) %>%
  summarise(
    total = sum(probability),
    abs_error = abs(total - 1),
    .groups = "drop"
  )

joint_contrasts <- joint_probabilities %>%
  filter(state != "remain_P2") %>%
  select(time_years, state, state_label, group, probability) %>%
  pivot_wider(names_from = group, values_from = probability) %>%
  mutate(
    risk_difference = first_hip - first_knee,
    risk_ratio = if_else(first_knee > 0, first_hip / first_knee, NA_real_)
  ) %>%
  arrange(state, time_years)

plot_joint_revision <- joint_probabilities %>%
  filter(state %in% c("hip_revision", "knee_revision")) %>%
  ggplot(
    aes(
      x = time_years,
      y = probability,
      linetype = group,
      group = group
    )
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ state_label, scales = "free_y") +
  labs(
    title = "Joint-specific revision probabilities after direct P2",