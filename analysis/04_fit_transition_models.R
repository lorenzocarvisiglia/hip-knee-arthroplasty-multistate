source("config.R")

# ============================================================
# 04_p1_multistate_final_models.R
# modelli finali p1
#
# main analysis:
#   age stratified
#   hip time-varying only when supported after age stratification
#
# sensitivity:
#   continuous age with beta(t)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(survival)
  library(broom)
  library(openxlsx)
})

#impostazioni
out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load(
  file.path(
    out_dir,
    "p1_multistate_analysis_objects.RData"
  )
)

times_hr <- c(1, 3, 5, 10, 15)
ph_alpha <- 0.05

main_transitions <- c(
  "P1_to_P2",
  "P1_to_Rpre",
  "P1_to_D"
)

extended_transitions <- c(
  "P2_to_Rpost",
  "P2_to_D",
  "Rpre_to_P2"
)

descriptive_transitions <- c(
  "Rpre_to_D",
  "Rpost_to_D"
)

model_transitions <- c(
  main_transitions,
  extended_transitions
)

#riferimento per centratura età
age_reference <- p1_ms_long %>%
  distinct(
    patient_side_id,
    age_p1
  ) %>%
  summarise(
    age_median = median(
      age_p1,
      na.rm = TRUE
    )
  ) %>%
  pull(age_median)

age_center_10 <- age_reference / 10

#funzioni
safe_cox <- function(
    formula,
    data,
    keep_model = TRUE,
    ...
) {
  tryCatch(
    coxph(
      formula,
      data = data,
      ties = "efron",
      x = TRUE,
      model = keep_model,
      ...
    ),
    error = function(e) e
  )
}

tidy_cox_safe <- function(
    fit,
    model_name,
    trans_label
) {
  if (inherits(fit, "error")) {
    return(
      tibble(
        model = model_name,
        transition = trans_label,
        term = NA_character_,
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        note = fit$message
      )
    )
  }
  
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      model = model_name,
      transition = trans_label,
      note = NA_character_
    ) %>%
    select(
      model,
      transition,
      term,
      estimate,
      conf.low,
      conf.high,
      p.value,
      note
    )
}

zph_safe <- function(
    fit,
    trans_label,
    model_name
) {
  if (inherits(fit, "error")) {
    return(
      tibble(
        model = model_name,
        transition = trans_label,
        term = NA_character_,
        chisq = NA_real_,
        df = NA_real_,
        p = NA_real_,
        note = fit$message
      )
    )
  }
  
  z <- tryCatch(
    cox.zph(fit),
    error = function(e) e
  )
  
  if (inherits(z, "error")) {
    return(
      tibble(
        model = model_name,
        transition = trans_label,
        term = NA_character_,
        chisq = NA_real_,
        df = NA_real_,
        p = NA_real_,
        note = z$message
      )
    )
  }
  
  as.data.frame(z$table) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(
      model = model_name,
      transition = trans_label,
      note = NA_character_
    ) %>%
    select(
      model,
      transition,
      term,
      chisq,
      df,
      p,
      note
    )
}

safe_min_p <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  min(x)
}

tt_log <- function(x, t, ...) {
  x * log(t + 1)
}

get_model_data <- function(tr) {
  p1_ms_long %>%
    filter(
      as.character(trans) == tr
    ) %>%
    mutate(
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
      
      age_p1_10 = age_p1 / 10,
      
      age_p1_10_c =
        age_p1_10 - age_center_10,
      
      age_cat = cut(
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
    filter(
      !is.na(Tstop),
      Tstop > 0,
      !is.na(status),
      !is.na(hip),
      !is.na(age_p1_10),
      !is.na(age_cat),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms),
      !is.na(CODPAT)
    )
}

model_size_table <- function(
    fits,
    model_name
) {
  map_dfr(
    fits,
    function(x) {
      tibble(
        model = model_name,
        transition = x$trans,
        n = x$n,
        events = x$events,
        hip_time_varying = x$hip_tv
      )
    }
  )
}

# ============================================================
# modello ph standard
# solo riferimento diagnostico
# ============================================================

fit_reference_ph <- function(tr) {
  d <- get_model_data(tr)
  
  fit <- safe_cox(
    Surv(Tstop, status) ~
      hip +
      age_p1_10 +
      sex_ms +
      bmi_p1_f_ms +
      cluster(CODPAT),
    data = d
  )
  
  list(
    trans = tr,
    fit = fit,
    n = nrow(d),
    events = sum(
      d$status == 1,
      na.rm = TRUE
    ),
    hip_tv = FALSE
  )
}

reference_ph_fits <- map(
  model_transitions,
  fit_reference_ph
)

reference_ph_results <- map_dfr(
  reference_ph_fits,
  function(x) {
    tidy_cox_safe(
      x$fit,
      "reference_ph",
      x$trans
    )
  }
)

reference_ph_zph <- map_dfr(
  reference_ph_fits,
  function(x) {
    zph_safe(
      x$fit,
      x$trans,
      "reference_ph"
    )
  }
)

reference_ph_size <- model_size_table(
  reference_ph_fits,
  "reference_ph"
)

reference_ph_diagnostics <- reference_ph_zph %>%
  group_by(transition) %>%
  summarise(
    hip_p_reference = safe_min_p(
      p[
        str_detect(term, "^hip") &
          term != "GLOBAL"
      ]
    ),
    
    age_p_reference = safe_min_p(
      p[
        str_detect(
          term,
          "^age_p1_10"
        )
      ]
    ),
    
    global_p_reference = safe_min_p(
      p[
        term == "GLOBAL"
      ]
    ),
    
    .groups = "drop"
  )

# ============================================================
# modello ph con età stratificata
# questo determina se hip richiede beta(t)
# ============================================================

fit_age_strat_ph <- function(tr) {
  d <- get_model_data(tr)
  
  fit <- safe_cox(
    Surv(Tstop, status) ~
      hip +
      sex_ms +
      bmi_p1_f_ms +
      strata(age_cat) +
      cluster(CODPAT),
    data = d
  )
  
  list(
    trans = tr,
    fit = fit,
    n = nrow(d),
    events = sum(
      d$status == 1,
      na.rm = TRUE
    ),
    hip_tv = FALSE
  )
}

age_strat_ph_fits <- map(
  model_transitions,
  fit_age_strat_ph
)

age_strat_ph_results <- map_dfr(
  age_strat_ph_fits,
  function(x) {
    tidy_cox_safe(
      x$fit,
      "age_stratified_ph_diagnostic",
      x$trans
    )
  }
)

age_strat_ph_zph <- map_dfr(
  age_strat_ph_fits,
  function(x) {
    zph_safe(
      x$fit,
      x$trans,
      "age_stratified_ph_diagnostic"
    )
  }
)

age_strat_ph_size <- model_size_table(
  age_strat_ph_fits,
  "age_stratified_ph_diagnostic"
)

age_strat_ph_diagnostics <- age_strat_ph_zph %>%
  group_by(transition) %>%
  summarise(
    hip_p_age_strat = safe_min_p(
      p[
        str_detect(term, "^hip") &
          term != "GLOBAL"
      ]
    ),
    
    global_p_age_strat = safe_min_p(
      p[
        term == "GLOBAL"
      ]
    ),
    
    .groups = "drop"
  )

# ============================================================
# decisione hip(t)
#
# importante:
# la decisione è basata solo sul modello con età stratificata
# ============================================================

ph_diagnostic_comparison <- reference_ph_diagnostics %>%
  left_join(
    age_strat_ph_diagnostics,
    by = "transition"
  ) %>%
  mutate(
    hip_nonph_reference =
      !is.na(hip_p_reference) &
      hip_p_reference < ph_alpha,
    
    age_nonph_reference =
      !is.na(age_p_reference) &
      age_p_reference < ph_alpha,
    
    hip_nonph_age_strat =
      !is.na(hip_p_age_strat) &
      hip_p_age_strat < ph_alpha,
    
    use_hip_time_varying =
      hip_nonph_age_strat,
    
    decision_basis =
      "PH diagnostic after age stratification"
  ) %>%
  arrange(
    match(
      transition,
      model_transitions
    )
  )

tv_hip_transitions <- ph_diagnostic_comparison %>%
  filter(
    use_hip_time_varying
  ) %>%
  pull(transition)

# ============================================================
# numerosità degli strata di età
# ============================================================

age_strata_counts <- map_dfr(
  model_transitions,
  function(tr) {
    d <- get_model_data(tr)
    
    d %>%
      group_by(age_cat) %>%
      summarise(
        transition = tr,
        n_risk_rows = n(),
        n_patients = n_distinct(CODPAT),
        events = sum(
          status == 1,
          na.rm = TRUE
        ),
        event_percent =
          100 * events / n_risk_rows,
        .groups = "drop"
      ) %>%
      select(
        transition,
        everything()
      )
  }
)

age_strata_counts_by_hip <- map_dfr(
  model_transitions,
  function(tr) {
    d <- get_model_data(tr)
    
    d %>%
      group_by(
        age_cat,
        hip
      ) %>%
      summarise(
        transition = tr,
        n_risk_rows = n(),
        n_patients = n_distinct(CODPAT),
        events = sum(
          status == 1,
          na.rm = TRUE
        ),
        event_percent =
          100 * events / n_risk_rows,
        .groups = "drop"
      ) %>%
      select(
        transition,
        everything()
      )
  }
)

# ============================================================
# modello principale
# età stratificata
# ============================================================

fit_final_main <- function(tr) {
  d <- get_model_data(tr)
  
  hip_tv <- tr %in% tv_hip_transitions
  
  if (hip_tv) {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip_binary +
        tt(hip_binary) +
        sex_ms +
        bmi_p1_f_ms +
        strata(age_cat) +
        cluster(CODPAT),
      data = d,
      keep_model = FALSE,
      tt = tt_log
    )
  } else {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip +
        sex_ms +
        bmi_p1_f_ms +
        strata(age_cat) +
        cluster(CODPAT),
      data = d
    )
  }
  
  list(
    trans = tr,
    fit = fit,
    n = nrow(d),
    events = sum(
      d$status == 1,
      na.rm = TRUE
    ),
    hip_tv = hip_tv
  )
}

final_main_fits <- map(
  model_transitions,
  fit_final_main
)

final_main_results <- map_dfr(
  final_main_fits,
  function(x) {
    tidy_cox_safe(
      x$fit,
      "final_main_age_stratified",
      x$trans
    )
  }
)

final_main_size <- model_size_table(
  final_main_fits,
  "final_main_age_stratified"
)

# ============================================================
# sensitivity
# età continua con beta(t)
#
# usa la stessa decisione su hip(t) del modello principale
# ============================================================

fit_age_tv_sensitivity <- function(tr) {
  d <- get_model_data(tr)
  
  hip_tv <- tr %in% tv_hip_transitions
  
  if (hip_tv) {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip_binary +
        tt(hip_binary) +
        age_p1_10_c +
        tt(age_p1_10_c) +
        sex_ms +
        bmi_p1_f_ms +
        cluster(CODPAT),
      data = d,
      keep_model = FALSE,
      tt = tt_log
    )
  } else {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip +
        age_p1_10_c +
        tt(age_p1_10_c) +
        sex_ms +
        bmi_p1_f_ms +
        cluster(CODPAT),
      data = d,
      keep_model = FALSE,
      tt = tt_log
    )
  }
  
  list(
    trans = tr,
    fit = fit,
    n = nrow(d),
    events = sum(
      d$status == 1,
      na.rm = TRUE
    ),
    hip_tv = hip_tv
  )
}

age_tv_sensitivity_fits <- map(
  model_transitions,
  fit_age_tv_sensitivity
)

age_tv_sensitivity_results <- map_dfr(
  age_tv_sensitivity_fits,
  function(x) {
    tidy_cox_safe(
      x$fit,
      "sensitivity_age_time_varying",
      x$trans
    )
  }
)

age_tv_sensitivity_size <- model_size_table(
  age_tv_sensitivity_fits,
  "sensitivity_age_time_varying"
)

# ============================================================
# funzione hr hip
# ============================================================

extract_hip_hr <- function(
    fit_obj,
    model_name,
    times = times_hr
) {
  tr <- fit_obj$trans
  fit <- fit_obj$fit
  hip_tv <- fit_obj$hip_tv
  
  if (inherits(fit, "error")) {
    return(
      tibble(
        model = model_name,
        transition = tr,
        hip_time_varying = hip_tv,
        time_years = times,
        hr = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        note = fit$message
      )
    )
  }
  
  b <- coef(fit)
  v <- vcov(fit)
  
  if (hip_tv) {
    term_main <- "hip_binary"
    term_tv <- "tt(hip_binary)"
    
    if (
      !(term_main %in% names(b)) |
      !(term_tv %in% names(b))
    ) {
      return(
        tibble(
          model = model_name,
          transition = tr,
          hip_time_varying = TRUE,
          time_years = times,
          hr = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          p.value = NA_real_,
          note = "hip time-varying terms not found"
        )
      )
    }
    
    return(
      map_dfr(
        times,
        function(t) {
          l <- log(t + 1)
          
          beta_t <-
            b[term_main] +
            b[term_tv] * l
          
          var_t <-
            v[term_main, term_main] +
            l^2 * v[term_tv, term_tv] +
            2 * l *
            v[term_main, term_tv]
          
          se_t <- sqrt(var_t)
          z_t <- beta_t / se_t
          
          tibble(
            model = model_name,
            transition = tr,
            hip_time_varying = TRUE,
            time_years = t,
            hr = exp(beta_t),
            conf.low = exp(
              beta_t - 1.96 * se_t
            ),
            conf.high = exp(
              beta_t + 1.96 * se_t
            ),
            p.value = 2 * pnorm(
              abs(z_t),
              lower.tail = FALSE
            ),
            note = NA_character_
          )
        }
      )
    )
  }
  
  hip_terms <- names(b)[
    str_detect(
      names(b),
      "^hip"
    ) &
      !str_detect(
        names(b),
        "^tt"
      )
  ]
  
  if (length(hip_terms) != 1) {
    return(
      tibble(
        model = model_name,
        transition = tr,
        hip_time_varying = FALSE,
        time_years = times,
        hr = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        note =
          "proportional hip coefficient not uniquely identified"
      )
    )
  }
  
  term_main <- hip_terms[1]
  
  beta <- b[term_main]
  
  se <- sqrt(
    v[
      term_main,
      term_main
    ]
  )
  
  z <- beta / se
  
  tibble(
    model = model_name,
    transition = tr,
    hip_time_varying = FALSE,
    time_years = times,
    hr = exp(beta),
    conf.low = exp(
      beta - 1.96 * se
    ),
    conf.high = exp(
      beta + 1.96 * se
    ),
    p.value = 2 * pnorm(
      abs(z),
      lower.tail = FALSE
    ),
    note = NA_character_
  )
}

final_main_hip_hr <- map_dfr(
  final_main_fits,
  extract_hip_hr,
  model_name =
    "final_main_age_stratified"
)

age_tv_sensitivity_hip_hr <- map_dfr(
  age_tv_sensitivity_fits,
  extract_hip_hr,
  model_name =
    "sensitivity_age_time_varying"
)

# ============================================================
# effetto età nel modello sensitivity
# ============================================================

extract_age_hr <- function(
    fit_obj,
    times = times_hr
) {
  tr <- fit_obj$trans
  fit <- fit_obj$fit
  
  if (inherits(fit, "error")) {
    return(
      tibble(
        transition = tr,
        time_years = times,
        hr_per_10_years = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        note = fit$message
      )
    )
  }
  
  b <- coef(fit)
  v <- vcov(fit)
  
  term_main <- "age_p1_10_c"
  term_tv <- "tt(age_p1_10_c)"
  
  if (
    !(term_main %in% names(b)) |
    !(term_tv %in% names(b))
  ) {
    return(
      tibble(
        transition = tr,
        time_years = times,
        hr_per_10_years = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        note =
          "age time-varying terms not found"
      )
    )
  }
  
  map_dfr(
    times,
    function(t) {
      l <- log(t + 1)
      
      beta_t <-
        b[term_main] +
        b[term_tv] * l
      
      var_t <-
        v[term_main, term_main] +
        l^2 * v[term_tv, term_tv] +
        2 * l *
        v[term_main, term_tv]
      
      se_t <- sqrt(var_t)
      z_t <- beta_t / se_t
      
      tibble(
        transition = tr,
        time_years = t,
        hr_per_10_years = exp(beta_t),
        conf.low = exp(
          beta_t - 1.96 * se_t
        ),
        conf.high = exp(
          beta_t + 1.96 * se_t
        ),
        p.value = 2 * pnorm(
          abs(z_t),
          lower.tail = FALSE
        ),
        note = NA_character_
      )
    }
  )
}

age_tv_sensitivity_age_hr <- map_dfr(
  age_tv_sensitivity_fits,
  extract_age_hr
)

# ============================================================
# test termini tempo-variabili
# ============================================================

hip_tv_term_tests <- bind_rows(
  final_main_results,
  age_tv_sensitivity_results
) %>%
  filter(
    term == "tt(hip_binary)"
  ) %>%
  arrange(
    transition,
    model
  )

age_tv_term_tests <- age_tv_sensitivity_results %>%
  filter(
    term %in% c(
      "age_p1_10_c",
      "tt(age_p1_10_c)"