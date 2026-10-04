source("config.R")

# ============================================================
# 05_p1_multistate_state_probabilities.R
# probabilità cumulative origin-specific clock-reset
#
# main:
# age stratified
# hip time-varying secondo file 04
#
# standardizzazione:
# stessa distribuzione di età, sesso e bmi nei due gruppi
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(survival)
  library(openxlsx)
  library(ggplot2)
})

#impostazioni
out_dir <- output_dir
fig_dir <- file.path(out_dir, "figures_p1")

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

load(
  file.path(
    out_dir,
    "p1_multistate_analysis_objects.RData"
  )
)

load(
  file.path(
    out_dir,
    "p1_multistate_final_models.RData"
  )
)

times_report <- c(
  1,
  3,
  5,
  10,
  15
)

times_curve <- sort(
  unique(
    c(
      0,
      seq(
        0,
        max(times_report),
        by = 0.25
      ),
      times_report
    )
  )
)

max_time <- max(times_report)

tol <- 1e-10

#specifiche transizioni
transition_map <- tibble(
  trans = c(
    "P1_to_P2",
    "P1_to_Rpre",
    "P1_to_D",
    "P2_to_Rpost",
    "P2_to_D",
    "Rpre_to_P2",
    "Rpre_to_D"
  ),
  
  origin = c(
    "P1",
    "P1",
    "P1",
    "P2",
    "P2",
    "R_pre",
    "R_pre"
  ),
  
  destination = c(
    "P2",
    "R_pre",
    "D",
    "R_post",
    "D",
    "P2_after_Rpre",
    "D"
  ),
  
  inferential_role = c(
    "final_main",
    "final_main",
    "final_main",
    "final_main",
    "final_main",
    "final_main",
    "supporting_competing_risk"
  )
)

#rpre -> d usa hip(t) come modello di supporto dopo la sensitivity 05b
supporting_tv_transitions <- c(
  "Rpre_to_D"
)

probability_tv_transitions <- union(
  tv_hip_transitions,
  supporting_tv_transitions
)

transition_map <- transition_map %>%
  mutate(
    hip_time_varying =
      trans %in% probability_tv_transitions,
    
    model_type = case_when(
      inferential_role ==
        "supporting_competing_risk" &
        hip_time_varying ~
        "supporting_age_stratified_hip_tv",
      
      inferential_role ==
        "supporting_competing_risk" ~
        "supporting_age_stratified_ph",
      
      hip_time_varying ~
        "final_age_stratified_hip_tv",
      
      TRUE ~
        "final_age_stratified_ph"
    )
  )

origin_order <- c(
  "P1",
  "P2",
  "R_pre"
)

#controllo specifica hip(t)
expected_main_tv <- c(
  "P1_to_P2",
  "P2_to_Rpost",
  "Rpre_to_P2"
)

expected_probability_tv <- union(
  expected_main_tv,
  supporting_tv_transitions
)

tv_specification_check <- transition_map %>%
  transmute(
    transition = trans,
    expected_main_tv =
      transition %in% expected_main_tv,
    selected_in_file04 =
      transition %in% tv_hip_transitions,
    supporting_tv_from_05b =
      transition %in% supporting_tv_transitions,
    expected_probability_tv =
      transition %in% expected_probability_tv,
    used_in_file05 =
      hip_time_varying,
    match =
      expected_probability_tv ==
      used_in_file05
  )

if (
  !all(
    tv_specification_check$match
  )
) {
  warning(
    paste(
      "la specificazione hip(t) usata nel file 05",
      "non coincide con quella attesa.",
      "controllare tv_specification_check."
    )
  )
}

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
      y = TRUE,
      model = keep_model,
      ...
    ),
    error = function(e) e
  )
}

tt_log <- function(
    x,
    t,
    ...
) {
  x * log(t + 1)
}

prepare_prob_data <- function(tr) {
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
      !is.na(age_cat),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms),
      !is.na(CODPAT)
    )
}

fit_probability_model <- function(tr) {
  meta <- transition_map %>%
    filter(
      trans == tr
    )
  
  d <- prepare_prob_data(tr)
  
  hip_tv <-
    meta$hip_time_varying[1]
  
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
    origin = meta$origin[1],
    destination =
      meta$destination[1],
    model_type =
      meta$model_type[1],
    inferential_role =
      meta$inferential_role[1],
    hip_tv = hip_tv,
    data = d,
    fit = fit,
    n = nrow(d),
    events = sum(
      d$status == 1,
      na.rm = TRUE
    )
  )
}

#linear predictor vettorializzato
calc_lp_times <- function(
    pattern,
    beta,
    times
) {
  if (
    !is.data.frame(pattern)
  ) {
    pattern <-
      as.data.frame(pattern)
  }
  
  n_pattern <- nrow(pattern)
  n_times <- length(times)
  
  n <- max(
    n_pattern,
    n_times
  )
  
  if (
    n_pattern == 1 &&
    n > 1
  ) {
    pattern <-
      pattern[
        rep(
          1,
          n
        ),
        ,
        drop = FALSE
      ]
  } else if (
    n_pattern != n
  ) {
    stop(
      "dimensioni incompatibili tra pattern e times"
    )
  }
  
  if (
    n_times == 1 &&
    n > 1
  ) {
    times <-
      rep(
        times,
        n
      )
  } else if (
    n_times != n
  ) {
    stop(
      "dimensioni incompatibili tra pattern e times"
    )
  }
  
  hip_value <-
    as.numeric(
      pattern$hip_binary
    )
  
  sex_value <-
    as.character(
      pattern$sex_ms
    )
  
  bmi_value <-
    as.character(
      pattern$bmi_p1_f_ms
    )
  
  lp <- rep(
    0,
    n
  )
  
  if (
    "hipfirst_hip" %in%
    names(beta)
  ) {
    lp <-
      lp +
      beta[
        "hipfirst_hip"
      ] *
      hip_value
  }
  
  if (
    "hip_binary" %in%
    names(beta)
  ) {
    lp <-
      lp +
      beta[
        "hip_binary"
      ] *
      hip_value
  }
  
  if (
    "tt(hip_binary)" %in%
    names(beta)
  ) {
    lp <-
      lp +
      beta[
        "tt(hip_binary)"
      ] *
      hip_value *
      log(
        times + 1
      )
  }
  
  if (
    "sex_msfemale" %in%
    names(beta)
  ) {
    lp <-
      lp +
      beta[
        "sex_msfemale"
      ] *
      as.numeric(
        sex_value ==
          "female"
      )
  }
  
  if (
    "bmi_p1_f_msoverweight" %in%
    names(beta)
  ) {
    lp <-
      lp +
      beta[
        "bmi_p1_f_msoverweight"
      ] *
      as.numeric(
        bmi_value ==
          "overweight"
      )
  }
  
  if (
    "bmi_p1_f_msobese" %in%
    names(beta)
  ) {
    lp <-
      lp +
      beta[
        "bmi_p1_f_msobese"
      ] *
      as.numeric(
        bmi_value ==
          "obese"
      )
  }
  
  if (
    "bmi_p1_f_msmissing" %in%
    names(beta)
  ) {
    lp <-
      lp +
      beta[
        "bmi_p1_f_msmissing"
      ] *
      as.numeric(
        bmi_value ==
          "missing"
      )
  }
  
  as.numeric(lp)
}

#baseline hazard efron per strato di età
make_efron_baseline <- function(
    model_obj
) {
  tr <- model_obj$trans
  d <- model_obj$data
  fit <- model_obj$fit
  
  if (
    inherits(
      fit,
      "error"
    )
  ) {
    return(
      tibble(
        trans = tr,
        age_stratum =
          NA_character_,
        time = NA_real_,
        dLambda0 = NA_real_,
        n_events = NA_integer_,
        n_risk = NA_integer_,
        risk_weight = NA_real_,
        event_weight = NA_real_,
        note = fit$message
      )
    )
  }
  
  beta <- coef(fit)
  
  if (
    any(
      !is.finite(beta)
    )
  ) {
    return(
      tibble(
        trans = tr,
        age_stratum =
          NA_character_,
        time = NA_real_,
        dLambda0 = NA_real_,
        n_events = NA_integer_,
        n_risk = NA_integer_,
        risk_weight = NA_real_,
        event_weight = NA_real_,
        note =
          "non-finite coefficient"
      )
    )
  }
  
  age_levels <-
    levels(
      d$age_cat
    )
  
  map_dfr(
    age_levels,
    function(age_level) {
      d_s <-
        d %>%
        filter(
          as.character(
            age_cat
          ) ==
            age_level
        )
      
      event_times <-
        sort(
          unique(
            d_s$Tstop[
              d_s$status == 1 &
                d_s$Tstop <=
                max_time
            ]
          )
        )
      
      if (
        length(
          event_times
        ) == 0
      ) {
        return(
          tibble(
            trans = character(),
            age_stratum =
              character(),
            time = numeric(),
            dLambda0 = numeric(),
            n_events = integer(),
            n_risk = integer(),
            risk_weight = numeric(),
            event_weight = numeric(),
            note = character()
          )
        )
      }
      
      map_dfr(
        event_times,
        function(tt) {
          risk_idx <-
            d_s$Tstop >=
            tt - tol
          
          event_idx <-
            d_s$status == 1 &
            abs(
              d_s$Tstop -
                tt
            ) < tol
          
          n_events <-
            sum(
              event_idx,
              na.rm = TRUE
            )
          
          n_risk <-
            sum(
              risk_idx,
              na.rm = TRUE
            )
          
          risk_pattern <-
            d_s[
              risk_idx,
              ,
              drop = FALSE
            ]
          
          event_pattern <-
            d_s[
              event_idx,
              ,
              drop = FALSE
            ]
          
          lp_risk <-
            calc_lp_times(
              pattern =
                risk_pattern,
              beta =
                beta,
              times =
                rep(
                  tt,
                  nrow(
                    risk_pattern
                  )
                )
            )
          
          lp_event <-
            calc_lp_times(
              pattern =
                event_pattern,
              beta =
                beta,
              times =
                rep(
                  tt,
                  nrow(
                    event_pattern
                  )
                )
            )
          
          risk_weight <-
            sum(
              exp(
                lp_risk
              )
            )
          
          event_weight <-
            sum(
              exp(
                lp_event
              )
            )
          
          efron_denominators <-
            risk_weight -
            (
              seq(
                0,
                n_events - 1
              ) /
                n_events
            ) *
            event_weight
          
          if (
            any(
              efron_denominators <=
              0
            )
          ) {
            return(
              tibble(
                trans = tr,
                age_stratum =
                  age_level,
                time = tt,
                dLambda0 =
                  NA_real_,
                n_events =
                  n_events,
                n_risk =
                  n_risk,
                risk_weight =
                  risk_weight,
                event_weight =
                  event_weight,
                note =
                  "non-positive Efron denominator"
              )
            )
          }
          
          dLambda0 <-
            sum(
              1 /
                efron_denominators
            )
          
          tibble(
            trans = tr,
            age_stratum =
              age_level,
            time = tt,
            dLambda0 =
              dLambda0,
            n_events =
              n_events,
            n_risk =
              n_risk,
            risk_weight =
              risk_weight,
            event_weight =
              event_weight,
            note =
              NA_character_
          )
        }
      )
    }
  )
}

make_pattern <- function(
    age_stratum,
    sex_value,
    bmi_value,
    hip_value
) {
  tibble(
    age_cat = factor(
      age_stratum,
      levels = c(
        "<60",
        "60-69",
        "70-79",
        "80+"
      )
    ),
    
    sex_ms = factor(
      sex_value,
      levels = c(
        "male",
        "female"
      )
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
      if_else(
        hip_value == 1,
        "first_hip",
        "first_knee"
      ),
      levels = c(
        "first_knee",
        "first_hip"
      )
    ),
    
    hip_binary =
      as.numeric(
        hip_value
      )
  )
}

get_origin_subjects <- function(
    origin_name
) {
  p1_ms_long %>%
    filter(
      as.character(from) ==
        origin_name
    ) %>%
    mutate(
      hip = factor(
        as.character(hip),
        levels = c(
          "first_knee",
          "first_hip"
        )
      ),
      
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
        as.character(
          bmi_p1_f_ms
        ),
        levels = c(
          "normal_underweight",
          "overweight",
          "obese",
          "missing"
        )
      )
    ) %>%
    distinct(
      patient_side_id,
      .keep_all = TRUE
    ) %>%
    filter(
      !is.na(hip),
      !is.na(age_cat),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms)
    )
}

make_pattern_weights <- function(
    origin_name
) {
  d <-
    get_origin_subjects(
      origin_name
    )
  
  n_total <-
    nrow(d)
  
  d %>%
    count(
      age_cat,
      sex_ms,
      bmi_p1_f_ms,
      name =
        "n_pattern",
      .drop = FALSE
    ) %>%
    filter(
      n_pattern > 0
    ) %>%
    mutate(
      origin =
        origin_name,
      
      weight =
        n_pattern /
        n_total,
      
      n_standardization =
        n_total
    )
}

#predizione per singolo pattern
predict_pattern <- function(
    origin_name,
    age_stratum,
    sex_value,
    bmi_value,
    hip_value,
    model_objects,
    baseline_increments,
    times = times_curve
) {
  meta_origin <-
    transition_map %>%
    filter(
      origin ==
        origin_name
    )
  
  dest_states <-
    unique(
      meta_origin$destination
    )
  
  pattern <-
    make_pattern(
      age_stratum =
        age_stratum,
      sex_value =
        sex_value,
      bmi_value =
        bmi_value,
      hip_value =
        hip_value
    )
  
  hazard_parts <-
    map_dfr(
      meta_origin$trans,
      function(tr) {
        inc <-
          baseline_increments %>%
          filter(
            .data$trans == .env$tr,
            .data$age_stratum == .env$age_stratum,
            !is.na(.data$time),
            !is.na(.data$dLambda0),
            .data$time <= max(times)
          )
        