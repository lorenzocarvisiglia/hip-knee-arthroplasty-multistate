source("config.R")

# ============================================================
# 05b_rpre_death_tv_sensitivity.R
# sensitivity analysis per il rischio concorrente Rpre -> D
#
# confronto:
# main: hip proporzionale su Rpre -> D
# sensitivity: hip + hip*log(t+1) su Rpre -> D
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(survival)
  library(broom)
  library(openxlsx)
})

#impostazioni
out_dir <- output_dir

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

load(
  file.path(
    out_dir,
    "p1_multistate_state_probabilities.RData"
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
        15,
        by = 0.25
      ),
      times_report
    )
  )
)

max_time <- max(times_report)
tol <- 1e-10

#funzioni
tt_log <- function(
    x,
    t,
    ...
) {
  x * log(t + 1)
}

prepare_rpre_data <- function(tr) {
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

calc_lp_times <- function(
    pattern,
    beta,
    times
) {
  if (!is.data.frame(pattern)) {
    pattern <- as.data.frame(pattern)
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
    pattern <- pattern[
      rep(1, n),
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
    times <- rep(
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
  
  hip_value <- as.numeric(
    pattern$hip_binary
  )
  
  sex_value <- as.character(
    pattern$sex_ms
  )
  
  bmi_value <- as.character(
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
    lp <- lp +
      beta["hipfirst_hip"] *
      hip_value
  }
  
  if (
    "hip_binary" %in%
    names(beta)
  ) {
    lp <- lp +
      beta["hip_binary"] *
      hip_value
  }
  
  if (
    "tt(hip_binary)" %in%
    names(beta)
  ) {
    lp <- lp +
      beta["tt(hip_binary)"] *
      hip_value *
      log(times + 1)
  }
  
  if (
    "sex_msfemale" %in%
    names(beta)
  ) {
    lp <- lp +
      beta["sex_msfemale"] *
      as.numeric(
        sex_value == "female"
      )
  }
  
  if (
    "bmi_p1_f_msoverweight" %in%
    names(beta)
  ) {
    lp <- lp +
      beta["bmi_p1_f_msoverweight"] *
      as.numeric(
        bmi_value == "overweight"
      )
  }
  
  if (
    "bmi_p1_f_msobese" %in%
    names(beta)
  ) {
    lp <- lp +
      beta["bmi_p1_f_msobese"] *
      as.numeric(
        bmi_value == "obese"
      )
  }
  
  if (
    "bmi_p1_f_msmissing" %in%
    names(beta)
  ) {
    lp <- lp +
      beta["bmi_p1_f_msmissing"] *
      as.numeric(
        bmi_value == "missing"
      )
  }
  
  as.numeric(lp)
}

make_efron_baseline <- function(
    tr,
    d,
    fit
) {
  beta <- coef(fit)
  
  age_levels <- levels(
    d$age_cat
  )
  
  map_dfr(
    age_levels,
    function(age_level) {
      d_s <- d %>%
        filter(
          as.character(age_cat) ==
            age_level
        )
      
      event_times <- sort(
        unique(
          d_s$Tstop[
            d_s$status == 1 &
              d_s$Tstop <= max_time
          ]
        )
      )
      
      if (
        length(event_times) == 0
      ) {
        return(
          tibble(
            trans = character(),
            age_stratum = character(),
            time = numeric(),
            dLambda0 = numeric(),
            n_events = integer(),
            n_risk = integer()
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
              d_s$Tstop - tt
            ) < tol
          
          n_events <- sum(
            event_idx
          )
          
          n_risk <- sum(
            risk_idx
          )
          
          risk_pattern <- d_s[
            risk_idx,
            ,
            drop = FALSE
          ]
          
          event_pattern <- d_s[
            event_idx,
            ,
            drop = FALSE
          ]
          
          lp_risk <- calc_lp_times(
            risk_pattern,
            beta,
            rep(
              tt,
              nrow(risk_pattern)
            )
          )
          
          lp_event <- calc_lp_times(
            event_pattern,
            beta,
            rep(
              tt,
              nrow(event_pattern)
            )
          )
          
          risk_weight <- sum(
            exp(lp_risk)
          )
          
          event_weight <- sum(
            exp(lp_event)
          )
          
          den <- risk_weight -
            (
              seq(
                0,
                n_events - 1
              ) /
                n_events
            ) *
            event_weight
          
          dLambda0 <- sum(
            1 / den
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
              n_risk
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

# ============================================================
# dati rpre
# ============================================================

d_rpre_p2 <- prepare_rpre_data(
  "Rpre_to_P2"
)

d_rpre_d <- prepare_rpre_data(
  "Rpre_to_D"
)

#controllo che i risk set coincidano
riskset_check <- tibble(
  n_rpre_p2 =
    nrow(d_rpre_p2),
  
  n_rpre_d =
    nrow(d_rpre_d),
  
  same_ids =
    setequal(
      d_rpre_p2$patient_side_id,
      d_rpre_d$patient_side_id
    )
)

# ============================================================
# modello Rpre -> P2
# invariato rispetto al main
# ============================================================

fit_rpre_p2 <- coxph(
  Surv(
    Tstop,
    status
  ) ~
    hip_binary +
    tt(hip_binary) +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_cat) +
    cluster(CODPAT),
  data = d_rpre_p2,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = FALSE,
  tt = tt_log
)

# ============================================================
# modello Rpre -> D
# main PH
# ============================================================

fit_rpre_d_ph <- coxph(
  Surv(
    Tstop,
    status
  ) ~
    hip +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_cat) +
    cluster(CODPAT),
  data = d_rpre_d,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE
)

# ============================================================
# modello Rpre -> D
# sensitivity hip(t)
# ============================================================

fit_rpre_d_tv <- coxph(
  Surv(
    Tstop,
    status
  ) ~
    hip_binary +
    tt(hip_binary) +
    sex_ms +
    bmi_p1_f_ms +
    strata(age_cat) +
    cluster(CODPAT),
  data = d_rpre_d,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = FALSE,
  tt = tt_log
)

#risultati modelli
model_results <- bind_rows(
  broom::tidy(
    fit_rpre_d_ph,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      model =
        "Rpre_to_D_PH"
    ),
  
  broom::tidy(
    fit_rpre_d_tv,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      model =
        "Rpre_to_D_hip_tv"
    )
) %>%
  select(
    model,
    everything()
  )

tv_term_test <- broom::tidy(
  fit_rpre_d_tv,
  exponentiate = FALSE,
  conf.int = TRUE
) %>%
  filter(
    term ==
      "tt(hip_binary)"
  )

# ============================================================
# baseline hazard
# ============================================================

base_rpre_p2 <- make_efron_baseline(
  tr = "Rpre_to_P2",
  d = d_rpre_p2,
  fit = fit_rpre_p2
)

base_rpre_d_ph <- make_efron_baseline(
  tr = "Rpre_to_D",
  d = d_rpre_d,
  fit = fit_rpre_d_ph
)

base_rpre_d_tv <- make_efron_baseline(
  tr = "Rpre_to_D",
  d = d_rpre_d,
  fit = fit_rpre_d_tv
)

# ============================================================
# standardizzazione
# ============================================================

standardization_patterns <- d_rpre_p2 %>%
  distinct(
    patient_side_id,
    .keep_all = TRUE
  ) %>%
  count(
    age_cat,
    sex_ms,
    bmi_p1_f_ms,
    name = "n_pattern",
    .drop = FALSE
  ) %>%
  filter(
    n_pattern > 0
  ) %>%
  mutate(
    weight =
      n_pattern /
      sum(n_pattern)
  )

# ============================================================
# predizione rpre per una specificazione del rischio di morte
# ============================================================

predict_one_pattern <- function(
    age_stratum,
    sex_value,
    bmi_value,
    hip_value,
    death_fit,
    death_baseline
) {
  pattern <- make_pattern(
    age_stratum,
    sex_value,
    bmi_value,
    hip_value
  )
  
  beta_p2 <- coef(
    fit_rpre_p2
  )
  
  beta_d <- coef(
    death_fit
  )
  
  inc_p2 <- base_rpre_p2 %>%
    filter(
      .data$age_stratum ==
        .env$age_stratum
    )
  
  inc_d <- death_baseline %>%
    filter(
      .data$age_stratum ==
        .env$age_stratum
    )
  
  hazard_parts <- bind_rows(
    tibble(
      time =
        inc_p2$time,
      
      destination =
        "P2_after_Rpre",
      
      dh =
        inc_p2$dLambda0 *
        exp(
          calc_lp_times(
            pattern,
            beta_p2,
            inc_p2$time
          )
        )
    ),
    
    tibble(
      time =
        inc_d$time,
      
      destination =
        "D",
      
      dh =
        inc_d$dLambda0 *
        exp(
          calc_lp_times(
            pattern,
            beta_d,
            inc_d$time
          )
        )
    )
  ) %>%
    group_by(
      time,
      destination
    ) %>%
    summarise(
      dh = sum(dh),
      .groups = "drop"
    )
  
  all_times <- sort(
    unique(
      hazard_parts$time
    )
  )
  
  S <- 1
  
  cif <- c(
    P2_after_Rpre = 0,
    D = 0
  )
  
  event_index <- 1
  
  out <- vector(
    "list",
    length(times_curve)
  )
  
  for (
    ii in seq_along(
      times_curve
    )
  ) {
    target <- times_curve[ii]
    
    while (
      event_index <=
      length(all_times) &&
      all_times[event_index] <=
      target + tol
    ) {
      tt <- all_times[
        event_index
      ]
      
      z <- hazard_parts %>%
        filter(
          abs(
            time - tt
          ) < tol
        )
      
      dH <- sum(
        z$dh
      )
      
      if (
        dH > 0
      ) {
        event_probability <-
          1 - exp(-dH)
        
        S_before <- S
        
        for (
          jj in seq_len(
            nrow(z)
          )
        ) {
          dest <-
            z$destination[jj]
          
          cif[dest] <-
            cif[dest] +
            S_before *
            event_probability *
            z$dh[jj] /
            dH
        }
        
        S <-
          S_before *
          exp(-dH)
      }
      
      event_index <-
        event_index + 1
    }
    
    out[[ii]] <- bind_rows(
      tibble(
        time_years =
          target,
        state =
          "R_pre",
        probability =
          S
      ),
      
      tibble(
        time_years =
          target,
        state =
          names(cif),
        probability =
          as.numeric(cif)
      )
    )
  }
  
  bind_rows(out)
}

predict_standardized <- function(
    hip_value,
    death_fit,
    death_baseline,
    model_label
) {
  group_name <- if_else(
    hip_value == 1,
    "first_hip",
    "first_knee"
  )
  
  pred <- map_dfr(
    seq_len(
      nrow(
        standardization_patterns
      )
    ),
    function(i) {
      pat <-
        standardization_patterns[
          i,
          ,
          drop = FALSE
        ]
      
      predict_one_pattern(
        age_stratum =
          as.character(
            pat$age_cat
          ),
        
        sex_value =
          as.character(
            pat$sex_ms
          ),
        
        bmi_value =
          as.character(
            pat$bmi_p1_f_ms
          ),
        
        hip_value =
          hip_value,
        
        death_fit =
          death_fit,
        
        death_baseline =
          death_baseline
      ) %>%
        mutate(
          weight =
            pat$weight
        )
    }
  )
  
  pred %>%
    group_by(
      time_years,
      state
    ) %>%
    summarise(
      probability =
        sum(
          probability *
            weight
        ),
      .groups =
        "drop"
    ) %>%
    mutate(
      group =
        group_name,
      model =
        model_label
    )
}

# ============================================================
# predizioni
# ============================================================

pred_main <- bind_rows(
  predict_standardized(
    hip_value = 0,
    death_fit =
      fit_rpre_d_ph,
    death_baseline =
      base_rpre_d_ph,
    model_label =
      "death_PH"
  ),
  
  predict_standardized(
    hip_value = 1,
    death_fit =
      fit_rpre_d_ph,
    death_baseline =
      base_rpre_d_ph,
    model_label =
      "death_PH"
  )
)

pred_tv <- bind_rows(
  predict_standardized(
    hip_value = 0,
    death_fit =
      fit_rpre_d_tv,
    death_baseline =
      base_rpre_d_tv,
    model_label =
      "death_hip_tv"
  ),
  
  predict_standardized(
    hip_value = 1,
    death_fit =
      fit_rpre_d_tv,
    death_baseline =
      base_rpre_d_tv,
    model_label =
      "death_hip_tv"
  )
)

pred_all <- bind_rows(
  pred_main,
  pred_tv
)

# ============================================================
# confronto ai tempi principali
# ============================================================

probability_comparison <- pred_all %>%
  filter(
    time_years %in%
      times_report,
    state ==
      "P2_after_Rpre"
  ) %>%
  select(
    time_years,
    group,
    model,
    probability
  ) %>%
  pivot_wider(
    names_from =
      model,
    values_from =