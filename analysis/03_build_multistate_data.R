# 03_make_p1_multistate_analysis.R
# analisi multi-state da p1, opzione 3

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(survival)
  library(broom)
  library(openxlsx)
})

#impostazioni
out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(clean_data_file)) {
  stop("Missing cleaned data file: ", clean_data_file, ". Run analysis/01_clean_data.R first.")
}
load(clean_data_file)

min_events_adjusted <- 100
min_events_simple <- 30

dat0 <- cohort_p1_clean

#funzioni
to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  
  if (is.numeric(x)) {
    return(as.Date(x, origin = "1899-12-30"))
  }
  
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c("", "-", "NA", "NaN", "NULL", "NULLA")] <- NA_character_
  
  out <- suppressWarnings(as.Date(as.numeric(x_chr), origin = "1899-12-30"))
  
  idx <- is.na(out) & !is.na(x_chr)
  if (any(idx)) out[idx] <- suppressWarnings(ymd(x_chr[idx]))
  
  idx <- is.na(out) & !is.na(x_chr)
  if (any(idx)) out[idx] <- suppressWarnings(dmy(x_chr[idx]))
  
  idx <- is.na(out) & !is.na(x_chr)
  if (any(idx)) out[idx] <- suppressWarnings(mdy(x_chr[idx]))
  
  out
}

years_between <- function(end, start) {
  as.numeric(difftime(as.Date(end), as.Date(start), units = "days")) / 365.25
}

pmin_na <- function(...) {
  m <- cbind(...)
  out <- apply(m, 1, function(z) {
    z <- z[!is.na(z)]
    if (length(z) == 0) NA_real_ else min(z)
  })
  as.numeric(out)
}

first_at_min <- function(time, min_time) {
  !is.na(time) & !is.na(min_time) & abs(time - min_time) < 1e-8
}

clean_bmi_cat <- function(x) {
  out <- case_when(
    x %in% c(
      "normal_underweight",
      "normal",
      "normal/underweight",
      "Normal or underweight BMI",
      "Sottopeso o Normopeso"
    ) ~ "normal_underweight",
    x %in% c(
      "overweight",
      "Overweight BMI",
      "Sovrappeso"
    ) ~ "overweight",
    x %in% c(
      "obese",
      "Obese BMI",
      "Obeso"
    ) ~ "obese",
    is.na(x) ~ "missing",
    TRUE ~ as.character(x)
  )
  
  factor(
    out,
    levels = c(
      "normal_underweight",
      "overweight",
      "obese",
      "missing"
    )
  )
}

clean_sex <- function(x) {
  case_when(
    x %in% c(
      "F",
      "female",
      "Female",
      "femmina",
      "Femmina",
      "1"
    ) ~ "female",
    x %in% c(
      "M",
      "male",
      "Male",
      "maschio",
      "Maschio",
      "0"
    ) ~ "male",
    TRUE ~ as.character(x)
  ) %>%
    factor(levels = c("male", "female"))
}

safe_cox <- function(formula, data) {
  tryCatch(
    coxph(
      formula,
      data = data,
      ties = "efron",
      x = TRUE,
      model = TRUE
    ),
    error = function(e) e
  )
}

tidy_cox_safe <- function(fit, model_name, trans_label) {
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

#preparazione base
dat <- dat0 %>%
  mutate(
    p1_date = to_date(p1_date),
    p2_date = to_date(p2_date),
    hip_primary_date = to_date(hip_primary_date),
    knee_primary_date = to_date(knee_primary_date),
    hip_revision_date = to_date(hip_revision_date),
    knee_revision_date = to_date(knee_revision_date),
    first_revision_after_p1_date = to_date(first_revision_after_p1_date),
    death_date = to_date(death_date),
    censor_date = to_date(censor_date),
    
    hip = as.numeric(hip),
    has_p2_by_admin = as.logical(has_p2_by_admin),
    
    order_group_ms = case_when(
      hip == 1 ~ "first_hip",
      hip == 0 ~ "first_knee",
      TRUE ~ "undefined"
    ),
    order_group_ms = factor(
      order_group_ms,
      levels = c(
        "first_knee",
        "first_hip",
        "undefined"
      )
    ),
    
    sex_ms = clean_sex(sex),
    bmi_p1_f_ms = clean_bmi_cat(bmi_p1_f),
    
    p1_year = year(p1_date),
    
    time_censor_from_p1 = years_between(
      censor_date,
      p1_date
    ),
    time_p2_from_p1 = years_between(
      p2_date,
      p1_date
    ),
    time_death_from_p1 = years_between(
      death_date,
      p1_date
    ),
    time_hip_revision_from_p1 = years_between(
      hip_revision_date,
      p1_date
    ),
    time_knee_revision_from_p1 = years_between(
      knee_revision_date,
      p1_date
    )
  ) %>%
  mutate(
    time_censor_from_p1 = if_else(
      time_censor_from_p1 > 0,
      time_censor_from_p1,
      NA_real_
    ),
    time_p2_from_p1 = if_else(
      time_p2_from_p1 > 0,
      time_p2_from_p1,
      NA_real_
    ),
    time_death_from_p1 = if_else(
      time_death_from_p1 > 0,
      time_death_from_p1,
      NA_real_
    ),
    time_hip_revision_from_p1 = if_else(
      time_hip_revision_from_p1 > 0,
      time_hip_revision_from_p1,
      NA_real_
    ),
    time_knee_revision_from_p1 = if_else(
      time_knee_revision_from_p1 > 0,
      time_knee_revision_from_p1,
      NA_real_
    ),
    
    time_first_revision_from_p1 = pmin_na(
      time_hip_revision_from_p1,
      time_knee_revision_from_p1
    ),
    
    first_revision_joint_from_p1 = case_when(
      first_at_min(
        time_hip_revision_from_p1,
        time_first_revision_from_p1
      ) ~ "hip_revision",
      first_at_min(
        time_knee_revision_from_p1,
        time_first_revision_from_p1
      ) ~ "knee_revision",
      TRUE ~ NA_character_
    ),
    
    time_first_from_p1 = pmin_na(
      time_p2_from_p1,
      time_first_revision_from_p1,
      time_death_from_p1,
      time_censor_from_p1
    ),
    
    first_event_from_p1_ms = case_when(
      first_at_min(
        time_p2_from_p1,
        time_first_from_p1
      ) ~ "p2",
      first_at_min(
        time_first_revision_from_p1,
        time_first_from_p1
      ) ~ "revision",
      first_at_min(
        time_death_from_p1,
        time_first_from_p1
      ) ~ "death",
      first_at_min(
        time_censor_from_p1,
        time_first_from_p1
      ) ~ "censor",
      TRUE ~ NA_character_
    )
  )

#esclusione ordine indefinito per modelli con hip
analysis_dat <- dat %>%
  filter(!is.na(hip)) %>%
  mutate(
    hip = factor(
      hip,
      levels = c(0, 1),
      labels = c(
        "first_knee",
        "first_hip"
      )
    ),
    hip_binary = if_else(
      hip == "first_hip",
      1,
      0
    )
  )

#tempi dopo p2
analysis_dat <- analysis_dat %>%
  mutate(
    time_hip_revision_from_p2 = years_between(
      hip_revision_date,
      p2_date
    ),
    time_knee_revision_from_p2 = years_between(
      knee_revision_date,
      p2_date
    ),
    time_death_from_p2 = years_between(
      death_date,
      p2_date
    ),
    time_censor_from_p2 = years_between(
      censor_date,
      p2_date
    ),
    
    time_hip_revision_from_p2 = if_else(
      time_hip_revision_from_p2 > 0,
      time_hip_revision_from_p2,
      NA_real_
    ),
    time_knee_revision_from_p2 = if_else(
      time_knee_revision_from_p2 > 0,
      time_knee_revision_from_p2,
      NA_real_
    ),
    time_death_from_p2 = if_else(
      time_death_from_p2 > 0,
      time_death_from_p2,
      NA_real_
    ),
    time_censor_from_p2 = if_else(
      time_censor_from_p2 > 0,
      time_censor_from_p2,
      NA_real_
    ),
    
    time_first_revision_from_p2 = pmin_na(
      time_hip_revision_from_p2,
      time_knee_revision_from_p2
    ),
    
    first_revision_joint_from_p2 = case_when(
      first_at_min(
        time_hip_revision_from_p2,
        time_first_revision_from_p2
      ) ~ "hip_revision",
      first_at_min(
        time_knee_revision_from_p2,
        time_first_revision_from_p2
      ) ~ "knee_revision",
      TRUE ~ NA_character_
    ),
    
    time_first_from_p2 = pmin_na(
      time_first_revision_from_p2,
      time_death_from_p2,
      time_censor_from_p2
    ),
    
    first_event_from_p2_ms = case_when(
      first_at_min(
        time_first_revision_from_p2,
        time_first_from_p2
      ) ~ "revision",
      first_at_min(
        time_death_from_p2,
        time_first_from_p2
      ) ~ "death",
      first_at_min(
        time_censor_from_p2,
        time_first_from_p2
      ) ~ "censor",
      TRUE ~ NA_character_
    )
  )

#tempi dopo revisione pre-p2
analysis_dat <- analysis_dat %>%
  mutate(
    revision_pre_p2_date = case_when(
      first_event_from_p1_ms == "revision" &
        first_revision_joint_from_p1 == "hip_revision" ~ hip_revision_date,
      first_event_from_p1_ms == "revision" &
        first_revision_joint_from_p1 == "knee_revision" ~ knee_revision_date,
      TRUE ~ as.Date(NA)
    ),
    
    time_p2_from_rpre = years_between(
      p2_date,
      revision_pre_p2_date
    ),
    time_death_from_rpre = years_between(
      death_date,
      revision_pre_p2_date
    ),
    time_censor_from_rpre = years_between(
      censor_date,
      revision_pre_p2_date
    ),
    
    time_p2_from_rpre = if_else(
      time_p2_from_rpre > 0,
      time_p2_from_rpre,
      NA_real_
    ),
    time_death_from_rpre = if_else(
      time_death_from_rpre > 0,
      time_death_from_rpre,
      NA_real_
    ),
    time_censor_from_rpre = if_else(
      time_censor_from_rpre > 0,
      time_censor_from_rpre,
      NA_real_
    ),
    
    time_first_from_rpre = pmin_na(
      time_p2_from_rpre,
      time_death_from_rpre,
      time_censor_from_rpre
    ),
    
    first_event_from_rpre_ms = case_when(
      first_at_min(
        time_p2_from_rpre,
        time_first_from_rpre
      ) ~ "p2",
      first_at_min(
        time_death_from_rpre,
        time_first_from_rpre
      ) ~ "death",
      first_at_min(
        time_censor_from_rpre,
        time_first_from_rpre
      ) ~ "censor",
      TRUE ~ NA_character_
    )
  )

#tempi dopo revisione post-p2
analysis_dat <- analysis_dat %>%
  mutate(
    revision_post_p2_date = case_when(
      first_event_from_p2_ms == "revision" &
        first_revision_joint_from_p2 == "hip_revision" ~ hip_revision_date,
      first_event_from_p2_ms == "revision" &
        first_revision_joint_from_p2 == "knee_revision" ~ knee_revision_date,
      TRUE ~ as.Date(NA)
    ),
    
    time_death_from_rpost = years_between(
      death_date,
      revision_post_p2_date
    ),
    time_censor_from_rpost = years_between(
      censor_date,
      revision_post_p2_date
    ),
    
    time_death_from_rpost = if_else(
      time_death_from_rpost > 0,
      time_death_from_rpost,
      NA_real_
    ),
    time_censor_from_rpost = if_else(
      time_censor_from_rpost > 0,
      time_censor_from_rpost,
      NA_real_
    ),
    
    time_first_from_rpost = pmin_na(
      time_death_from_rpost,
      time_censor_from_rpost
    ),
    
    first_event_from_rpost_ms = case_when(
      first_at_min(
        time_death_from_rpost,
        time_first_from_rpost
      ) ~ "death",
      first_at_min(
        time_censor_from_rpost,
        time_first_from_rpost
      ) ~ "censor",
      TRUE ~ NA_character_
    )
  )

#controlli wide
p1_ms_wide <- analysis_dat %>%
  mutate(
    p1_first_event_check = first_event_from_p1_ms,
    p2_first_event_check = first_event_from_p2_ms,
    rpre_first_event_check = first_event_from_rpre_ms,
    rpost_first_event_check = first_event_from_rpost_ms
  )

qc_p1_first_event <- p1_ms_wide %>%
  count(
    order_group_ms,
    first_event_from_p1_ms
  ) %>%
  group_by(order_group_ms) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

qc_p2_first_event <- p1_ms_wide %>%
  filter(
    first_event_from_p1_ms == "p2"
  ) %>%
  count(
    order_group_ms,
    first_event_from_p2_ms
  ) %>%
  group_by(order_group_ms) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

qc_rpre_first_event <- p1_ms_wide %>%
  filter(
    first_event_from_p1_ms == "revision"
  ) %>%
  count(
    order_group_ms,
    first_event_from_rpre_ms
  ) %>%
  group_by(order_group_ms) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

qc_rpost_first_event <- p1_ms_wide %>%
  filter(
    first_event_from_p1_ms == "p2",
    first_event_from_p2_ms == "revision"
  ) %>%
  count(
    order_group_ms,
    first_event_from_rpost_ms
  ) %>%
  group_by(order_group_ms) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

#costruzione long
base_cols <- c(
  "patient_side_id",
  "CODPAT",
  "side",
  "hip",
  "hip_binary",
  "order_group_ms",
  "sex_ms",
  "age_p1",
  "bmi_p1",
  "bmi_p1_f_ms",
  "center_p1_label",
  "p1_year"
)

base_cols <- base_cols[
  base_cols %in% names(p1_ms_wide)
]

make_rows_p1 <- function(d) {
  d %>%
    transmute(
      across(all_of(base_cols)),
      from = "P1",
      to = "P2",
      trans = "P1_to_P2",
      Tstart = 0,
      Tstop = time_first_from_p1,
      status = as.integer(
        first_event_from_p1_ms == "p2"
      )
    ) %>%
    bind_rows(
      d %>%
        transmute(
          across(all_of(base_cols)),
          from = "P1",
          to = "R_pre",
          trans = "P1_to_Rpre",
          Tstart = 0,
          Tstop = time_first_from_p1,
          status = as.integer(
            first_event_from_p1_ms == "revision"
          )
        )
    ) %>%
    bind_rows(
      d %>%
        transmute(
          across(all_of(base_cols)),
          from = "P1",
          to = "D",
          trans = "P1_to_D",
          Tstart = 0,
          Tstop = time_first_from_p1,
          status = as.integer(
            first_event_from_p1_ms == "death"
          )
        )
    )
}

make_rows_p2 <- function(d) {
  d %>%
    filter(
      first_event_from_p1_ms == "p2"
    ) %>%
    transmute(
      across(all_of(base_cols)),
      from = "P2",
      to = "R_post",
      trans = "P2_to_Rpost",
      Tstart = 0,
      Tstop = time_first_from_p2,
      status = as.integer(
        first_event_from_p2_ms == "revision"
      )
    ) %>%
    bind_rows(
      d %>%
        filter(
          first_event_from_p1_ms == "p2"
        ) %>%
        transmute(
          across(all_of(base_cols)),
          from = "P2",
          to = "D",
          trans = "P2_to_D",
          Tstart = 0,
          Tstop = time_first_from_p2,
          status = as.integer(
            first_event_from_p2_ms == "death"
          )
        )
    )
}

make_rows_rpre <- function(d) {
  d %>%
    filter(
      first_event_from_p1_ms == "revision"
    ) %>%
    transmute(
      across(all_of(base_cols)),
      from = "R_pre",
      to = "P2_after_Rpre",
      trans = "Rpre_to_P2",
      Tstart = 0,
      Tstop = time_first_from_rpre,
      status = as.integer(
        first_event_from_rpre_ms == "p2"
      )
    ) %>%
    bind_rows(
      d %>%
        filter(
          first_event_from_p1_ms == "revision"
        ) %>%
        transmute(
          across(all_of(base_cols)),
          from = "R_pre",
          to = "D",
          trans = "Rpre_to_D",
          Tstart = 0,
          Tstop = time_first_from_rpre,
          status = as.integer(
            first_event_from_rpre_ms == "death"
          )
        )
    )
}

make_rows_rpost <- function(d) {
  d %>%
    filter(
      first_event_from_p1_ms == "p2",
      first_event_from_p2_ms == "revision"
    ) %>%
    transmute(
      across(all_of(base_cols)),
      from = "R_post",
      to = "D",
      trans = "Rpost_to_D",
      Tstart = 0,
      Tstop = time_first_from_rpost,
      status = as.integer(
        first_event_from_rpost_ms == "death"
      )
    )
}

p1_ms_long <- bind_rows(
  make_rows_p1(p1_ms_wide),
  make_rows_p2(p1_ms_wide),
  make_rows_rpre(p1_ms_wide),
  make_rows_rpost(p1_ms_wide)
) %>%
  filter(
    !is.na(Tstop),
    Tstop > 0
  ) %>%
  mutate(
    trans = factor(
      trans,
      levels = c(
        "P1_to_P2",
        "P1_to_Rpre",
        "P1_to_D",
        "P2_to_Rpost",
        "P2_to_D",
        "Rpre_to_P2",
        "Rpre_to_D",
        "Rpost_to_D"
      )
    ),
    from = factor(
      from,
      levels = c(
        "P1",
        "P2",
        "R_pre",
        "R_post"
      )
    ),
    to = factor(
      to,
      levels = c(
        "P2",
        "P2_after_Rpre",
        "R_pre",
        "R_post",
        "D"
      )
    )
  )

#conteggi transizioni
transition_counts <- p1_ms_long %>%
  group_by(
    trans,
    from,
    to
  ) %>%
  summarise(
    n_risk_rows = n(),
    n_events = sum(
      status == 1,
      na.rm = TRUE
    ),
    n_censored_or_competing = sum(
      status == 0,
      na.rm = TRUE
    ),
    event_rate_percent = 100 * n_events / n_risk_rows,
    .groups = "drop"
  ) %>%
  mutate(
    model_recommendation = case_when(
      n_events >= min_events_adjusted ~ "adjusted_model",
      n_events >= min_events_simple ~ "simple_or_descriptive",
      TRUE ~ "descriptive_only"
    )
  )

transition_counts_by_order <- p1_ms_long %>%
  group_by(
    trans,
    from,
    to,
    hip
  ) %>%
  summarise(
    n_risk_rows = n(),
    n_events = sum(
      status == 1,
      na.rm = TRUE
    ),
    n_censored_or_competing = sum(
      status == 0,
      na.rm = TRUE
    ),
    event_rate_percent = 100 * n_events / n_risk_rows,
    .groups = "drop"
  )

#matrice transizioni descrittiva
transition_matrix <- tibble(
  transition_id = 1:8,
  trans = c(
    "P1_to_P2",
    "P1_to_Rpre",
    "P1_to_D",
    "P2_to_Rpost",
    "P2_to_D",
    "Rpre_to_P2",
    "Rpre_to_D",
    "Rpost_to_D"
  ),
  from = c(
    "P1",
    "P1",
    "P1",
    "P2",
    "P2",
    "R_pre",
    "R_pre",
    "R_post"
  ),
  to = c(
    "P2",
    "R_pre",
    "D",
    "R_post",
    "D",
    "P2_after_Rpre",
    "D",
    "D"
  )
) %>%
  left_join(
    transition_counts %>%
      select(
        trans,
        n_risk_rows,
        n_events,
        model_recommendation
      ),
    by = "trans"
  )

#dataset modellabile
model_long <- p1_ms_long %>%
  mutate(
    age_p1_10 = age_p1 / 10,
    p1_year_c = p1_year - median(
      p1_year,
      na.rm = TRUE
    ),
    bmi_p1_f_ms = relevel(
      bmi_p1_f_ms,
      ref = "normal_underweight"
    ),
    sex_ms = relevel(
      sex_ms,
      ref = "male"
    )
  )

#transizioni per modelli
adjusted_transitions <- transition_counts %>%
  filter(
    model_recommendation == "adjusted_model"
  ) %>%
  pull(trans) %>%
  as.character()

simple_transitions <- transition_counts %>%
  filter(
    model_recommendation %in% c(
      "adjusted_model",
      "simple_or_descriptive"
    )
  ) %>%
  pull(trans) %>%
  as.character()

#modelli aggiustati
fit_adjusted_one <- function(tr) {
  d <- model_long %>%
    filter(
      trans == tr
    ) %>%
    filter(
      !is.na(Tstop),
      Tstop > 0,
      !is.na(status),
      !is.na(hip),
      !is.na(age_p1_10),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms)
    )
  
  n_ev <- sum(
    d$status == 1,
    na.rm = TRUE
  )
  
  if (n_ev < min_events_adjusted) {
    return(
      list(
        trans = tr,
        fit = structure(
          list(
            message = paste0(
              "not enough events for adjusted model: ",
              n_ev
            )
          ),
          class = "error"
        )
      )
    )
  }
  
  f <- Surv(Tstop, status) ~
    hip +
    age_p1_10 +
    sex_ms +
    bmi_p1_f_ms +
    cluster(CODPAT)
  
  list(
    trans = tr,
    fit = safe_cox(f, d)
  )
}

adjusted_fits <- map(
  adjusted_transitions,
  fit_adjusted_one
)

adjusted_results <- map_dfr(
  adjusted_fits,
  function(x) {
    tidy_cox_safe(
      x$fit,
      model_name = "adjusted",
      trans_label = x$trans
    )
  }
)

#modelli semplici solo hip
fit_simple_one <- function(tr) {
  d <- model_long %>%
    filter(
      trans == tr
    ) %>%
    filter(
      !is.na(Tstop),
      Tstop > 0,
      !is.na(status),
      !is.na(hip)
    )
  
  n_ev <- sum(
    d$status == 1,
    na.rm = TRUE
  )
  
  if (n_ev < min_events_simple) {
    return(
      list(
        trans = tr,
        fit = structure(
          list(
            message = paste0(
              "not enough events for simple model: ",
              n_ev