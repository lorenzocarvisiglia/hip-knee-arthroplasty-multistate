source("config.R")

# ============================================================
# 11_p1_year_center_sensitivity_updated.R
# sensitivity to calendar year and P1 center
# aligned with final file 04 specification
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(survival)
  library(openxlsx)
})

# ============================================================
# settings
# ============================================================

out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
final_models_file <- file.path(out_dir, "p1_multistate_final_models.RData")

if (!file.exists(analysis_file)) {
  stop("non trovo p1_multistate_analysis_objects.RData. lancia prima lo script 03.")
}

if (!file.exists(final_models_file)) {
  stop("non trovo p1_multistate_final_models.RData. lancia prima lo script 04 finale.")
}

min_center_n <- 100
times_report <- c(1, 3, 5, 10, 15)
tol <- 1e-8

adjusted_transitions <- c(
  "P1_to_P2",
  "P1_to_Rpre",
  "P1_to_D",
  "P2_to_Rpost",
  "P2_to_D",
  "Rpre_to_P2"
)

#final main hip(t) specification selected in file 04
final_tv_transitions <- c(
  "P1_to_P2",
  "P2_to_Rpost",
  "Rpre_to_P2"
)

model_types <- c(
  "base",
  "plus_p1_year",
  "plus_p1_year_center"
)

# ============================================================
# load objects
# ============================================================

analysis_env <- new.env()
load(analysis_file, envir = analysis_env)

if (!exists("p1_ms_long", envir = analysis_env)) {
  stop("p1_ms_long non trovato nel file 03.")
}

p1_ms_long <- analysis_env$p1_ms_long

final_env <- new.env()
load(final_models_file, envir = final_env)

if (!exists("final_main_hip_hr", envir = final_env)) {
  stop("final_main_hip_hr non trovato nel file 04 finale.")
}

final_main_hip_hr <- final_env$final_main_hip_hr

# ============================================================
# helpers
# ============================================================

standardize_order <- function(x) {
  x_chr <- str_to_lower(as.character(x))
  case_when(
    x_chr %in% c("1", "first_hip", "hip", "h", "anca", "first hip") ~ "first_hip",
    x_chr %in% c("0", "first_knee", "knee", "k", "ginocchio", "first knee") ~ "first_knee",
    TRUE ~ NA_character_
  )
}

standardize_sex <- function(x) {
  x_chr <- str_to_lower(as.character(x))
  case_when(
    x_chr %in% c("male", "m", "maschio", "1") ~ "male",
    x_chr %in% c("female", "f", "femmina", "0") ~ "female",
    TRUE ~ NA_character_
  )
}

standardize_bmi <- function(x) {
  x_chr <- str_to_lower(as.character(x))
  case_when(
    is.na(x) ~ "missing",
    str_detect(x_chr, "missing|unknown|na|non") ~ "missing",
    str_detect(x_chr, "obes") ~ "obese",
    str_detect(x_chr, "over|sovra") ~ "overweight",
    str_detect(x_chr, "normal|under|normo|sotto") ~ "normal_underweight",
    TRUE ~ "missing"
  )
}

extract_year_safe <- function(x) {
  if (inherits(x, "Date")) {
    return(as.integer(format(x, "%Y")))
  }
  
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) {
    return(as.integer(format(as.Date(x), "%Y")))
  }
  
  if (is.numeric(x)) {
    if (all(is.na(x) | (x > 1900 & x < 2100))) {
      return(as.integer(x))
    }
    
    out <- suppressWarnings(as.Date(x, origin = "1970-01-01"))
    return(as.integer(format(out, "%Y")))
  }
  
  out <- suppressWarnings(as.Date(x))
  as.integer(format(out, "%Y"))
}

clean_center <- function(x) {
  out <- str_squish(as.character(x))
  out[is.na(out) | out == ""] <- "Missing_center"
  out
}

tt_log <- function(x, t, ...) {
  x * log(t + 1)
}

safe_cox <- function(formula, data, hip_tv = FALSE) {
  warnings <- character()
  
  fit <- tryCatch(
    withCallingHandlers(
      coxph(
        formula,
        data = data,
        ties = "efron",
        x = TRUE,
        y = TRUE,
        model = !hip_tv,
        tt = if (hip_tv) tt_log else NULL,
        control = coxph.control(iter.max = 50)
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  list(
    fit = fit,
    warnings = unique(warnings)
  )
}

# ============================================================
# prepare common data
# ============================================================

dat <- p1_ms_long %>%
  mutate(
    hip_group = standardize_order(hip),
    hip = factor(hip_group, levels = c("first_knee", "first_hip")),
    hip_binary = case_when(
      hip_group == "first_hip" ~ 1,
      hip_group == "first_knee" ~ 0,
      TRUE ~ NA_real_
    ),
    age_cat = cut(
      age_p1,
      breaks = c(-Inf, 60, 70, 80, Inf),
      right = FALSE,
      labels = c("<60", "60-69", "70-79", "80+")
    )
  )

if ("sex_ms" %in% names(dat)) {
  dat <- dat %>% mutate(sex_std = standardize_sex(sex_ms))
} else if ("sex" %in% names(dat)) {
  dat <- dat %>% mutate(sex_std = standardize_sex(sex))
} else {
  stop("non trovo sex_ms né sex.")
}

dat <- dat %>%
  mutate(
    sex_std = factor(sex_std, levels = c("male", "female"))
  )

if ("bmi_p1_f_ms" %in% names(dat)) {
  dat <- dat %>% mutate(bmi_std = standardize_bmi(bmi_p1_f_ms))
} else if ("bmi_p1_f" %in% names(dat)) {
  dat <- dat %>% mutate(bmi_std = standardize_bmi(bmi_p1_f))
} else {
  stop("non trovo bmi_p1_f_ms né bmi_p1_f.")
}

dat <- dat %>%
  mutate(
    bmi_std = factor(
      bmi_std,
      levels = c("normal_underweight", "overweight", "obese", "missing")
    )
  )

#calendar year of P1
if ("p1_year" %in% names(dat)) {
  dat <- dat %>% mutate(p1_year = as.integer(p1_year))
} else if ("p1_date" %in% names(dat)) {
  dat <- dat %>% mutate(p1_year = extract_year_safe(p1_date))
} else {
  dat$p1_year <- NA_integer_
}

if (all(is.na(dat$p1_year))) {
  stop("p1_year non disponibile: impossibile eseguire la sensitivity per anno.")
}

median_p1_year <- median(dat$p1_year, na.rm = TRUE)

dat <- dat %>%
  mutate(
    p1_year_5 = (p1_year - median_p1_year) / 5
  )

#P1 center
center_available <- FALSE

if ("center_p1_label" %in% names(dat)) {
  center_available <- TRUE
  dat <- dat %>% mutate(center_p1_raw = clean_center(center_p1_label))
} else if ("center_p1" %in% names(dat)) {
  center_available <- TRUE
  dat <- dat %>% mutate(center_p1_raw = clean_center(center_p1))
} else {
  dat$center_p1_raw <- "Missing_center"
}

if (!center_available) {
  stop("center_p1 non disponibile: impossibile eseguire la sensitivity per centro.")
}

p1_base_units <- dat %>%
  filter(as.character(trans) == "P1_to_P2") %>%
  distinct(patient_side_id, .keep_all = TRUE)

center_lumping_table <- p1_base_units %>%
  count(center_p1_raw, name = "n_patient_sides") %>%
  arrange(desc(n_patient_sides)) %>%
  mutate(
    keep_as_separate_center =
      center_p1_raw != "Missing_center" &
      n_patient_sides >= min_center_n,
    center_p1_lumped = if_else(
      keep_as_separate_center,
      center_p1_raw,
      "Other_small_or_missing_centers"
    )
  )

centers_to_keep <- center_lumping_table %>%
  filter(keep_as_separate_center) %>%
  pull(center_p1_raw)

dat <- dat %>%
  mutate(
    center_p1_lumped = if_else(
      center_p1_raw %in% centers_to_keep,
      center_p1_raw,
      "Other_small_or_missing_centers"
    ),
    center_p1_lumped = factor(center_p1_lumped)
  )

# ============================================================
# model fitting
# ============================================================

prepare_transition <- function(tr, model_type) {
  d <- dat %>%
    filter(as.character(trans) == tr) %>%
    filter(
      !is.na(Tstop),
      Tstop > 0,
      !is.na(status),
      !is.na(CODPAT),
      !is.na(hip),
      !is.na(age_cat),
      !is.na(sex_std),
      !is.na(bmi_std)
    )
  
  if (model_type %in% c("plus_p1_year", "plus_p1_year_center")) {
    d <- d %>% filter(!is.na(p1_year_5))
  }
  
  if (model_type == "plus_p1_year_center") {
    d <- d %>% filter(!is.na(center_p1_lumped))
  }
  
  d
}

make_formula <- function(tr, model_type) {
  hip_tv <- tr %in% final_tv_transitions
  
  hip_terms <- if (hip_tv) {
    "hip_binary + tt(hip_binary)"
  } else {
    "hip"
  }
  
  rhs <- paste(
    hip_terms,
    "+ sex_std + bmi_std + strata(age_cat)"
  )
  
  if (model_type == "plus_p1_year") {
    rhs <- paste(rhs, "+ p1_year_5")
  }
  
  if (model_type == "plus_p1_year_center") {
    rhs <- paste(rhs, "+ p1_year_5 + strata(center_p1_lumped)")
  }
  
  as.formula(
    paste0(
      "Surv(Tstop, status) ~ ",
      rhs,
      " + cluster(CODPAT)"
    )
  )
}

extract_hip_effect <- function(fit, tr, model_type) {
  hip_tv <- tr %in% final_tv_transitions
  
  if (inherits(fit, "error")) {
    return(tibble(
      transition = tr,
      model_type = model_type,
      hip_time_varying = hip_tv,
      time_years = if (hip_tv) times_report else NA_real_,
      hr = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p = NA_real_,
      error = conditionMessage(fit)
    ))
  }
  
  b <- coef(fit)
  v <- vcov(fit)
  
  if (hip_tv) {
    needed <- c("hip_binary", "tt(hip_binary)")
    
    if (!all(needed %in% names(b))) {
      return(tibble(
        transition = tr,
        model_type = model_type,
        hip_time_varying = TRUE,
        time_years = times_report,
        hr = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p = NA_real_,
        error = "hip time-varying terms not found"
      ))
    }
    
    return(map_dfr(times_report, function(tt) {
      l <- log(tt + 1)
      est <- b["hip_binary"] + l * b["tt(hip_binary)"]
      var_est <-
        v["hip_binary", "hip_binary"] +
        l^2 * v["tt(hip_binary)", "tt(hip_binary)"] +
        2 * l * v["hip_binary", "tt(hip_binary)"]
      se <- sqrt(var_est)
      z <- est / se
      
      tibble(
        transition = tr,
        model_type = model_type,
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
      model_type = model_type,
      hip_time_varying = FALSE,
      time_years = NA_real_,
      hr = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p = NA_real_,
      error = "single hip PH term not found"
    ))
  }
  
  est <- b[hip_term]
  se <- sqrt(v[hip_term, hip_term])
  z <- est / se
  
  tibble(
    transition = tr,
    model_type = model_type,
    hip_time_varying = FALSE,
    time_years = NA_real_,
    hr = exp(est),
    conf_low = exp(est - 1.96 * se),
    conf_high = exp(est + 1.96 * se),
    p = 2 * pnorm(abs(z), lower.tail = FALSE),
    error = NA_character_
  )
}

fit_one <- function(tr, model_type) {
  d <- prepare_transition(tr, model_type)
  hip_tv <- tr %in% final_tv_transitions
  form <- make_formula(tr, model_type)
  res <- safe_cox(form, d, hip_tv = hip_tv)
  
  fit <- res$fit
  
  fit_status <- if (inherits(fit, "error")) {
    "FAILED"
  } else if (length(res$warnings) > 0) {
    "WARNING"
  } else {
    "OK"
  }
  
  fit_note <- if (inherits(fit, "error")) {
    conditionMessage(fit)
  } else if (length(res$warnings) > 0) {
    paste(res$warnings, collapse = " | ")
  } else {
    ""
  }
  
  info <- tibble(
    transition = tr,
    model_type = model_type,
    hip_time_varying = hip_tv,
    n_rows = nrow(d),
    n_patient_sides = n_distinct(d$patient_side_id),
    n_patients = n_distinct(d$CODPAT),
    n_events = sum(d$status == 1, na.rm = TRUE),
    n_first_knee = sum(d$hip == "first_knee", na.rm = TRUE),
    n_first_hip = sum(d$hip == "first_hip", na.rm = TRUE),
    events_first_knee = sum(d$status == 1 & d$hip == "first_knee", na.rm = TRUE),
    events_first_hip = sum(d$status == 1 & d$hip == "first_hip", na.rm = TRUE),
    status = fit_status,
    note = fit_note
  )
  
  hip_effect <- extract_hip_effect(fit, tr, model_type)
  
  all_coef <- if (inherits(fit, "error")) {
    tibble()
  } else {
    b <- coef(fit)
    v <- vcov(fit)
    se <- sqrt(diag(v))
    
    tibble(
      transition = tr,
      model_type = model_type,
      term = names(b),
      estimate = as.numeric(b),
      robust_se = as.numeric(se),
      hr = exp(as.numeric(b)),
      conf_low = exp(as.numeric(b) - 1.96 * as.numeric(se)),
      conf_high = exp(as.numeric(b) + 1.96 * as.numeric(se)),
      p = 2 * pnorm(abs(as.numeric(b) / as.numeric(se)), lower.tail = FALSE)
    )
  }
  
  list(
    fit = fit,
    data = d,
    info = info,
    hip_effect = hip_effect,
    all_coef = all_coef
  )
}

model_grid <- tidyr::crossing(
  transition = adjusted_transitions,
  model_type = factor(model_types, levels = model_types)
) %>%
  mutate(model_type = as.character(model_type))

fits <- pmap(
  model_grid,
  function(transition, model_type) {
    message("fit: ", transition, " | ", model_type)
    fit_one(transition, model_type)
  }
)

names(fits) <- paste(model_grid$transition, model_grid$model_type, sep = "__")

model_info <- map_dfr(fits, "info")
hip_effects <- map_dfr(fits, "hip_effect")
all_coefficients <- map_dfr(fits, "all_coef")

# ============================================================
# validate base model against file 04
# ============================================================

file04_reference <- final_main_hip_hr %>%
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
  filter(transition %in% adjusted_transitions) %>%
  distinct()

base_validation <- hip_effects %>%
  filter(model_type == "base") %>%
  select(
    transition,
    hip_time_varying,
    time_years,
    hr_refit = hr,
    low_refit = conf_low,
    high_refit = conf_high
  ) %>%
  left_join(
    file04_reference,
    by = c("transition", "hip_time_varying", "time_years")
  ) %>%
  mutate(
    abs_diff_hr = abs(hr_refit - hr_file04),
    match = !is.na(abs_diff_hr) & abs_diff_hr < tol
  )

if (!all(base_validation$match)) {
  stop("il modello base del file 11 non replica il file 04; controllare base_validation.")
}

# ============================================================
# compare hip effects across sensitivity models
# ============================================================

base_effects <- hip_effects %>%
  filter(model_type == "base") %>%
  select(
    transition,
    hip_time_varying,
    time_years,
    hr_base = hr
  )

hip_effect_comparison <- hip_effects %>%
  left_join(
    base_effects,
    by = c("transition", "hip_time_varying", "time_years")
  ) %>%
  mutate(
    hr_ratio_vs_base = hr / hr_base,
    percent_difference_vs_base = 100 * (hr / hr_base - 1),
    abs_log_hr_difference = abs(log(hr) - log(hr_base))
  ) %>%
  arrange(
    match(transition, adjusted_transitions),
    match(model_type, model_types),
    time_years
  )

# ============================================================
# descriptive year/center summaries
# ============================================================

year_by_order <- p1_base_units %>%
  filter(hip_group %in% c("first_knee", "first_hip")) %>%
  group_by(hip_group) %>%
  summarise(
    n_patient_sides = n_distinct(patient_side_id),
    min_p1_year = min(p1_year, na.rm = TRUE),
    q1_p1_year = quantile(p1_year, 0.25, na.rm = TRUE),
    median_p1_year = median(p1_year, na.rm = TRUE),
    mean_p1_year = mean(p1_year, na.rm = TRUE),
    q3_p1_year = quantile(p1_year, 0.75, na.rm = TRUE),
    max_p1_year = max(p1_year, na.rm = TRUE),
    .groups = "drop"
  )

center_by_order <- p1_base_units %>%
  mutate(
    center_p1_lumped = if_else(
      center_p1_raw %in% centers_to_keep,
      center_p1_raw,
      "Other_small_or_missing_centers"
    )
  ) %>%
  filter(hip_group %in% c("first_knee", "first_hip")) %>%
  count(center_p1_lumped, hip_group, name = "n_patient_sides") %>%
  group_by(hip_group) %>%
  mutate(
    pct_within_order = 100 * n_patient_sides / sum(n_patient_sides)
  ) %>%
  ungroup() %>%
  arrange(center_p1_lumped, hip_group)

center_age_strata_counts <- dat %>%
  filter(as.character(trans) %in% adjusted_transitions) %>%
  group_by(trans, center_p1_lumped, age_cat) %>%
  summarise(
    n_rows = n(),
    n_events = sum(status == 1, na.rm = TRUE),
    n_patients = n_distinct(CODPAT),
    .groups = "drop"
  ) %>%
  arrange(trans, center_p1_lumped, age_cat)

# ============================================================
# QC
# ============================================================

qc_sensitivity <- tibble(
  check = c(
    "base model reproduces file 04",
    "all 18 models fitted",
    "failed models",
    "models with warnings",
    "P1 year available",
    "P1 center available",
    "number of retained center strata",
    "median P1 year used for centering"
  ),
  value = c(
    as.character(all(base_validation$match)),
    as.character(nrow(model_info)),
    as.character(sum(model_info$status == "FAILED")),
    as.character(sum(model_info$status == "WARNING")),
    as.character(!all(is.na(dat$p1_year))),
    as.character(center_available),
    as.character(length(unique(dat$center_p1_lumped))),
    as.character(median_p1_year)
  ),
  status = c(
    ifelse(all(base_validation$match), "OK", "FAIL"),
    ifelse(nrow(model_info) == length(adjusted_transitions) * length(model_types), "OK", "FAIL"),
    ifelse(sum(model_info$status == "FAILED") == 0, "OK", "FAIL"),
    ifelse(sum(model_info$status == "WARNING") == 0, "OK", "CHECK"),
    "OK",
    "OK",
    "INFO",
    "INFO"
  )
)

# ============================================================
# save outputs
# ============================================================

out_xlsx <- file.path(out_dir, "p1_year_center_sensitivity_updated.xlsx")
out_rdata <- file.path(out_dir, "p1_year_center_sensitivity_updated.RData")

notes <- tibble(
  item = c(
    "base specification",
    "calendar year",
    "center",
    "time-varying order effects"
  ),
  description = c(
    "The base model exactly reproduces the final file 04 specification: age strata, sex, BMI missing category, patient clustering, Efron ties.",
    paste0("P1 calendar year is entered linearly per 5 years and centered at ", median_p1_year, "."),
    paste0("Centers with at least ", min_center_n, " P1 patient-sides are retained; smaller or missing centers are combined. Center is added as a Cox stratum together with age strata."),
    "hip(t) is retained only for P1_to_P2, P2_to_Rpost, and Rpre_to_P2, exactly as in the final file 04."
  )
)

wb <- createWorkbook()

tabs <- list(
  notes = notes,
  qc = qc_sensitivity,
  base_validation = base_validation,
  hip_effect_comparison = hip_effect_comparison,
  model_info = model_info,
  all_coefficients = all_coefficients,
  year_by_order = year_by_order,
  center_lumping = center_lumping_table,
  center_by_order = center_by_order,
  center_age_strata = center_age_strata_counts
)

for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tabs[[nm]])
  setColWidths(wb, nm, cols = 1:100, widths = "auto")
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

save(
  dat,
  fits,
  model_grid,
  model_info,
  hip_effects,
  hip_effect_comparison,
  all_coefficients,
  base_validation,
  qc_sensitivity,
  year_by_order,
  center_lumping_table,
  center_by_order,
  center_age_strata_counts,
  adjusted_transitions,
  final_tv_transitions,
  model_types,
  min_center_n,
  median_p1_year,
  file = out_rdata
)

cat("\n==============================\n")
cat("QC YEAR/CENTER SENSITIVITY\n")
cat("==============================\n")
print(qc_sensitivity)

cat("\n==============================\n")
cat("BASE VALIDATION VS FILE 04\n")
cat("==============================\n")
print(base_validation)

cat("\n==============================\n")
cat("HIP EFFECT COMPARISON\n")
cat("==============================\n")
print(hip_effect_comparison)

cat("\noutput salvati in:\n")
cat(out_xlsx, "\n")
cat(out_rdata, "\n")