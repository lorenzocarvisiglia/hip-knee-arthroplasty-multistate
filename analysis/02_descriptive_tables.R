source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(survival)
  library(openxlsx)
})

#parametri
out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(clean_data_file)) {
  stop("Missing cleaned data file: ", clean_data_file, ". Run analysis/01_clean_data.R first.")
}
load(clean_data_file)

dat0 <- cohort_p1_clean

#funzioni
to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x_chr <- str_squish(as.character(x))
  x_chr[x_chr %in% c("", "-", "NA", "NaN", "NULL")] <- NA_character_
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

fmt_n_pct <- function(n, den) {
  ifelse(
    is.na(den) | den == 0,
    paste0(n, " (NA)"),
    paste0(n, " (", sprintf("%.1f", 100 * n / den), "%)")
  )
}

q25 <- function(x) quantile(x, 0.25, na.rm = TRUE, names = FALSE)
q75 <- function(x) quantile(x, 0.75, na.rm = TRUE, names = FALSE)

iqr_text <- function(x) {
  paste0(sprintf("%.2f", q25(x)), " to ", sprintf("%.2f", q75(x)))
}

make_km_estimates <- function(dat, time_col, event_col, group_col, endpoint_label, times = c(5, 10, 15)) {
  d <- dat %>%
    filter(
      !is.na(.data[[time_col]]),
      .data[[time_col]] > 0,
      !is.na(.data[[event_col]]),
      !is.na(.data[[group_col]])
    )
  
  if (nrow(d) == 0 || sum(d[[event_col]], na.rm = TRUE) == 0) {
    return(tibble(endpoint = endpoint_label, note = "no events or no data"))
  }
  
  fit <- survfit(
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", group_col)),
    data = d
  )
  
  ss <- summary(fit, times = times, extend = TRUE)
  
  tibble(
    endpoint = endpoint_label,
    group = sub(paste0(group_col, "="), "", ss$strata),
    time_years = ss$time,
    survival = ss$surv,
    event_probability = 1 - ss$surv,
    lower_event_probability = 1 - ss$upper,
    upper_event_probability = 1 - ss$lower
  )
}

make_logrank <- function(dat, time_col, event_col, group_col, endpoint_label) {
  d <- dat %>%
    filter(
      !is.na(.data[[time_col]]),
      .data[[time_col]] > 0,
      !is.na(.data[[event_col]]),
      !is.na(.data[[group_col]])
    )
  
  if (nrow(d) == 0 || sum(d[[event_col]], na.rm = TRUE) == 0) {
    return(tibble(endpoint = endpoint_label, chisq = NA_real_, df = NA_integer_, p_value = NA_real_))
  }
  
  lr <- survdiff(
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", group_col)),
    data = d
  )
  
  df <- length(lr$n) - 1
  
  tibble(
    endpoint = endpoint_label,
    chisq = unname(lr$chisq),
    df = df,
    p_value = pchisq(lr$chisq, df = df, lower.tail = FALSE)
  )
}

#preparazione dati
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
    order_group_clean = case_when(
      hip == 1 ~ "First hip",
      hip == 0 ~ "First knee",
      TRUE ~ "Undefined"
    ),
    order_group_clean = factor(order_group_clean, levels = c("First hip", "First knee", "Undefined")),
    sex_clean = case_when(
      sex %in% c("F", "female", "Female", "1") ~ "Female",
      sex %in% c("M", "male", "Male", "0") ~ "Male",
      TRUE ~ as.character(sex)
    ),
    bmi_p1_f_clean = case_when(
      bmi_p1_f %in% c("normal_underweight", "normal", "normal/underweight", "Sottopeso o Normopeso") ~ "Normal or underweight BMI",
      bmi_p1_f %in% c("overweight", "Sovrappeso") ~ "Overweight BMI",
      bmi_p1_f %in% c("obese", "Obeso") ~ "Obese BMI",
      is.na(bmi_p1_f) ~ "Missing BMI",
      TRUE ~ as.character(bmi_p1_f)
    ),
    bmi_p1_f_clean = factor(
      bmi_p1_f_clean,
      levels = c("Normal or underweight BMI", "Overweight BMI", "Obese BMI", "Missing BMI")
    ),
    time_censor_from_p1 = years_between(censor_date, p1_date),
    time_p2_from_p1 = years_between(p2_date, p1_date),
    time_death_from_p1 = years_between(death_date, p1_date),
    time_hip_revision_from_p1 = years_between(hip_revision_date, p1_date),
    time_knee_revision_from_p1 = years_between(knee_revision_date, p1_date),
    time_hip_revision_from_p1 = if_else(time_hip_revision_from_p1 <= 0, NA_real_, time_hip_revision_from_p1),
    time_knee_revision_from_p1 = if_else(time_knee_revision_from_p1 <= 0, NA_real_, time_knee_revision_from_p1),
    time_death_from_p1 = if_else(time_death_from_p1 <= 0, NA_real_, time_death_from_p1),
    time_p2_from_p1 = if_else(time_p2_from_p1 <= 0, NA_real_, time_p2_from_p1),
    time_censor_from_p1 = if_else(time_censor_from_p1 <= 0, NA_real_, time_censor_from_p1),
    time_first_revision_from_p1 = pmin_na(time_hip_revision_from_p1, time_knee_revision_from_p1),
    first_revision_joint_from_p1 = case_when(
      !is.na(time_hip_revision_from_p1) &
        (is.na(time_knee_revision_from_p1) | time_hip_revision_from_p1 <= time_knee_revision_from_p1) ~ "Hip revision",
      !is.na(time_knee_revision_from_p1) &
        (is.na(time_hip_revision_from_p1) | time_knee_revision_from_p1 < time_hip_revision_from_p1) ~ "Knee revision",
      TRUE ~ NA_character_
    ),
    time_first_from_p1_recalc = pmin_na(
      time_p2_from_p1,
      time_first_revision_from_p1,
      time_death_from_p1,
      time_censor_from_p1
    ),
    first_event_from_p1_recalc = case_when(
      !is.na(time_p2_from_p1) &
        abs(time_p2_from_p1 - time_first_from_p1_recalc) < 1e-8 ~ "p2",
      !is.na(time_first_revision_from_p1) &
        abs(time_first_revision_from_p1 - time_first_from_p1_recalc) < 1e-8 ~ "revision",
      !is.na(time_death_from_p1) &
        abs(time_death_from_p1 - time_first_from_p1_recalc) < 1e-8 ~ "death",
      !is.na(time_censor_from_p1) &
        abs(time_censor_from_p1 - time_first_from_p1_recalc) < 1e-8 ~ "censor",
      TRUE ~ NA_character_
    ),
    has_p2_by_admin = as.logical(has_p2_by_admin),
    any_revision_after_p1 = !is.na(time_first_revision_from_p1),
    any_death_within_censoring = !is.na(time_death_from_p1) &
      !is.na(time_censor_from_p1) &
      time_death_from_p1 <= time_censor_from_p1
  )

#variabili post-p2 descrittive
dat <- dat %>%
  mutate(
    time_hip_revision_from_p2 = years_between(hip_revision_date, p2_date),
    time_knee_revision_from_p2 = years_between(knee_revision_date, p2_date),
    time_death_from_p2 = years_between(death_date, p2_date),
    time_censor_from_p2 = years_between(censor_date, p2_date),
    time_hip_revision_from_p2 = if_else(has_p2_by_admin & time_hip_revision_from_p2 > 0, time_hip_revision_from_p2, NA_real_),
    time_knee_revision_from_p2 = if_else(has_p2_by_admin & time_knee_revision_from_p2 > 0, time_knee_revision_from_p2, NA_real_),
    time_death_from_p2 = if_else(has_p2_by_admin & time_death_from_p2 > 0, time_death_from_p2, NA_real_),
    time_censor_from_p2 = if_else(has_p2_by_admin & time_censor_from_p2 > 0, time_censor_from_p2, NA_real_),
    time_first_from_p2 = pmin_na(
      time_hip_revision_from_p2,
      time_knee_revision_from_p2,
      time_death_from_p2,
      time_censor_from_p2
    ),
    first_event_from_p2 = case_when(
      !is.na(time_knee_revision_from_p2) &
        abs(time_knee_revision_from_p2 - time_first_from_p2) < 1e-8 ~ "first knee revision after P2",
      !is.na(time_hip_revision_from_p2) &
        abs(time_hip_revision_from_p2 - time_first_from_p2) < 1e-8 ~ "first hip revision after P2",
      !is.na(time_death_from_p2) &
        abs(time_death_from_p2 - time_first_from_p2) < 1e-8 ~ "death after P2 before revision",
      !is.na(time_censor_from_p2) &
        abs(time_censor_from_p2 - time_first_from_p2) < 1e-8 ~ "censor after P2",
      TRUE ~ NA_character_
    )
  )

analysis_dat <- dat %>%
  filter(!is.na(hip))

#tabella 1: riepilogo generale
general_qc <- tibble(
  Quantity = c(
    "Patient-side units",
    "Patients",
    "Patients contributing both sides",
    "First hip group",
    "First knee group",
    "Undefined order coding",
    "Patient-side units reaching P2 within censoring",
    "Patient-side units not reaching P2 within censoring",
    "Any first revision after P1",
    "First revision before P2 or without P2",
    "Deaths before P2 or without P2",
    "Deaths within censoring",
    "Median follow-up from P1, years",
    "IQR follow-up from P1, years",
    "Maximum follow-up from P1, years"
  ),
  Value = c(
    nrow(dat),
    n_distinct(dat$CODPAT),
    sum(table(dat$CODPAT) > 1),
    sum(dat$hip == 1, na.rm = TRUE),
    sum(dat$hip == 0, na.rm = TRUE),
    sum(is.na(dat$hip)),
    sum(dat$has_p2_by_admin, na.rm = TRUE),
    sum(!dat$has_p2_by_admin, na.rm = TRUE),
    sum(dat$any_revision_after_p1, na.rm = TRUE),
    sum(dat$first_event_from_p1_recalc == "revision", na.rm = TRUE),
    sum(dat$first_event_from_p1_recalc == "death", na.rm = TRUE),
    sum(dat$any_death_within_censoring, na.rm = TRUE),
    sprintf("%.2f", median(dat$time_censor_from_p1, na.rm = TRUE)),
    iqr_text(dat$time_censor_from_p1),
    sprintf("%.2f", max(dat$time_censor_from_p1, na.rm = TRUE))
  )
)

#tabella 2: popolazione analisi p1
p1_population <- tibble(
  Quantity = c(
    "P1 patient-side units with defined order",
    "Patients",
    "Patients contributing both sides",
    "First hip group",
    "First knee group",
    "P2 as first observed transition from P1",
    "Revision as first observed transition from P1",
    "Death as first observed transition from P1",
    "Censoring as first observed transition from P1",
    "Median follow-up from P1, years",
    "IQR follow-up from P1, years",
    "Maximum follow-up from P1, years"
  ),
  Value = c(
    nrow(analysis_dat),
    n_distinct(analysis_dat$CODPAT),
    sum(table(analysis_dat$CODPAT) > 1),
    sum(analysis_dat$hip == 1, na.rm = TRUE),
    sum(analysis_dat$hip == 0, na.rm = TRUE),
    sum(analysis_dat$first_event_from_p1_recalc == "p2", na.rm = TRUE),
    sum(analysis_dat$first_event_from_p1_recalc == "revision", na.rm = TRUE),
    sum(analysis_dat$first_event_from_p1_recalc == "death", na.rm = TRUE),
    sum(analysis_dat$first_event_from_p1_recalc == "censor", na.rm = TRUE),
    sprintf("%.2f", median(analysis_dat$time_censor_from_p1, na.rm = TRUE)),
    iqr_text(analysis_dat$time_censor_from_p1),
    sprintf("%.2f", max(analysis_dat$time_censor_from_p1, na.rm = TRUE))
  )
)

#tabella 3: continui per ordine
continuous_by_order <- analysis_dat %>%
  group_by(order_group_clean) %>%
  summarise(
    `Patient-side units` = n(),
    Patients = n_distinct(CODPAT),
    `Mean age at P1, years` = mean(age_p1, na.rm = TRUE),
    `SD age at P1, years` = sd(age_p1, na.rm = TRUE),
    `Median age at P1, years` = median(age_p1, na.rm = TRUE),
    `Median t12 among those reaching P2, years` = median(t12_by_admin, na.rm = TRUE),
    `IQR t12 among those reaching P2, years` = iqr_text(t12_by_admin),
    `Median follow-up from P1, years` = median(time_censor_from_p1, na.rm = TRUE),
    `IQR follow-up from P1, years` = iqr_text(time_censor_from_p1),
    .groups = "drop"
  )

#tabella 4: sesso e bmi per ordine
den_by_order <- analysis_dat %>%
  count(order_group_clean, name = "den")

sex_bmi_by_order <- bind_rows(
  analysis_dat %>%
    count(order_group_clean, category = sex_clean, name = "n") %>%
    left_join(den_by_order, by = "order_group_clean") %>%
    mutate(section = "Sex", value = fmt_n_pct(n, den)) %>%
    select(section, category, order_group_clean, value),
  analysis_dat %>%
    count(order_group_clean, category = bmi_p1_f_clean, name = "n") %>%
    left_join(den_by_order, by = "order_group_clean") %>%
    mutate(section = "BMI at P1", value = fmt_n_pct(n, den)) %>%
    select(section, category, order_group_clean, value)
) %>%
  pivot_wider(names_from = order_group_clean, values_from = value) %>%
  arrange(section, category)

#tabella 5: primo evento da p1 per ordine
first_event_by_order <- analysis_dat %>%
  count(order_group_clean, first_event_from_p1_recalc, name = "n") %>%
  group_by(order_group_clean) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(value = paste0(n, " (", sprintf("%.1f", percent), "%)")) %>%
  select(first_event_from_p1_recalc, order_group_clean, value) %>%
  pivot_wider(names_from = order_group_clean, values_from = value)

#tabella 6: eventi osservati per ordine
events_by_order <- analysis_dat %>%
  group_by(order_group_clean) %>%
  summarise(
    `Patient-side units` = n(),
    `P1 to P2 as first event` = sum(first_event_from_p1_recalc == "p2", na.rm = TRUE),
    `P1 to first revision as first event` = sum(first_event_from_p1_recalc == "revision", na.rm = TRUE),
    `P1 to death as first event` = sum(first_event_from_p1_recalc == "death", na.rm = TRUE),
    `P1 to censoring as first event` = sum(first_event_from_p1_recalc == "censor", na.rm = TRUE),
    `Reached P2 within censoring` = sum(has_p2_by_admin, na.rm = TRUE),
    `Did not reach P2 within censoring` = sum(!has_p2_by_admin, na.rm = TRUE),
    `First hip revision after P1` = sum(first_revision_joint_from_p1 == "Hip revision", na.rm = TRUE),
    `First knee revision after P1` = sum(first_revision_joint_from_p1 == "Knee revision", na.rm = TRUE),
    `First hip revision after P2` = sum(first_event_from_p2 == "first hip revision after P2", na.rm = TRUE),
    `First knee revision after P2` = sum(first_event_from_p2 == "first knee revision after P2", na.rm = TRUE),
    `Death after P2 before revision` = sum(first_event_from_p2 == "death after P2 before revision", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-order_group_clean, names_to = "Event", values_to = "n") %>%
  pivot_wider(names_from = order_group_clean, values_from = n)

#tabella 7: subset con p2
p2_subset_summary <- analysis_dat %>%
  filter(has_p2_by_admin) %>%
  group_by(order_group_clean) %>%
  summarise(
    `Patient-side units reaching P2` = n(),
    Patients = n_distinct(CODPAT),
    `Mean age at P2, years` = mean(age_p2, na.rm = TRUE),
    `SD age at P2, years` = sd(age_p2, na.rm = TRUE),
    `Median t12, years` = median(t12_by_admin, na.rm = TRUE),
    `IQR t12, years` = iqr_text(t12_by_admin),
    `Median follow-up from P2, years` = median(time_censor_from_p2, na.rm = TRUE),
    `IQR follow-up from P2, years` = iqr_text(time_censor_from_p2),
    .groups = "drop"
  )

#tabella 8: missingness
missingness <- analysis_dat %>%
  group_by(order_group_clean) %>%
  summarise(
    n = n(),
    missing_age_p1 = sum(is.na(age_p1)),
    missing_bmi_p1 = sum(is.na(bmi_p1)),
    missing_bmi_p1_f = sum(is.na(bmi_p1_f_clean)),
    missing_sex = sum(is.na(sex_clean)),
    missing_center_p1 = sum(is.na(center_p1_label)),
    missing_knee_prosthesis_type = sum(is.na(knee_prosthesis_label)),
    missing_p1_date = sum(is.na(p1_date)),
    missing_censor_date = sum(is.na(censor_date)),
    .groups = "drop"
  ) %>%
  mutate(
    pct_missing_age_p1 = 100 * missing_age_p1 / n,
    pct_missing_bmi_p1 = 100 * missing_bmi_p1 / n,
    pct_missing_sex = 100 * missing_sex / n,
    pct_missing_center_p1 = 100 * missing_center_p1 / n,
    pct_missing_knee_prosthesis_type = 100 * missing_knee_prosthesis_type / n
  )

#tabella 9: centri p1
center_p1 <- analysis_dat %>%
  count(order_group_clean, center_p1_label, sort = TRUE) %>%
  group_by(order_group_clean) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

#tabella 10: tipo protesi ginocchio
knee_type <- analysis_dat %>%
  count(order_group_clean, knee_prosthesis_label, sort = TRUE) %>%
  group_by(order_group_clean) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

#km esplorative da p1
km_dat <- analysis_dat %>%
  mutate(
    time_any_revision_p1 = pmin_na(time_first_revision_from_p1, time_death_from_p1, time_censor_from_p1),
    event_any_revision_p1 = as.integer(
      !is.na(time_first_revision_from_p1) &
        abs(time_first_revision_from_p1 - time_any_revision_p1) < 1e-8
    ),
    time_p2_p1 = pmin_na(time_p2_from_p1, time_first_revision_from_p1, time_death_from_p1, time_censor_from_p1),
    event_p2_p1 = as.integer(
      !is.na(time_p2_from_p1) &
        abs(time_p2_from_p1 - time_p2_p1) < 1e-8
    ),
    time_death_p1 = pmin_na(time_death_from_p1, time_censor_from_p1),
    event_death_p1 = as.integer(
      !is.na(time_death_from_p1) &
        abs(time_death_from_p1 - time_death_p1) < 1e-8
    )
  )

km_logrank <- bind_rows(
  make_logrank(km_dat, "time_any_revision_p1", "event_any_revision_p1", "order_group_clean", "Any first revision after P1"),
  make_logrank(km_dat, "time_p2_p1", "event_p2_p1", "order_group_clean", "Reaching P2 after P1"),
  make_logrank(km_dat, "time_death_p1", "event_death_p1", "order_group_clean", "Death after P1")
)

km_estimates <- bind_rows(
  make_km_estimates(km_dat, "time_any_revision_p1", "event_any_revision_p1", "order_group_clean", "Any first revision after P1"),
  make_km_estimates(km_dat, "time_p2_p1", "event_p2_p1", "order_group_clean", "Reaching P2 after P1"),
  make_km_estimates(km_dat, "time_death_p1", "event_death_p1", "order_group_clean", "Death after P1")
)

#salvataggio excel
wb <- createWorkbook()

tabs <- list(
  general_qc = general_qc,
  p1_population = p1_population,
  continuous_by_order = continuous_by_order,
  sex_bmi_by_order = sex_bmi_by_order,
  first_event_by_order = first_event_by_order,
  events_by_order = events_by_order,
  p2_subset_summary = p2_subset_summary,
  missingness = missingness,
  center_p1 = center_p1,
  knee_type = knee_type,
  km_logrank = km_logrank,
  km_estimates = km_estimates
)

for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tabs[[nm]])
  setColWidths(wb, nm, cols = 1:80, widths = "auto")
}

saveWorkbook(
  wb,
  file.path(out_dir, "p1_descriptive_tables_like_postp2.xlsx"),
  overwrite = TRUE
)

#salvataggio oggetti
save(
  dat,
  analysis_dat,
  general_qc,
  p1_population,
  continuous_by_order,
  sex_bmi_by_order,
  first_event_by_order,
  events_by_order,
  p2_subset_summary,
  missingness,
  center_p1,
  knee_type,
  km_logrank,
  km_estimates,
  file = file.path(out_dir, "p1_descriptive_tables_like_postp2.RData")
)

print(general_qc)
print(p1_population)
print(continuous_by_order)
print(first_event_by_order)

cat(
  "\noutput salvati in:\n",
  file.path(out_dir, "p1_descriptive_tables_like_postp2.xlsx"),
  "\n",
  file.path(out_dir, "p1_descriptive_tables_like_postp2.RData"),
  "\n",
  sep = ""
)