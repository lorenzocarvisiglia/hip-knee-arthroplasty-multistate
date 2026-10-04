source("config.R")

#implant-specific secondary analysis using time since each joint implantation

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(survival)
  library(openxlsx)
  library(ggplot2)
  library(scales)
})

out_dir <- output_dir
fig_dir <- file.path(out_dir, "figures_p1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
clean_candidates <- c(
  file.path(out_dir, "p1_multistate_option3_clean_data.RData"),
  clean_data_file,
  file.path("Excel Data", "analisi_P1", "p1_multistate_option3_clean_data.RData"),
  "p1_multistate_option3_clean_data.RData"
)
out_xlsx <- file.path(out_dir, "p1_implant_specific_revision_secondary.xlsx")
out_rdata <- file.path(out_dir, "p1_implant_specific_revision_secondary.RData")
out_plot <- file.path(fig_dir, "p1_implant_specific_revision_secondary.png")

report_times <- c(1, 3, 5, 10, 15)
curve_times <- seq(0, 15, by = 0.25)
tol <- 1e-10

load_cohort <- function() {
  candidates <- unique(c(analysis_file, clean_candidates))
  candidates <- candidates[file.exists(candidates)]
  
  for (f in candidates) {
    e <- new.env(parent = globalenv())
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
    if (!ok) next
    
    if (exists("cohort_p1_clean", envir = e, inherits = FALSE) &&
        is.data.frame(e$cohort_p1_clean)) {
      dat <- as_tibble(e$cohort_p1_clean)
      message("cohort_p1_clean loaded from: ", normalizePath(f, winslash = "/", mustWork = FALSE))
      return(list(data = dat, source = f))
    }
  }
  
  stop(
    "cohort_p1_clean not found. Expected the clean file produced by script 01, usually at: ",
    file.path("output_rizzoli", "full", "p1_multistate_option3_cleaning", "p1_multistate_option3_clean_data.RData")
  )
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  suppressWarnings(as.Date(x))
}

years_between <- function(end, start) {
  as.numeric(as_date_safe(end) - as_date_safe(start)) / 365.25
}

standardize_sex <- function(x) {
  z <- str_to_lower(as.character(x))
  case_when(
    z %in% c("male", "m", "maschio", "1") ~ "male",
    z %in% c("female", "f", "femmina", "0", "2") ~ "female",
    TRUE ~ NA_character_
  )
}

standardize_bmi <- function(x) {
  z <- str_to_lower(as.character(x))
  case_when(
    is.na(x) ~ "missing",
    str_detect(z, "missing|unknown|non|^na$") ~ "missing",
    str_detect(z, "obes") ~ "obese",
    str_detect(z, "over|sovra") ~ "overweight",
    str_detect(z, "normal|under|normo|sotto") ~ "normal_underweight",
    TRUE ~ "missing"
  )
}

pmin_date <- function(...) {
  xs <- list(...)
  m <- do.call(cbind, lapply(xs, function(x) as.numeric(as_date_safe(x))))
  z <- apply(m, 1, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
  as.Date(z, origin = "1970-01-01")
}

aj_at_times <- function(time, status, eval_times) {
  d <- tibble(time = as.numeric(time), status = as.integer(status)) %>%
    filter(is.finite(time), time >= 0, status %in% 0:2)
  
  if (nrow(d) == 0) {
    return(tibble(time_years = eval_times, revision = NA_real_, death = NA_real_, survival = NA_real_))
  }
  
  ev_times <- sort(unique(d$time[d$status %in% c(1L, 2L)]))
  ev_times <- ev_times[ev_times <= max(eval_times)]
  S <- 1
  F_rev <- 0
  F_death <- 0
  step <- tibble(time = 0, revision = 0, death = 0, survival = 1)
  
  for (u in ev_times) {
    nr <- sum(d$time >= u)
    if (nr <= 0) next
    d1 <- sum(abs(d$time - u) < tol & d$status == 1L)
    d2 <- sum(abs(d$time - u) < tol & d$status == 2L)
    S0 <- S
    F_rev <- F_rev + S0 * d1 / nr
    F_death <- F_death + S0 * d2 / nr
    S <- S0 * (1 - (d1 + d2) / nr)
    step <- bind_rows(step, tibble(time = u, revision = F_rev, death = F_death, survival = S))
  }
  
  map_dfr(eval_times, function(tt) {
    idx <- max(which(step$time <= tt))
    step[idx, ] %>% transmute(time_years = tt, revision, death, survival)
  })
}

fit_revision_cox <- function(d) {
  dd <- d %>%
    mutate(
      event_revision = as.integer(status == 1L),
      implant_first = as.integer(implant_order == "first"),
      age_at_implant_10 = age_at_implant / 10,
      sex_f = factor(sex_std, levels = c("male", "female")),
      bmi_f = factor(bmi_std, levels = c("normal_underweight", "overweight", "obese", "missing"))
    ) %>%
    filter(is.finite(time_years), time_years > 0, !is.na(implant_first), !is.na(CODPAT))
  
  rhs <- c("implant_first")
  if (sum(is.finite(dd$age_at_implant_10)) > 20 && sd(dd$age_at_implant_10, na.rm = TRUE) > 0) rhs <- c(rhs, "age_at_implant_10")
  if (n_distinct(dd$sex_f[!is.na(dd$sex_f)]) >= 2) rhs <- c(rhs, "sex_f")
  if (n_distinct(dd$bmi_f[!is.na(dd$bmi_f)]) >= 2) rhs <- c(rhs, "bmi_f")
  if (sum(is.finite(dd$implant_year)) > 20 && sd(dd$implant_year, na.rm = TRUE) > 0) rhs <- c(rhs, "implant_year")
  rhs <- c(rhs, "cluster(CODPAT)")
  
  f <- as.formula(paste0("Surv(time_years, event_revision) ~ ", paste(rhs, collapse = " + ")))
  fit <- tryCatch(coxph(f, data = dd, ties = "efron", x = TRUE, y = TRUE), error = function(e) e)
  list(fit = fit, data = dd, formula = paste(deparse(f), collapse = ""))
}

extract_cox <- function(obj, population, joint) {
  fit <- obj$fit
  if (inherits(fit, "error")) {
    return(tibble(analysis_population = population, joint = joint, term = "error", hr = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_, n = nrow(obj$data), events = sum(obj$data$event_revision), formula = obj$formula, note = conditionMessage(fit)))
  }
  s <- summary(fit)
  co <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  pcol <- grep("Pr", names(co), value = TRUE)[1]
  tibble(
    analysis_population = population,
    joint = joint,
    term = rownames(co),
    hr = ci[["exp(coef)"]],
    lower_95 = ci[["lower .95"]],
    upper_95 = ci[["upper .95"]],
    p_value = co[[pcol]],
    n = fit$n,
    events = fit$nevent,
    formula = obj$formula,
    note = NA_character_
  )
}

extract_ph <- function(obj, population, joint) {
  fit <- obj$fit
  if (inherits(fit, "error")) return(tibble(analysis_population = population, joint = joint, term = "implant_first", chisq = NA_real_, p_value = NA_real_, note = "model error"))
  z <- tryCatch(cox.zph(fit), error = function(e) e)
  if (inherits(z, "error")) return(tibble(analysis_population = population, joint = joint, term = "implant_first", chisq = NA_real_, p_value = NA_real_, note = conditionMessage(z)))
  tab <- as.data.frame(z$table)
  if (!"implant_first" %in% rownames(tab)) return(tibble(analysis_population = population, joint = joint, term = "implant_first", chisq = NA_real_, p_value = NA_real_, note = "term not found"))
  tibble(analysis_population = population, joint = joint, term = "implant_first", chisq = tab["implant_first", "chisq"], p_value = tab["implant_first", "p"], note = NA_character_)
}

src <- load_cohort()
cohort <- src$data
needed <- c("patient_side_id", "CODPAT", "p1_date", "p1_joint", "p2_date", "p2_joint", "hip_primary_date", "knee_primary_date", "hip_revision_date", "knee_revision_date", "death_date", "censor_date", "age_p1", "sex", "bmi_p1_f")
missing_needed <- setdiff(needed, names(cohort))
if (length(missing_needed) > 0) stop("missing columns in cohort_p1_clean: ", paste(missing_needed, collapse = ", "))

base <- cohort %>%
  transmute(
    patient_side_id = as.character(patient_side_id),
    CODPAT = as.character(CODPAT),
    p1_date = as_date_safe(p1_date),
    p1_joint = as.character(p1_joint),
    p2_date = as_date_safe(p2_date),
    p2_joint = as.character(p2_joint),
    hip_primary_date = as_date_safe(hip_primary_date),
    knee_primary_date = as_date_safe(knee_primary_date),
    hip_revision_date = as_date_safe(hip_revision_date),
    knee_revision_date = as_date_safe(knee_revision_date),
    death_date = as_date_safe(death_date),
    censor_date = as_date_safe(censor_date),
    age_p1 = as.numeric(age_p1),
    sex_std = standardize_sex(sex),
    bmi_std = standardize_bmi(bmi_p1_f)
  ) %>%
  filter(!is.na(patient_side_id), !is.na(CODPAT), !is.na(p1_date), p1_joint %in% c("hip", "knee"))

make_joint_rows <- function(d, joint_name) {
  if (joint_name == "hip") {
    implant_date <- d$hip_primary_date
    revision_date <- d$hip_revision_date
  } else {
    implant_date <- d$knee_primary_date
    revision_date <- d$knee_revision_date
  }
  
  tibble(
    patient_side_id = d$patient_side_id,
    CODPAT = d$CODPAT,
    joint = joint_name,
    implant_date = as_date_safe(implant_date),
    revision_date = as_date_safe(revision_date),
    death_date = d$death_date,
    censor_date = d$censor_date,
    p1_date = d$p1_date,
    p1_joint = d$p1_joint,
    p2_date = d$p2_date,
    p2_joint = d$p2_joint,
    age_p1 = d$age_p1,
    sex_std = d$sex_std,
    bmi_std = d$bmi_std
  ) %>%
    mutate(
      implant_order = case_when(
        joint == p1_joint & !is.na(implant_date) ~ "first",
        joint == p2_joint & !is.na(p2_date) & !is.na(implant_date) ~ "second",
        TRUE ~ NA_character_
      ),
      both_primaries_observed = !is.na(p2_date) & p2_joint %in% c("hip", "knee"),
      revision_valid = if_else(!is.na(revision_date) & revision_date > implant_date & revision_date <= censor_date, revision_date, as.Date(NA)),
      death_valid = if_else(!is.na(death_date) & death_date > implant_date & death_date <= censor_date, death_date, as.Date(NA)),
      event_date = pmin_date(revision_valid, death_valid, censor_date),
      status = case_when(
        !is.na(revision_valid) & revision_valid == event_date & (is.na(death_valid) | revision_valid <= death_valid) ~ 1L,
        !is.na(death_valid) & death_valid == event_date & (is.na(revision_valid) | death_valid < revision_valid) ~ 2L,
        TRUE ~ 0L
      ),
      time_years = years_between(event_date, implant_date),
      age_at_implant = age_p1 + years_between(implant_date, p1_date),
      implant_year = as.integer(format(implant_date, "%Y"))
    ) %>%
    filter(implant_order %in% c("first", "second"), !is.na(time_years), time_years > 0)
}

implant_raw <- bind_rows(make_joint_rows(base, "hip"), make_joint_rows(base, "knee")) %>%
  mutate(implant_order = factor(implant_order, levels = c("second", "first")))

implant_data <- bind_rows(
  implant_raw %>% mutate(analysis_population = "all_observed_implants"),
  implant_raw %>% filter(both_primaries_observed) %>% mutate(analysis_population = "both_primaries_observed")
) %>%
  mutate(
    analysis_population = factor(analysis_population, levels = c("both_primaries_observed", "all_observed_implants")),
    joint = factor(joint, levels = c("hip", "knee"))
  )

counts <- implant_data %>%
  group_by(analysis_population, joint, implant_order) %>%
  summarise(n_implants = n(), n_patients = n_distinct(CODPAT), revisions = sum(status == 1L), deaths = sum(status == 2L), median_followup = median(time_years), max_followup = max(time_years), .groups = "drop")

number_at_risk <- implant_data %>%
  crossing(time_years_report = report_times) %>%
  group_by(analysis_population, joint, implant_order, time_years_report) %>%
  summarise(n_at_risk = sum(time_years >= time_years_report), .groups = "drop")

cif_report <- implant_data %>%
  group_by(analysis_population, joint, implant_order) %>%
  group_modify(~ aj_at_times(.x$time_years, .x$status, report_times)) %>%
  ungroup()

cif_curve <- implant_data %>%
  group_by(analysis_population, joint, implant_order) %>%
  group_modify(~ aj_at_times(.x$time_years, .x$status, curve_times)) %>%
  ungroup()

revision_difference <- cif_report %>%
  select(analysis_population, joint, implant_order, time_years, revision) %>%
  pivot_wider(names_from = implant_order, values_from = revision) %>%
  mutate(diff_first_minus_second = first - second, diff_percentage_points = 100 * diff_first_minus_second)

fit_grid <- implant_data %>% distinct(analysis_population, joint)
fit_objects <- pmap(fit_grid, function(analysis_population, joint) {
  pop_i <- as.character(analysis_population)
  joint_i <- as.character(joint)
  d <- implant_data %>%
    filter(
      as.character(.data$analysis_population) == .env$pop_i,
      as.character(.data$joint) == .env$joint_i
    )
  fit_revision_cox(d)
})
names(fit_objects) <- paste(fit_grid$analysis_population, fit_grid$joint, sep = "::")

cox_results <- map2_dfr(seq_len(nrow(fit_grid)), fit_objects, function(i, obj) extract_cox(obj, as.character(fit_grid$analysis_population[i]), as.character(fit_grid$joint[i])))
ph_results <- map2_dfr(seq_len(nrow(fit_grid)), fit_objects, function(i, obj) extract_ph(obj, as.character(fit_grid$analysis_population[i]), as.character(fit_grid$joint[i])))

cox_fit_validation <- counts %>%
  group_by(analysis_population, joint) %>%
  summarise(expected_n = sum(n_implants), expected_events = sum(revisions), .groups = "drop") %>%
  left_join(
    cox_results %>%
      group_by(analysis_population, joint) %>%
      summarise(fit_n = first(n), fit_events = first(events), .groups = "drop"),
    by = c("analysis_population", "joint")
  ) %>%
  mutate(
    n_match = expected_n == fit_n,
    events_match = expected_events == fit_events,
    status = if_else(n_match & events_match, "PASS", "CHECK")
  )

if (any(cox_fit_validation$status != "PASS")) {
  print(cox_fit_validation)
  stop("joint-specific Cox risk-set validation failed")
}

qc <- tibble(
  check = c("source cohort is cohort_p1_clean", "no same-day order included", "all follow-up positive", "revision CIF valid", "death CIF valid"),
  value = c(basename(src$source), sum(base$p1_joint == "same_day", na.rm = TRUE), min(implant_data$time_years), range(cif_report$revision, na.rm = TRUE) |> paste(collapse = " to "), range(cif_report$death, na.rm = TRUE) |> paste(collapse = " to ")),
  status = c("PASS", if_else(sum(base$p1_joint == "same_day", na.rm = TRUE) == 0, "PASS", "CHECK"), if_else(min(implant_data$time_years) > 0, "PASS", "CHECK"), if_else(all(cif_report$revision >= -tol & cif_report$revision <= 1 + tol), "PASS", "CHECK"), if_else(all(cif_report$death >= -tol & cif_report$death <= 1 + tol), "PASS", "CHECK"))
)

notes <- tibble(note = c(
  "Secondary descriptive analysis, not a causal analysis.",
  "Time zero is implantation of the specific hip or knee prosthesis.",
  "Revision of the same joint is the event of interest and death is a competing event for cumulative incidence.",
  "The all-observed-implants comparison is structurally asymmetric because a second implant can only be observed after P2.",
  "The both-primaries-observed comparison conditions on reaching P2 and is therefore also subject to post-baseline selection.",
  "Cox models are supporting cause-specific models with robust sandwich inference clustered by patient.",
  "No bootstrap is run here because this analysis is secondary; nonparametric CIFs are descriptive point estimates."
))

plot_dat <- cif_curve %>%
  mutate(
    order_label = recode(as.character(implant_order), first = "Joint implanted first", second = "Joint implanted second"),
    joint_label = recode(as.character(joint), hip = "Hip implant", knee = "Knee implant"),
    population_label = recode(as.character(analysis_population), both_primaries_observed = "Both primaries observed", all_observed_implants = "All observed implants")
  )

p <- ggplot(plot_dat, aes(time_years, revision, linetype = order_label, group = order_label)) +
  geom_step(linewidth = 0.9) +
  facet_grid(population_label ~ joint_label, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = report_times) +
  labs(x = "Years since implantation of the specific joint", y = "Cumulative incidence of same-joint revision", linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(out_plot, p, width = 9, height = 6.5, dpi = 300)

wb <- createWorkbook()
tabs <- list(notes = notes, qc = qc, cox_fit_validation = cox_fit_validation, counts = counts, number_at_risk = number_at_risk, cif_report = cif_report, revision_difference = revision_difference, cox_results = cox_results, ph_results = ph_results, cif_curve = cif_curve)
for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tabs[[nm]])
  setColWidths(wb, nm, cols = 1:50, widths = "auto")
}
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
save(implant_data, counts, number_at_risk, cif_report, cif_curve, revision_difference, cox_results, ph_results, cox_fit_validation, qc, notes, file = out_rdata)

cat("\n16 complete\n")
print(qc)
cat(out_xlsx, "\n")