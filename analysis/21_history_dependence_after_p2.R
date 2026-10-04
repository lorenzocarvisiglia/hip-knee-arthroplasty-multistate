source("config.R")

# ============================================================
# 20_p1_history_dependence_after_direct_P2.R
#
# Sensitivity analysis for history dependence after direct P2.
#
# Question:
# Does the post-P2 transition hazard depend on the time taken
# to reach P2 from P1, beyond time since entry into P2?
#
# This is NOT a causal adjustment for implant order.
# Time P1->P2 is post-baseline and is used here only to assess
# whether the clock-reset / Markov-renewal representation leaves
# residual dependence on the previous transition history.
#
# Models:
#   1) P2_direct -> Rpost
#      frozen main specification:
#      hip + hip*log(time+1), sex, BMI, age strata, patient cluster
#
#   2) P2_direct -> death
#      frozen main specification:
#      proportional hip effect, sex, BMI, age strata, patient cluster
#
# Sensitivity:
#   A) add time P1->P2 linearly, per 5 years
#   B) replace the linear history term by a 3-df natural spline
#
# Outputs:
#   analisi_P1/p1_history_dependence_after_direct_P2.xlsx
#   analisi_P1/p1_history_dependence_after_direct_P2.RData
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

out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
final_models_file <- file.path(out_dir, "p1_multistate_final_models.RData")
out_xlsx <- file.path(out_dir, "p1_history_dependence_after_direct_P2.xlsx")
out_rdata <- file.path(out_dir, "p1_history_dependence_after_direct_P2.RData")

times_report <- c(1, 3, 5, 10, 15)
tol <- 1e-8

if (!file.exists(analysis_file)) stop("Missing: ", analysis_file, ". Run scripts 03 and 04 first.")
if (!file.exists(final_models_file)) stop("Missing: ", final_models_file, ". Run script 04 first.")

env03 <- new.env(parent = globalenv())
env04 <- new.env(parent = globalenv())
load(analysis_file, envir = env03)
load(final_models_file, envir = env04)

required_03 <- c("p1_ms_long", "p1_ms_wide")
required_04 <- c("final_main_hip_hr", "final_main_size")
missing_03 <- required_03[!vapply(required_03, exists, logical(1), envir = env03)]
missing_04 <- required_04[!vapply(required_04, exists, logical(1), envir = env04)]
if (length(missing_03) > 0) stop("Objects missing from script 03 RData: ", paste(missing_03, collapse = ", "))
if (length(missing_04) > 0) stop("Objects missing from script 04 RData: ", paste(missing_04, collapse = ", "))

p1_ms_long <- as_tibble(env03$p1_ms_long)
p1_ms_wide <- as_tibble(env03$p1_ms_wide)
final_main_hip_hr <- as_tibble(env04$final_main_hip_hr)
final_main_size <- as_tibble(env04$final_main_size)

tt_log <- function(x, t, ...) x * log(t + 1)

safe_cox <- function(formula, data, use_tt = FALSE) {
  tryCatch(
    {
      if (use_tt) {
        coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE,
              model = FALSE, tt = tt_log)
      } else {
        coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE,
              model = TRUE)
      }
    },
    error = function(e) e
  )
}

assert_fit <- function(fit, label) {
  if (inherits(fit, "error")) stop(label, " failed: ", fit$message)
  invisible(TRUE)
}

extract_tv_hip_hr <- function(fit, model_label, transition = "P2_to_Rpost",
                              times = times_report) {
  b <- coef(fit); V <- vcov(fit)
  term_main <- "hip_binary"; term_tv <- "tt(hip_binary)"
  if (!all(c(term_main, term_tv) %in% names(b))) stop("Cannot find hip TV coefficients in ", model_label)
  map_dfr(times, function(tt) {
    l <- log(tt + 1)
    log_hr <- b[term_main] + b[term_tv] * l
    var_log_hr <- V[term_main, term_main] + l^2 * V[term_tv, term_tv] +
      2 * l * V[term_main, term_tv]
    se <- sqrt(var_log_hr)
    tibble(transition = transition, model = model_label, time_years = tt,
           hr = exp(log_hr), conf_low = exp(log_hr - 1.96 * se),
           conf_high = exp(log_hr + 1.96 * se))
  })
}

extract_ph_hip_hr <- function(fit, model_label, transition = "P2_to_D") {
  b <- coef(fit); V <- vcov(fit)
  hip_terms <- names(b)[str_detect(names(b), "^hip") & !str_detect(names(b), "^tt")]
  if (length(hip_terms) != 1) stop("Could not identify unique proportional hip term in ", model_label)
  term <- hip_terms[1]; beta <- b[term]; se <- sqrt(V[term, term])
  tibble(transition = transition, model = model_label, time_years = NA_real_,
         hr = exp(beta), conf_low = exp(beta - 1.96 * se),
         conf_high = exp(beta + 1.96 * se))
}

extract_linear_history <- function(fit, transition, term = "time_p1_p2_5y_c") {
  b <- coef(fit); V <- vcov(fit)
  if (!(term %in% names(b))) stop("History term not found in ", transition)
  beta <- b[term]; se <- sqrt(V[term, term]); z <- beta / se
  tibble(transition = transition,
         history_scale = "per 5 additional years from P1 to P2",
         log_hr = beta, hr = exp(beta),
         conf_low = exp(beta - 1.96 * se),
         conf_high = exp(beta + 1.96 * se),
         p_value = 2 * pnorm(abs(z), lower.tail = FALSE))
}

spline_joint_wald <- function(fit, transition) {
  b <- coef(fit); V <- vcov(fit)
  spline_terms <- grep("ns\\(time_p1_p2_years", names(b), value = TRUE)
  if (length(spline_terms) == 0) spline_terms <- grep("splines::ns\\(time_p1_p2_years", names(b), value = TRUE)
  if (length(spline_terms) == 0) stop("Could not identify spline terms for ", transition)
  bs <- b[spline_terms]
  Vs <- V[spline_terms, spline_terms, drop = FALSE]
  stat <- as.numeric(t(bs) %*% solve(Vs, bs))
  df <- length(bs)
  tibble(transition = transition,
         test = "robust joint Wald test for 3-df spline history term",
         df = df, chi_square = stat,
         p_value = pchisq(stat, df = df, lower.tail = FALSE))
}

coef_table <- function(fit, model_label) {
  b <- coef(fit); V <- vcov(fit); se <- sqrt(diag(V))
  tibble(model = model_label, term = names(b), estimate = as.numeric(b),
         robust_se = se, hr = exp(as.numeric(b)),
         conf_low = exp(as.numeric(b) - 1.96 * se),
         conf_high = exp(as.numeric(b) + 1.96 * se))
}

needed_wide <- c("patient_side_id", "first_event_from_p1_ms", "time_p2_from_p1")
if (!all(needed_wide %in% names(p1_ms_wide))) {
  stop("p1_ms_wide missing: ", paste(setdiff(needed_wide, names(p1_ms_wide)), collapse = ", "))
}

history_lookup <- p1_ms_wide %>%
  transmute(patient_side_id = as.character(patient_side_id),
            reached_direct_p2 = first_event_from_p1_ms == "p2",
            time_p1_p2_years = as.numeric(time_p2_from_p1)) %>%
  filter(reached_direct_p2) %>%
  distinct(patient_side_id, .keep_all = TRUE)

if (any(!is.finite(history_lookup$time_p1_p2_years))) stop("Non-finite P1->P2 times among direct-P2 units.")
if (any(history_lookup$time_p1_p2_years <= 0)) stop("Non-positive P1->P2 times among direct-P2 units.")
history_median <- median(history_lookup$time_p1_p2_years, na.rm = TRUE)

prepare_p2_transition <- function(tr) {
  p1_ms_long %>%
    filter(as.character(trans) == tr) %>%
    mutate(
      patient_side_id = as.character(patient_side_id),
      hip = factor(as.character(hip), levels = c("first_knee", "first_hip")),
      hip_binary = if_else(hip == "first_hip", 1, 0),
      age_cat = cut(age_p1, breaks = c(-Inf, 60, 70, 80, Inf),
                    right = FALSE, labels = c("<60", "60-69", "70-79", "80+")),
      sex_ms = factor(as.character(sex_ms), levels = c("male", "female")),
      bmi_p1_f_ms = factor(as.character(bmi_p1_f_ms),
                           levels = c("normal_underweight", "overweight", "obese", "missing"))
    ) %>%
    left_join(history_lookup %>% select(patient_side_id, time_p1_p2_years),
              by = "patient_side_id") %>%
    mutate(time_p1_p2_5y_c = (time_p1_p2_years - history_median) / 5) %>%
    filter(!is.na(Tstop), Tstop > 0, !is.na(status), !is.na(hip),
           !is.na(age_cat), !is.na(sex_ms), !is.na(bmi_p1_f_ms),
           !is.na(CODPAT), is.finite(time_p1_p2_years), time_p1_p2_years > 0)
}

d_rpost <- prepare_p2_transition("P2_to_Rpost")
d_death <- prepare_p2_transition("P2_to_D")

qc_observed <- tibble(
  transition = c("P2_to_Rpost", "P2_to_D"),
  n = c(nrow(d_rpost), nrow(d_death)),
  events = c(sum(d_rpost$status == 1, na.rm = TRUE),
             sum(d_death$status == 1, na.rm = TRUE))
)

qc_expected <- final_main_size %>%
  filter(transition %in% c("P2_to_Rpost", "P2_to_D"),
         model == "final_main_age_stratified") %>%
  select(transition, expected_n = n, expected_events = events)

qc_riskset <- qc_observed %>%
  left_join(qc_expected, by = "transition") %>%
  mutate(n_match = n == expected_n,
         events_match = events == expected_events,
         status = if_else(n_match & events_match, "PASS", "CHECK"))

if (nrow(qc_riskset) != 2 || any(qc_riskset$status != "PASS")) {
  print(qc_riskset)
  stop("P2 risk-set reconstruction does not match frozen script 04.")
}

history_distribution <- d_rpost %>%
  distinct(patient_side_id, hip, time_p1_p2_years) %>%
  group_by(hip) %>%
  summarise(
    n = n(),
    mean_years = mean(time_p1_p2_years),
    sd_years = sd(time_p1_p2_years),
    p05 = quantile(time_p1_p2_years, 0.05, names = FALSE),
    p25 = quantile(time_p1_p2_years, 0.25, names = FALSE),
    median = median(time_p1_p2_years),
    p75 = quantile(time_p1_p2_years, 0.75, names = FALSE),
    p95 = quantile(time_p1_p2_years, 0.95, names = FALSE),
    .groups = "drop"
  )

fit_rpost_main_refit <- safe_cox(
  Surv(Tstop, status) ~ hip_binary + tt(hip_binary) + sex_ms +
    bmi_p1_f_ms + strata(age_cat) + cluster(CODPAT),
  data = d_rpost, use_tt = TRUE
)

fit_death_main_refit <- safe_cox(
  Surv(Tstop, status) ~ hip + sex_ms + bmi_p1_f_ms +
    strata(age_cat) + cluster(CODPAT),
  data = d_death, use_tt = FALSE
)

assert_fit(fit_rpost_main_refit, "Refit frozen P2->Rpost model")
assert_fit(fit_death_main_refit, "Refit frozen P2->death model")

refit_rpost_hr <- extract_tv_hip_hr(fit_rpost_main_refit, "main_refit")
frozen_rpost_hr <- final_main_hip_hr %>%
  filter(transition == "P2_to_Rpost", time_years %in% times_report) %>%
  transmute(transition, model = "frozen_script04", time_years, hr,
            conf_low = conf.low, conf_high = conf.high)

main_reproduction <- frozen_rpost_hr %>%
  select(time_years, frozen_hr = hr) %>%
  left_join(refit_rpost_hr %>% select(time_years, refit_hr = hr), by = "time_years") %>%
  mutate(abs_difference = abs(frozen_hr - refit_hr),
         status = if_else(abs_difference < tol, "PASS", "CHECK"))

if (any(main_reproduction$status != "PASS")) {
  print(main_reproduction)
  stop("The refitted P2->Rpost model does not reproduce the frozen model.")
}

fit_rpost_history_linear <- safe_cox(
  Surv(Tstop, status) ~ hip_binary + tt(hip_binary) + time_p1_p2_5y_c +
    sex_ms + bmi_p1_f_ms + strata(age_cat) + cluster(CODPAT),
  data = d_rpost, use_tt = TRUE
)

fit_death_history_linear <- safe_cox(
  Surv(Tstop, status) ~ hip + time_p1_p2_5y_c + sex_ms +
    bmi_p1_f_ms + strata(age_cat) + cluster(CODPAT),
  data = d_death, use_tt = FALSE
)

assert_fit(fit_rpost_history_linear, "Linear-history P2->Rpost model")
assert_fit(fit_death_history_linear, "Linear-history P2->death model")

history_linear_effects <- bind_rows(
  extract_linear_history(fit_rpost_history_linear, "P2_to_Rpost"),
  extract_linear_history(fit_death_history_linear, "P2_to_D")
)

history_rpost_hr <- extract_tv_hip_hr(fit_rpost_history_linear, "plus_time_P1_to_P2")
main_vs_history_rpost <- bind_rows(frozen_rpost_hr, history_rpost_hr) %>%
  arrange(time_years, model)

refit_death_hr <- extract_ph_hip_hr(fit_death_main_refit, "main_refit", "P2_to_D")
history_death_hr <- extract_ph_hip_hr(fit_death_history_linear, "plus_time_P1_to_P2", "P2_to_D")
main_vs_history_death <- bind_rows(refit_death_hr, history_death_hr)

fit_rpost_history_spline <- safe_cox(
  Surv(Tstop, status) ~ hip_binary + tt(hip_binary) +
    ns(time_p1_p2_years, df = 3) + sex_ms + bmi_p1_f_ms +
    strata(age_cat) + cluster(CODPAT),
  data = d_rpost, use_tt = TRUE
)

fit_death_history_spline <- safe_cox(
  Surv(Tstop, status) ~ hip + ns(time_p1_p2_years, df = 3) +
    sex_ms + bmi_p1_f_ms + strata(age_cat) + cluster(CODPAT),
  data = d_death, use_tt = FALSE
)

assert_fit(fit_rpost_history_spline, "Spline-history P2->Rpost model")
assert_fit(fit_death_history_spline, "Spline-history P2->death model")

spline_joint_tests <- bind_rows(
  spline_joint_wald(fit_rpost_history_spline, "P2_to_Rpost"),
  spline_joint_wald(fit_death_history_spline, "P2_to_D")
)

spline_rpost_hr <- extract_tv_hip_hr(fit_rpost_history_spline, "plus_spline_time_P1_to_P2")
spline_death_hr <- extract_ph_hip_hr(fit_death_history_spline, "plus_spline_time_P1_to_P2", "P2_to_D")

rpost_hr_comparison <- bind_rows(
  frozen_rpost_hr,
  history_rpost_hr,
  spline_rpost_hr
) %>% arrange(time_years, model)

death_hr_comparison <- bind_rows(
  refit_death_hr,
  history_death_hr,
  spline_death_hr
)

model_fit_summary <- tibble(
  transition = c("P2_to_Rpost", "P2_to_Rpost", "P2_to_Rpost",
                 "P2_to_D", "P2_to_D", "P2_to_D"),
  model = c("main_refit", "linear_history", "spline_history",
            "main_refit", "linear_history", "spline_history"),
  n = c(fit_rpost_main_refit$n, fit_rpost_history_linear$n, fit_rpost_history_spline$n,
        fit_death_main_refit$n, fit_death_history_linear$n, fit_death_history_spline$n),
  events = c(fit_rpost_main_refit$nevent, fit_rpost_history_linear$nevent, fit_rpost_history_spline$nevent,
             fit_death_main_refit$nevent, fit_death_history_linear$nevent, fit_death_history_spline$nevent),
  partial_loglik = c(as.numeric(logLik(fit_rpost_main_refit)),
                     as.numeric(logLik(fit_rpost_history_linear)),
                     as.numeric(logLik(fit_rpost_history_spline)),
                     as.numeric(logLik(fit_death_main_refit)),
                     as.numeric(logLik(fit_death_history_linear)),
                     as.numeric(logLik(fit_death_history_spline)))
)

notes <- tibble(note = c(
  "Sensitivity analysis for dependence on previous transition history, not a causal adjustment for implant order.",
  "Time P1->P2 is post-baseline and is evaluated only among units that reached direct P2.",
  "The linear history HR is per 5 additional years between P1 and direct P2.",
  "The 3-df spline checks whether a linear history term is too restrictive.",
  "P2->Rpost retains the frozen time-varying implant-order effect from script 04.",
  "P2->death retains the frozen proportional implant-order specification from script 04.",
  "Robust sandwich inference remains clustered by CODPAT, as in the main analysis.",
  "A material history association would indicate residual dependence on the path to P2 and motivate a more explicit history-dependent or frailty model; it would not by itself identify frailty."
))

coefficient_audit <- bind_rows(
  coef_table(fit_rpost_main_refit, "P2_to_Rpost_main"),
  coef_table(fit_rpost_history_linear, "P2_to_Rpost_linear_history"),
  coef_table(fit_rpost_history_spline, "P2_to_Rpost_spline_history"),
  coef_table(fit_death_main_refit, "P2_to_D_main"),
  coef_table(fit_death_history_linear, "P2_to_D_linear_history"),
  coef_table(fit_death_history_spline, "P2_to_D_spline_history")
)

wb <- createWorkbook()
tabs <- list(
  notes = notes,
  qc_riskset = qc_riskset,
  main_reproduction = main_reproduction,
  history_distribution = history_distribution,
  linear_history_effects = history_linear_effects,
  rpost_hr_comparison = rpost_hr_comparison,
  death_hr_comparison = death_hr_comparison,
  spline_joint_tests = spline_joint_tests,
  model_fit_summary = model_fit_summary,
  coefficient_audit = coefficient_audit
)
for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tabs[[nm]])
  setColWidths(wb, nm, cols = 1:50, widths = "auto")
}
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

save(
  d_rpost, d_death, history_lookup, history_distribution,
  qc_riskset, main_reproduction,
  fit_rpost_main_refit, fit_rpost_history_linear, fit_rpost_history_spline,
  fit_death_main_refit, fit_death_history_linear, fit_death_history_spline,
  history_linear_effects, rpost_hr_comparison, death_hr_comparison,
  spline_joint_tests, model_fit_summary, coefficient_audit, notes,
  file = out_rdata
)

cat("\n============================================================\n")
cat("History-dependence sensitivity after direct P2\n")
cat("============================================================\n\n")
cat("Risk-set QC:\n"); print(qc_riskset)
cat("\nP1->P2 duration by observed implant order:\n"); print(history_distribution)
cat("\nLinear history effects (per +5 years P1->P2):\n"); print(history_linear_effects)
cat("\nFirst-hip vs first-knee HR after P2: main vs history model\n"); print(rpost_hr_comparison)
cat("\nSpline joint Wald tests for prior-history term:\n"); print(spline_joint_tests)
cat("\nWritten:\n"); cat(out_xlsx, "\n"); cat(out_rdata, "\n")