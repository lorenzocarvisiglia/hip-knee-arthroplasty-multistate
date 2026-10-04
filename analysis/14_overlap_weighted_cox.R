source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(survival)
  library(openxlsx)
})

out_dir <- output_dir
analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
ps_file <- file.path(out_dir, "p1_propensity_score_feasibility_updated.RData")
out_xlsx <- file.path(out_dir, "p1_weighted_cox_overlap_updated.xlsx")
out_rdata <- file.path(out_dir, "p1_weighted_cox_overlap_updated.RData")
times_report <- c(1, 3, 5, 10, 15)

if (!file.exists(analysis_file)) stop("missing P1 analysis objects")
if (!file.exists(ps_file)) stop("run updated script 12 first")

env <- new.env(); load(analysis_file, envir = env)
if (!exists("p1_ms_long", envir = env)) stop("p1_ms_long not found")
dat <- as_tibble(env$p1_ms_long)
psenv <- new.env(); load(ps_file, envir = psenv)
if (!exists("ps_dat", envir = psenv)) stop("ps_dat not found")
ps_dat <- as_tibble(psenv$ps_dat)

standardize_order <- function(x) {
  z <- str_to_lower(as.character(x))
  case_when(z %in% c("1", "first_hip", "hip", "h", "anca", "first hip") ~ "first_hip", z %in% c("0", "first_knee", "knee", "k", "ginocchio", "first knee") ~ "first_knee", TRUE ~ NA_character_)
}
standardize_sex <- function(x) {
  z <- str_to_lower(as.character(x)); case_when(z %in% c("male", "m", "maschio", "1") ~ "male", z %in% c("female", "f", "femmina", "0") ~ "female", TRUE ~ NA_character_)
}
standardize_bmi <- function(x) {
  z <- str_to_lower(as.character(x)); case_when(is.na(x) ~ "missing", str_detect(z, "missing|unknown|^na$|non") ~ "missing", str_detect(z, "obes") ~ "obese", str_detect(z, "over|sovra") ~ "overweight", str_detect(z, "normal|under|normo|sotto") ~ "normal_underweight", TRUE ~ "missing")
}

if ("Tstop" %in% names(dat) && "Tstart" %in% names(dat)) dat$analysis_time <- dat$Tstop - dat$Tstart else if ("Tstop" %in% names(dat)) dat$analysis_time <- dat$Tstop else stop("time variables not found")
dat$hip_group <- standardize_order(dat$hip)
dat$hip_binary <- ifelse(dat$hip_group == "first_hip", 1, ifelse(dat$hip_group == "first_knee", 0, NA))
if (!"age_p1" %in% names(dat) && "age_p1_10" %in% names(dat)) dat$age_p1 <- 10 * dat$age_p1_10
if (!"age_p1" %in% names(dat)) stop("age not found")
dat$age_stratum <- cut(dat$age_p1, breaks = c(-Inf, 60, 70, 80, Inf), right = FALSE, labels = c("<60", "60-69", "70-79", "80+"))
if ("sex_ms" %in% names(dat)) dat$sex_std <- standardize_sex(dat$sex_ms) else if ("sex" %in% names(dat)) dat$sex_std <- standardize_sex(dat$sex) else stop("sex not found")
dat$sex_std <- factor(dat$sex_std, levels = c("male", "female"))
if ("bmi_p1_f_ms" %in% names(dat)) dat$bmi_std <- standardize_bmi(dat$bmi_p1_f_ms) else if ("bmi_p1_f" %in% names(dat)) dat$bmi_std <- standardize_bmi(dat$bmi_p1_f) else stop("BMI factor not found")
dat$bmi_std <- factor(dat$bmi_std, levels = c("normal_underweight", "overweight", "obese", "missing"))

dat <- dat %>% left_join(ps_dat %>% select(patient_side_id, w_overlap, ps), by = "patient_side_id")

transitions <- c("P1_to_P2", "P1_to_Rpre", "P1_to_D")

safe_fit <- function(tr, model_type) {
  d <- dat %>% filter(as.character(trans) == tr, analysis_time > 0, hip_binary %in% c(0, 1), !is.na(status), !is.na(CODPAT)) %>% droplevels()
  tv <- tr == "P1_to_P2"
  if (model_type == "unweighted_adjusted") {
    if (tv) f <- Surv(analysis_time, status) ~ hip_binary + tt(hip_binary) + sex_std + bmi_std + strata(age_stratum) + cluster(CODPAT)
    else f <- Surv(analysis_time, status) ~ hip_binary + sex_std + bmi_std + strata(age_stratum) + cluster(CODPAT)
    w <- NULL
  } else if (model_type == "overlap_weighted_marginal") {
    d <- d %>% filter(is.finite(w_overlap), w_overlap > 0)
    if (tv) f <- Surv(analysis_time, status) ~ hip_binary + tt(hip_binary) + cluster(CODPAT)
    else f <- Surv(analysis_time, status) ~ hip_binary + cluster(CODPAT)
    w <- d$w_overlap
  } else stop("unknown model type")
  
  warns <- character()
  fit <- tryCatch(withCallingHandlers(
    coxph(f, data = d, weights = w, robust = TRUE, ties = "efron", tt = function(x, t, ...) x * log(t + 1), x = TRUE, y = TRUE),
    warning = function(z) { warns <<- c(warns, conditionMessage(z)); invokeRestart("muffleWarning") }
  ), error = function(e) e)
  list(fit = fit, warnings = unique(warns), n = nrow(d), nevent = sum(d$status == 1), n_patients = n_distinct(d$CODPAT))
}

extract_result <- function(obj, tr, model_type) {
  if (inherits(obj$fit, "error")) return(tibble(transition = tr, model = model_type, time_years = NA_real_, hr = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_, n = obj$n, nevent = obj$nevent, n_patients = obj$n_patients, status = "FAILED", note = conditionMessage(obj$fit)))
  b <- coef(obj$fit); V <- vcov(obj$fit)
  if (tr == "P1_to_P2") {
    nm_tv <- grep("tt\\(hip_binary\\)", names(b), value = TRUE)[1]
    if (is.na(nm_tv) || !"hip_binary" %in% names(b)) stop("time-varying hip coefficients not found")
    map <- lapply(times_report, function(tt) {
      l <- log(tt + 1); est <- b["hip_binary"] + l * b[nm_tv]
      se <- sqrt(V["hip_binary", "hip_binary"] + l^2 * V[nm_tv, nm_tv] + 2 * l * V["hip_binary", nm_tv])
      tibble(time_years = tt, hr = exp(est), lower_95 = exp(est - 1.96 * se), upper_95 = exp(est + 1.96 * se), p_value = 2 * pnorm(-abs(est / se)))
    })
    out <- bind_rows(map)
  } else {
    est <- b["hip_binary"]; se <- sqrt(V["hip_binary", "hip_binary"])
    out <- tibble(time_years = NA_real_, hr = exp(est), lower_95 = exp(est - 1.96 * se), upper_95 = exp(est + 1.96 * se), p_value = 2 * pnorm(-abs(est / se)))
  }
  out %>% mutate(transition = tr, model = model_type, n = obj$n, nevent = obj$nevent, n_patients = obj$n_patients, status = ifelse(length(obj$warnings) == 0, "OK", "WARNING"), note = paste(obj$warnings, collapse = " | ")) %>% select(transition, model, time_years, hr, lower_95, upper_95, p_value, n, nevent, n_patients, status, note)
}

fits <- list(); results <- list(); k <- 1
for (tr in transitions) for (mt in c("unweighted_adjusted", "overlap_weighted_marginal")) {
  cat("fit:", tr, "|", mt, "\n")
  obj <- safe_fit(tr, mt); fits[[paste(tr, mt, sep = "__")]] <- obj$fit; results[[k]] <- extract_result(obj, tr, mt); k <- k + 1
}
results <- bind_rows(results)

comparison <- results %>% select(transition, model, time_years, hr) %>% tidyr::pivot_wider(names_from = model, values_from = hr) %>% mutate(relative_change_overlap_vs_main = overlap_weighted_marginal / unweighted_adjusted - 1)

qc <- tibble(
  check = c("PS weights merged for all P1 rows", "six model blocks returned", "no failed models", "post-P1 transitions excluded from PS-weighted Cox"),
  value = c(sum(is.na(dat$w_overlap) & as.character(dat$trans) %in% transitions) == 0, n_distinct(paste(results$transition, results$model)) == 6, sum(results$status == "FAILED") == 0, TRUE),
  status = if_else(value, "PASS", "CHECK")
)

scope_note <- tibble(note = c(
  "Overlap-weighted Cox models are restricted to P1-origin transitions because the propensity score is defined from baseline P1 covariates.",
  "P2 and Rpre conditional risk sets are post-exposure selected populations; baseline overlap weighting is not used there as a causal adjustment.",
  "Weighted hazard ratios are supporting summaries. The overlap-weighted P1 cumulative probabilities in script 14 are the main PS sensitivity quantities."
))

wb <- createWorkbook(); tabs <- list(scope_note = scope_note, qc = qc, hip_hr = results, comparison = comparison)
for (nm in names(tabs)) { addWorksheet(wb, nm); writeData(wb, nm, tabs[[nm]]); setColWidths(wb, nm, cols = 1:30, widths = "auto") }
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
save(fits, results, comparison, qc, scope_note, file = out_rdata)
cat("\n13 complete\n"); print(qc); cat(out_xlsx, "\n")