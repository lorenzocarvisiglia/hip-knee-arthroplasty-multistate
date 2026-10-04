source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

out_dir <- output_dir
analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
ps_file <- file.path(out_dir, "p1_propensity_score_feasibility_updated.RData")
out_xlsx <- file.path(out_dir, "p1_overlap_weighted_P1_probabilities_B500.xlsx")
out_rdata <- file.path(out_dir, "p1_overlap_weighted_P1_probabilities_B500.RData")
checkpoint_file <- file.path(out_dir, "p1_overlap_weighted_P1_probabilities_B500_checkpoint.rds")
B <- 500
master_seed <- 20260830
checkpoint_every <- 25
resume_from_checkpoint <- FALSE 
times_report <- c(1, 3, 5, 10, 15)
eps <- 1e-6
tol <- 1e-10
causes <- c("P1_to_P2", "P1_to_Rpre", "P1_to_D")

if (!file.exists(analysis_file)) stop("missing P1 analysis objects")
if (!file.exists(ps_file)) stop("run updated script 12 first")

env <- new.env(); load(analysis_file, envir = env)
if (!exists("p1_ms_long", envir = env)) stop("p1_ms_long not found")
long <- as_tibble(env$p1_ms_long)
psenv <- new.env(); load(ps_file, envir = psenv)
for (nm in c("ps_dat", "ps_formula")) if (!exists(nm, envir = psenv)) stop(paste(nm, "not found in updated script 12 output"))
ps_dat <- as_tibble(psenv$ps_dat)
ps_formula <- psenv$ps_formula

if (!all(c("Tstart", "Tstop", "status", "trans", "patient_side_id", "CODPAT") %in% names(long))) stop("required long-data columns missing")

p1_long <- long %>%
  filter(as.character(trans) %in% causes) %>%
  mutate(time_from_P1 = Tstop - Tstart, trans_chr = as.character(trans))

time_check <- p1_long %>% group_by(patient_side_id) %>% summarise(n_times = n_distinct(round(time_from_P1, 12)), n_events = sum(status == 1, na.rm = TRUE), .groups = "drop")
if (any(time_check$n_times != 1)) stop("P1 outgoing rows do not share the same follow-up time")
if (any(time_check$n_events > 1)) stop("more than one P1 first event found for a patient-side")

outcomes <- p1_long %>%
  arrange(patient_side_id, trans_chr) %>%
  group_by(patient_side_id) %>%
  summarise(
    CODPAT_long = first(CODPAT),
    followup = first(time_from_P1),
    event = ifelse(any(status == 1), trans_chr[which(status == 1)[1]], "censor"),
    .groups = "drop"
  )

analysis_dat <- ps_dat %>%
  left_join(outcomes, by = "patient_side_id") %>%
  mutate(CODPAT = as.character(CODPAT), CODPAT_long = as.character(CODPAT_long))
if (any(is.na(analysis_dat$followup))) stop("missing P1 outcome after PS merge")
if (any(analysis_dat$CODPAT != analysis_dat$CODPAT_long)) stop("CODPAT mismatch after merge")
if (any(!analysis_dat$event %in% c("censor", causes))) stop("unexpected P1 event label")

weighted_aj <- function(dat, weight_var, method, replicate_id = NA_integer_) {
  out <- list(); kk <- 1
  for (g in c("first_knee", "first_hip")) {
    d <- dat %>% filter(hip_group == g, is.finite(.data[[weight_var]]), .data[[weight_var]] > 0, is.finite(followup), followup >= 0)
    if (nrow(d) == 0) next
    ev_times <- sort(unique(d$followup[d$event != "censor" & d$followup > 0]))
    S <- 1; F <- setNames(rep(0, length(causes)), causes)
    curve <- tibble(time = 0, survival = 1, P1_to_P2 = 0, P1_to_Rpre = 0, P1_to_D = 0)
    if (length(ev_times) > 0) {
      for (tt in ev_times) {
        risk_w <- sum(d[[weight_var]][d$followup >= tt], na.rm = TRUE)
        if (!is.finite(risk_w) || risk_w <= 0) next
        dA <- setNames(rep(0, length(causes)), causes)
        for (cc in causes) dA[cc] <- sum(d[[weight_var]][d$followup == tt & d$event == cc], na.rm = TRUE) / risk_w
        dA_total <- sum(dA)
        if (!is.finite(dA_total) || dA_total < -tol || dA_total > 1 + tol) stop("invalid weighted AJ hazard increment")
        dA_total <- min(max(dA_total, 0), 1)
        S_before <- S
        F <- F + S_before * dA
        S <- S_before * (1 - dA_total)
        if (S < -tol || any(F < -tol) || S + sum(F) > 1 + 1e-8) stop("invalid weighted AJ probability")
        curve <- bind_rows(curve, tibble(time = tt, survival = S, P1_to_P2 = F["P1_to_P2"], P1_to_Rpre = F["P1_to_Rpre"], P1_to_D = F["P1_to_D"]))
      }
    }
    for (hh in times_report) {
      row <- curve %>% filter(time <= hh) %>% slice_tail(n = 1)
      for (cc in causes) {
        out[[kk]] <- tibble(replicate = replicate_id, method = method, group = g, time_years = hh, transition = cc, probability = as.numeric(row[[cc]]), survival = row$survival, n_units = nrow(d), weighted_n = sum(d[[weight_var]]), effective_sample_size = sum(d[[weight_var]])^2 / sum(d[[weight_var]]^2)); kk <- kk + 1
      }
    }
  }
  bind_rows(out)
}

point_unweighted <- weighted_aj(analysis_dat, "w_unweighted", "unweighted_AJ")
point_overlap <- weighted_aj(analysis_dat, "w_overlap", "overlap_weighted_AJ")
point_probabilities <- bind_rows(point_unweighted, point_overlap)
point_contrasts <- point_probabilities %>% select(method, time_years, transition, group, probability) %>% pivot_wider(names_from = group, values_from = probability) %>% mutate(risk_difference = first_hip - first_knee, risk_difference_pp = 100 * risk_difference)

set.seed(master_seed)
bootstrap_seeds <- sample.int(.Machine$integer.max, B)
probability_results <- vector("list", B)
error_results <- vector("list", B)
warning_results <- vector("list", B)
done <- rep(FALSE, B)

if (resume_from_checkpoint && file.exists(checkpoint_file)) {
  ck <- readRDS(checkpoint_file)
  if (!identical(ck$B, B) || !identical(ck$master_seed, master_seed) || !identical(ck$bootstrap_seeds, bootstrap_seeds)) stop("incompatible checkpoint")
  probability_results <- ck$probability_results; error_results <- ck$error_results; warning_results <- ck$warning_results; done <- ck$done
  message("checkpoint loaded: ", sum(done), "/", B)
  if (all(done)) message("all bootstrap replicates are already complete: running post-processing only")
}

patient_ids <- unique(as.character(analysis_dat$CODPAT))
pb <- txtProgressBar(min = 0, max = B, initial = sum(done), style = 3)

for (b in seq_len(B)) {
  if (done[b]) next
  set.seed(bootstrap_seeds[b])
  draws <- sample(patient_ids, length(patient_ids), replace = TRUE)
  pieces <- lapply(seq_along(draws), function(j) {
    analysis_dat %>% filter(CODPAT == draws[j]) %>% mutate(boot_cluster = paste0("draw_", j), boot_unit_id = paste0("draw_", j, "::", patient_side_id))
  })
  bd <- bind_rows(pieces)
  warns <- character()
  psfit <- tryCatch(withCallingHandlers(
    glm(ps_formula, data = bd, family = binomial()),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  ), error = function(e) e)
  if (inherits(psfit, "error")) {
    error_results[[b]] <- tibble(replicate = b, message = conditionMessage(psfit)); done[b] <- TRUE
  } else {
    ps <- pmin(pmax(as.numeric(predict(psfit, type = "response")), eps), 1 - eps)
    bd$w_unweighted <- 1
    bd$w_overlap <- ifelse(bd$hip_binary == 1, 1 - ps, ps)
    res <- tryCatch(bind_rows(weighted_aj(bd, "w_unweighted", "unweighted_AJ", b), weighted_aj(bd, "w_overlap", "overlap_weighted_AJ", b)), error = function(e) e)
    if (inherits(res, "error")) error_results[[b]] <- tibble(replicate = b, message = conditionMessage(res)) else probability_results[[b]] <- res
    if (length(warns)) warning_results[[b]] <- tibble(replicate = b, message = paste(unique(warns), collapse = " | "))
    done[b] <- TRUE
  }
  setTxtProgressBar(pb, sum(done))
  if (b %% checkpoint_every == 0 || all(done)) {
    saveRDS(list(B = B, master_seed = master_seed, bootstrap_seeds = bootstrap_seeds, probability_results = probability_results, error_results = error_results, warning_results = warning_results, done = done), checkpoint_file)
    message(sprintf("checkpoint saved: %d/%d replicates (%.0f%%)", sum(done), B, 100 * sum(done) / B))
  }
}
close(pb)

bootstrap_probabilities <- bind_rows(probability_results)
bootstrap_errors <- bind_rows(error_results)
bootstrap_warnings <- bind_rows(warning_results)
if (nrow(bootstrap_errors) == 0) bootstrap_errors <- tibble(replicate = integer(), message = character())
if (nrow(bootstrap_warnings) == 0) bootstrap_warnings <- tibble(replicate = integer(), message = character())
if (nrow(bootstrap_probabilities) == 0) stop("no successful bootstrap results")

bootstrap_probability_ci <- bootstrap_probabilities %>%
  group_by(method, group, time_years, transition) %>%
  summarise(n_boot = n(), boot_mean = mean(probability), ci_low = quantile(probability, .025), ci_high = quantile(probability, .975), .groups = "drop") %>%
  left_join(point_probabilities %>% select(method, group, time_years, transition, point_estimate = probability), by = c("method", "group", "time_years", "transition"))

boot_contrasts <- bootstrap_probabilities %>% select(replicate, method, time_years, transition, group, probability) %>% pivot_wider(names_from = group, values_from = probability) %>% mutate(risk_difference = first_hip - first_knee)
bootstrap_contrast_ci <- boot_contrasts %>%
  group_by(method, time_years, transition) %>%
  summarise(n_boot = n(), boot_mean = mean(risk_difference), ci_low = quantile(risk_difference, .025), ci_high = quantile(risk_difference, .975), .groups = "drop") %>%
  left_join(point_contrasts %>% select(method, time_years, transition, point_estimate = risk_difference), by = c("method", "time_years", "transition")) %>%
  mutate(point_pp = 100 * point_estimate, ci_low_pp = 100 * ci_low, ci_high_pp = 100 * ci_high)

success_by_method <- bootstrap_probabilities %>% distinct(replicate, method) %>% count(method, name = "successful_replicates")
qc <- tibble(
  check = c("P1 outcome rows equal PS rows", "one P1 event maximum per unit", "point probabilities valid", "all 500 bootstrap iterations completed", "at least 450 successful overlap replicates"),
  value = c(nrow(analysis_dat) == nrow(ps_dat), all(time_check$n_events <= 1), all(point_probabilities$probability >= -tol & point_probabilities$probability <= 1 + tol), all(done), sum(unique(bootstrap_probabilities$replicate[bootstrap_probabilities$method == "overlap_weighted_AJ"]) %in% seq_len(B)) >= 450),
  status = if_else(value, "PASS", "CHECK")
)

scope_note <- tibble(note = c(
  "Overlap-weighted probabilities are restricted to the three first transitions from P1.",
  "The propensity score uses baseline P1 covariates. The same baseline weights are not used to claim causal effects conditional on reaching P2 or Rpre because those are post-exposure selected risk sets.",
  "Within each implant-order group, cumulative incidence is estimated by a weighted Aalen-Johansen first-event calculation with P2, Rpre and death as competing P1 outcomes.",
  "Patient-level bootstrap resamples patients and retains both sides; the propensity-score model and overlap weights are re-estimated in every replicate.",
  "These estimates are a measured-confounding sensitivity analysis and remain observational."
))

wb <- createWorkbook(); tabs <- list(scope_note = scope_note, qc = qc, point_probabilities = point_probabilities, point_contrasts = point_contrasts, bootstrap_probability_ci = bootstrap_probability_ci, bootstrap_contrast_ci = bootstrap_contrast_ci, success_by_method = success_by_method, bootstrap_errors = bootstrap_errors, bootstrap_warnings = bootstrap_warnings)
for (nm in names(tabs)) { addWorksheet(wb, substr(nm, 1, 31)); writeData(wb, substr(nm, 1, 31), tabs[[nm]]); setColWidths(wb, substr(nm, 1, 31), cols = 1:50, widths = "auto") }
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
save(point_probabilities, point_contrasts, bootstrap_probabilities, bootstrap_probability_ci, bootstrap_contrast_ci, bootstrap_errors, bootstrap_warnings, success_by_method, qc, scope_note, file = out_rdata)
cat("\n14 complete\n"); print(qc); print(success_by_method); cat(out_xlsx, "\n")