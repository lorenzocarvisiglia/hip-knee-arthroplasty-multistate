source("config.R")

#whole-pathway point probability P1 to direct P2 to Rpost

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(survival)
  library(openxlsx)
  library(ggplot2)
  library(scales)
})

out_dir <- output_dir
fig_dir <- file.path(out_dir, "figures_p1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

prob_file <- file.path(out_dir, "p1_multistate_state_probabilities.RData")
out_xlsx <- file.path(out_dir, "p1_pathway_P1_P2_Rpost_point.xlsx")
out_rdata <- file.path(out_dir, "p1_pathway_P1_P2_Rpost_point.RData")
out_plot <- file.path(fig_dir, "p1_pathway_P1_P2_Rpost_point.png")
out_diff_plot <- file.path(fig_dir, "p1_pathway_P1_P2_Rpost_point_difference.png")

if (!file.exists(prob_file)) stop("p1_multistate_state_probabilities.RData not found; run final script 05 first")
load(prob_file)

needed_objects <- c("transition_map", "prob_model_objects", "prob_model_summary", "standardization_patterns", "baseline_increments", "adjusted_probabilities_report")
missing_objects <- needed_objects[!vapply(needed_objects, exists, logical(1))]
if (length(missing_objects) > 0) stop("missing objects from file 05: ", paste(missing_objects, collapse = ", "))

report_times <- c(1, 3, 5, 10, 15)
curve_times <- sort(unique(c(0, seq(0, 15, by = 0.25), report_times)))
tol <- 1e-10

needed_transitions <- c("P1_to_P2", "P1_to_Rpre", "P1_to_D", "P2_to_Rpost", "P2_to_D")
expected_tv <- c(P1_to_P2 = TRUE, P1_to_Rpre = FALSE, P1_to_D = FALSE, P2_to_Rpost = TRUE, P2_to_D = FALSE)

spec_check <- prob_model_summary %>%
  filter(trans %in% needed_transitions) %>%
  transmute(trans, origin, destination, model_type, hip_time_varying, model_ok, events, expected_tv = unname(expected_tv[trans]), match_tv = hip_time_varying == expected_tv)
if (nrow(spec_check) != length(needed_transitions) || !all(spec_check$model_ok) || !all(spec_check$match_tv)) {
  print(spec_check)
  stop("file 05 model specification does not match the frozen pathway specification")
}

make_pattern <- function(age_stratum, sex_value, bmi_value, hip_value) {
  tibble(
    age_cat = factor(age_stratum, levels = c("<60", "60-69", "70-79", "80+")),
    sex_ms = factor(sex_value, levels = c("male", "female")),
    bmi_p1_f_ms = factor(bmi_value, levels = c("normal_underweight", "overweight", "obese", "missing")),
    hip_binary = as.numeric(hip_value)
  )
}

calc_lp_times <- function(pattern, beta, times) {
  if (!is.data.frame(pattern)) pattern <- as.data.frame(pattern)
  n <- max(nrow(pattern), length(times))
  if (nrow(pattern) == 1 && n > 1) pattern <- pattern[rep(1, n), , drop = FALSE]
  if (length(times) == 1 && n > 1) times <- rep(times, n)
  if (nrow(pattern) != n || length(times) != n) stop("incompatible pattern/time dimensions")
  
  hip_value <- as.numeric(pattern$hip_binary)
  sex_value <- as.character(pattern$sex_ms)
  bmi_value <- as.character(pattern$bmi_p1_f_ms)
  lp <- rep(0, n)
  
  if ("hipfirst_hip" %in% names(beta)) lp <- lp + beta["hipfirst_hip"] * hip_value
  if ("hip_binary" %in% names(beta)) lp <- lp + beta["hip_binary"] * hip_value
  if ("tt(hip_binary)" %in% names(beta)) lp <- lp + beta["tt(hip_binary)"] * hip_value * log(times + 1)
  if ("sex_msfemale" %in% names(beta)) lp <- lp + beta["sex_msfemale"] * as.numeric(sex_value == "female")
  if ("bmi_p1_f_msoverweight" %in% names(beta)) lp <- lp + beta["bmi_p1_f_msoverweight"] * as.numeric(bmi_value == "overweight")
  if ("bmi_p1_f_msobese" %in% names(beta)) lp <- lp + beta["bmi_p1_f_msobese"] * as.numeric(bmi_value == "obese")
  if ("bmi_p1_f_msmissing" %in% names(beta)) lp <- lp + beta["bmi_p1_f_msmissing"] * as.numeric(bmi_value == "missing")
  as.numeric(lp)
}

hazard_table_for_pattern <- function(transitions, age_stratum, sex_value, bmi_value, hip_value, max_time = 15) {
  pat <- make_pattern(age_stratum, sex_value, bmi_value, hip_value)
  map_dfr(transitions, function(tr) {
    inc <- baseline_increments %>%
      filter(.data$trans == .env$tr, .data$age_stratum == .env$age_stratum, is.finite(time), is.finite(dLambda0), time <= max_time)
    if (nrow(inc) == 0) return(tibble(time = numeric(), trans = character(), dH = numeric()))
    beta <- coef(prob_model_objects[[tr]]$fit)
    lp <- calc_lp_times(pat, beta, inc$time)
    tibble(time = inc$time, trans = tr, dH = inc$dLambda0 * exp(lp))
  })
}

build_p2_rpost_step <- function(age_stratum, sex_value, bmi_value, hip_value, max_time = 15) {
  h <- hazard_table_for_pattern(c("P2_to_Rpost", "P2_to_D"), age_stratum, sex_value, bmi_value, hip_value, max_time)
  times <- sort(unique(h$time))
  if (length(times) == 0) return(tibble(time = 0, cif_rpost = 0, survival_p2 = 1))
  
  S <- 1
  F <- 0
  out <- tibble(time = 0, cif_rpost = 0, survival_p2 = 1)
  for (tt in times) {
    z <- h %>% filter(abs(time - tt) < tol)
    d_r <- sum(z$dH[z$trans == "P2_to_Rpost"])
    d_d <- sum(z$dH[z$trans == "P2_to_D"])
    d_total <- d_r + d_d
    if (!is.finite(d_total) || d_total < -tol) stop("invalid P2 cumulative-hazard increment")
    S0 <- S
    if (d_total > tol) {
      event_prob <- 1 - exp(-d_total)
      F <- F + S0 * event_prob * d_r / d_total
      S <- S0 * exp(-d_total)
    }
    out <- bind_rows(out, tibble(time = tt, cif_rpost = F, survival_p2 = S))
  }
  out
}

eval_step <- function(step, x, value_col) {
  x <- pmax(as.numeric(x), 0)
  idx <- findInterval(x, step$time)
  idx[idx < 1] <- 1
  step[[value_col]][idx]
}

build_p1_entry_mass <- function(age_stratum, sex_value, bmi_value, hip_value, max_time = 15) {
  h <- hazard_table_for_pattern(c("P1_to_P2", "P1_to_Rpre", "P1_to_D"), age_stratum, sex_value, bmi_value, hip_value, max_time)
  times <- sort(unique(h$time))
  if (length(times) == 0) return(tibble(time = numeric(), enter_p2_mass = numeric(), survival_p1_after = numeric()))
  
  S <- 1
  out <- vector("list", length(times))
  for (i in seq_along(times)) {
    tt <- times[i]
    z <- h %>% filter(abs(time - tt) < tol)
    d_p2 <- sum(z$dH[z$trans == "P1_to_P2"])
    d_r <- sum(z$dH[z$trans == "P1_to_Rpre"])
    d_d <- sum(z$dH[z$trans == "P1_to_D"])
    d_total <- d_p2 + d_r + d_d
    if (!is.finite(d_total) || d_total < -tol) stop("invalid P1 cumulative-hazard increment")
    S0 <- S
    enter <- 0
    if (d_total > tol) {
      event_prob <- 1 - exp(-d_total)
      enter <- S0 * event_prob * d_p2 / d_total
      S <- S0 * exp(-d_total)
    }
    out[[i]] <- tibble(time = tt, enter_p2_mass = enter, survival_p1_after = S)
  }
  bind_rows(out)
}

predict_pattern_pathway <- function(age_stratum, sex_value, bmi_value, hip_value, horizons = curve_times) {
  p2_step <- build_p2_rpost_step(age_stratum, sex_value, bmi_value, hip_value, max(horizons))
  p1_mass <- build_p1_entry_mass(age_stratum, sex_value, bmi_value, hip_value, max(horizons))
  
  map_dfr(horizons, function(hh) {
    eligible <- p1_mass %>% filter(time <= hh + tol, enter_p2_mass > 0)
    if (nrow(eligible) == 0) return(tibble(time_years = hh, probability = 0))
    residual <- hh - eligible$time
    cond_rpost <- eval_step(p2_step, residual, "cif_rpost")
    tibble(time_years = hh, probability = sum(eligible$enter_p2_mass * cond_rpost))
  })
}

p1_patterns <- standardization_patterns %>%
  filter(as.character(origin) == "P1", weight > 0) %>%
  transmute(age_cat = as.character(age_cat), sex_ms = as.character(sex_ms), bmi_p1_f_ms = as.character(bmi_p1_f_ms), n_pattern, weight, n_standardization)
if (nrow(p1_patterns) == 0 || abs(sum(p1_patterns$weight) - 1) > 1e-8) stop("invalid P1 standardization patterns")

predict_group <- function(hip_value) {
  g <- if_else(hip_value == 1, "first_hip", "first_knee")
  pred <- map_dfr(seq_len(nrow(p1_patterns)), function(i) {
    p <- p1_patterns[i, ]
    predict_pattern_pathway(p$age_cat, p$sex_ms, p$bmi_p1_f_ms, hip_value, curve_times) %>% mutate(weight = p$weight)
  })
  pred %>% group_by(time_years) %>% summarise(probability = sum(probability * weight), .groups = "drop") %>% mutate(group = g, n_standardization = unique(p1_patterns$n_standardization))
}

message("computing whole-pathway point estimates from frozen file 05 models")
pathway_curve <- bind_rows(predict_group(0), predict_group(1)) %>% arrange(group, time_years)
pathway_report <- pathway_curve %>% filter(time_years %in% report_times)

pathway_contrast <- pathway_report %>%
  select(time_years, group, probability) %>%
  pivot_wider(names_from = group, values_from = probability) %>%
  mutate(risk_difference = first_hip - first_knee, risk_difference_pp = 100 * risk_difference)

p1_p2_reference <- adjusted_probabilities_report %>%
  filter(as.character(origin) == "P1", as.character(state) == "P2", time_years %in% report_times) %>%
  transmute(time_years, group = as.character(group), p1_to_p2_probability = probability)

validation <- pathway_report %>%
  left_join(p1_p2_reference, by = c("time_years", "group")) %>%
  mutate(
    within_bounds = is.finite(probability) & probability >= -tol & probability <= 1 + tol,
    below_p1_to_p2 = probability <= p1_to_p2_probability + tol
  )

monotonicity <- pathway_curve %>%
  group_by(group) %>%
  summarise(min_increment = min(diff(probability)), monotone = min_increment >= -1e-10, .groups = "drop")

qc <- tibble(
  check = c("five frozen model specifications match", "P1 standardization weights sum to one", "all pathway probabilities valid", "pathway never exceeds P1 to P2 probability", "pathway curves nondecreasing"),
  value = c(all(spec_check$match_tv), abs(sum(p1_patterns$weight) - 1) < 1e-8, all(validation$within_bounds), all(validation$below_p1_to_p2), all(monotonicity$monotone)),
  status = if_else(value, "PASS", "CHECK")
)
if (any(qc$status != "PASS")) { print(qc); stop("pathway QC failed") }

notes <- tibble(note = c(
  "Point estimate only. No additional bootstrap is run in this script.",
  "The pathway is P1 -> direct P2 -> Rpost and excludes P2 reached after Rpre.",
  "The calculation uses the frozen final file 05 Cox models, Efron age-stratified baseline increments, and final hip time-varying specifications.",
  "The two clocks are respected: P1 hazards use time since P1 and P2 hazards use residual time since entry into P2.",
  "At each event time cumulative-hazard increments are converted to event probabilities using 1-exp(-dH_total), with cause allocation proportional to dH.",
  "Both hypothetical implant-order groups are standardized to the same empirical P1 age-stratum, sex, and BMI distribution.",
  "If this pathway remains a main inferential result, uncertainty can be added later with a dedicated patient-level bootstrap."
))

plot_dat <- pathway_curve %>% mutate(group_label = recode(group, first_knee = "First knee", first_hip = "First hip"))
p <- ggplot(plot_dat, aes(time_years, probability, linetype = group_label, group = group_label)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = report_times) +
  labs(x = "Years from P1", y = "Whole-pathway probability", linetype = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(out_plot, p, width = 7.5, height = 4.8, dpi = 300)

pd <- pathway_contrast
p2 <- ggplot(pd, aes(time_years, risk_difference_pp)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_line(linewidth = 0.8) + geom_point(size = 2) + scale_x_continuous(breaks = report_times) + labs(x = "Years from P1", y = "First hip minus first knee (percentage points)") + theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(out_diff_plot, p2, width = 7, height = 4.5, dpi = 300)

wb <- createWorkbook()
tabs <- list(notes = notes, qc = qc, model_specification = spec_check, standardization_patterns = p1_patterns, pathway_report = pathway_report, pathway_contrast = pathway_contrast, validation = validation, monotonicity = monotonicity, pathway_curve = pathway_curve)
for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tabs[[nm]])
  setColWidths(wb, nm, cols = 1:50, widths = "auto")
}
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
save(pathway_curve, pathway_report, pathway_contrast, validation, monotonicity, qc, notes, spec_check, p1_patterns, file = out_rdata)

cat("\n18 complete\n")
print(qc)
print(pathway_contrast)
cat(out_xlsx, "\n")