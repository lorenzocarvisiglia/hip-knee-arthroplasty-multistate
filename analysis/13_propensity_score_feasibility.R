source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(openxlsx)
})

out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
out_xlsx <- file.path(out_dir, "p1_propensity_score_feasibility_updated.xlsx")
out_rdata <- file.path(out_dir, "p1_propensity_score_feasibility_updated.RData")
plot_ps_overlap <- file.path(out_dir, "p1_propensity_score_overlap_updated.png")
plot_weight_distribution <- file.path(out_dir, "p1_propensity_score_weights_updated.png")
plot_love <- file.path(out_dir, "p1_balance_love_plot_updated.png")
min_center_n <- 100
eps <- 1e-6

if (!file.exists(analysis_file)) stop("Missing analysis file: ", analysis_file)

env <- new.env()
load(analysis_file, envir = env)
if (!exists("p1_ms_long", envir = env)) stop("p1_ms_long not found")
p1_ms_long <- env$p1_ms_long

standardize_order <- function(x) {
  z <- str_to_lower(as.character(x))
  case_when(
    z %in% c("1", "first_hip", "hip", "h", "anca", "first hip") ~ "first_hip",
    z %in% c("0", "first_knee", "knee", "k", "ginocchio", "first knee") ~ "first_knee",
    TRUE ~ NA_character_
  )
}

standardize_sex <- function(x) {
  z <- str_to_lower(as.character(x))
  case_when(
    z %in% c("male", "m", "maschio", "1") ~ "male",
    z %in% c("female", "f", "femmina", "0") ~ "female",
    TRUE ~ NA_character_
  )
}

standardize_bmi <- function(x) {
  z <- str_to_lower(as.character(x))
  case_when(
    is.na(x) ~ "missing",
    str_detect(z, "missing|unknown|^na$|non") ~ "missing",
    str_detect(z, "obes") ~ "obese",
    str_detect(z, "over|sovra") ~ "overweight",
    str_detect(z, "normal|under|normo|sotto") ~ "normal_underweight",
    TRUE ~ "missing"
  )
}

make_bmi_numeric <- function(x) {
  case_when(
    is.na(x) ~ "missing",
    x < 25 ~ "normal_underweight",
    x < 30 ~ "overweight",
    TRUE ~ "obese"
  )
}

parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    d1 <- suppressWarnings(as.Date(x, origin = "1970-01-01"))
    d2 <- suppressWarnings(as.Date(x, origin = "1899-12-30"))
    y1 <- suppressWarnings(as.integer(format(d1, "%Y")))
    y2 <- suppressWarnings(as.integer(format(d2, "%Y")))
    ok1 <- mean(!is.na(y1) & y1 > 1980 & y1 < 2100) > .8
    ok2 <- mean(!is.na(y2) & y2 > 1980 & y2 < 2100) > .8
    if (ok2 && !ok1) return(d2)
    return(d1)
  }
  z <- as.character(x)
  out <- rep(as.Date(NA), length(z))
  for (fmt in c("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d", "%m/%d/%Y")) {
    ii <- is.na(out)
    out[ii] <- suppressWarnings(as.Date(z[ii], format = fmt))
  }
  out
}

extract_year_safe <- function(x) {
  if (is.numeric(x) && all(is.na(x) | (x > 1900 & x < 2100))) return(as.integer(x))
  as.integer(format(parse_date_safe(x), "%Y"))
}

clean_center <- function(x) {
  z <- str_squish(as.character(x))
  z[is.na(z) | z == ""] <- "Missing_center"
  z
}

clean_side <- function(x) {
  z <- str_to_lower(as.character(x))
  case_when(
    z %in% c("right", "r", "dx", "destra", "d", "2") ~ "right",
    z %in% c("left", "l", "sx", "sinistra", "s", "1") ~ "left",
    TRUE ~ "missing_side"
  )
}

weighted_mean2 <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok) || sum(w[ok]) <= 0) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_var2 <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (sum(ok) <= 1 || sum(w[ok]) <= 0) return(NA_real_)
  m <- weighted_mean2(x[ok], w[ok])
  sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
}

smd_one <- function(x, a, w) {
  ok <- !is.na(x) & !is.na(a) & !is.na(w) & is.finite(w)
  x <- x[ok]; a <- a[ok]; w <- w[ok]
  if (n_distinct(a) < 2) return(NA_real_)
  m1 <- weighted_mean2(x[a == 1], w[a == 1])
  m0 <- weighted_mean2(x[a == 0], w[a == 0])
  v1 <- weighted_var2(x[a == 1], w[a == 1])
  v0 <- weighted_var2(x[a == 0], w[a == 0])
  den <- sqrt((v1 + v0) / 2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  (m1 - m0) / den
}

effective_sample_size <- function(w) {
  w <- w[is.finite(w) & w >= 0]
  if (length(w) == 0 || sum(w^2) == 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

make_balance_matrix <- function(dat) {
  out <- tibble(
    age_p1 = as.numeric(dat$age_p1),
    p1_year = as.numeric(dat$p1_year)
  )
  for (v in c("sex_ms_std", "bmi_p1_f_ms_std", "side_std", "center_p1_lumped")) {
    f <- factor(dat[[v]])
    for (lv in levels(f)) out[[paste0(v, "=", lv)]] <- as.numeric(f == lv)
  }
  out
}

balance_table <- function(dat, mat, weight_var, label) {
  map_dfr(names(mat), function(v) {
    s <- smd_one(mat[[v]], dat$hip_binary, dat[[weight_var]])
    tibble(variable = v, weighting = label, smd = s, abs_smd = abs(s))
  })
}

p1_base <- p1_ms_long %>%
  filter(as.character(trans) == "P1_to_P2") %>%
  distinct(patient_side_id, .keep_all = TRUE) %>%
  mutate(
    hip_group = standardize_order(hip),
    hip_binary = case_when(hip_group == "first_hip" ~ 1, hip_group == "first_knee" ~ 0, TRUE ~ NA_real_)
  ) %>%
  filter(hip_binary %in% c(0, 1))

if (!"age_p1" %in% names(p1_base) && "age_p1_10" %in% names(p1_base)) p1_base$age_p1 <- 10 * p1_base$age_p1_10
if (!"age_p1" %in% names(p1_base)) stop("age_p1/age_p1_10 not found")
p1_base$age_p1_10 <- p1_base$age_p1 / 10

if ("sex_ms" %in% names(p1_base)) p1_base$sex_ms_std <- standardize_sex(p1_base$sex_ms) else if ("sex" %in% names(p1_base)) p1_base$sex_ms_std <- standardize_sex(p1_base$sex) else stop("sex not found")
p1_base$sex_ms_std <- factor(p1_base$sex_ms_std, levels = c("male", "female"))

if ("bmi_p1_f_ms" %in% names(p1_base)) {
  p1_base$bmi_p1_f_ms_std <- standardize_bmi(p1_base$bmi_p1_f_ms)
} else if ("bmi_p1_f" %in% names(p1_base)) {
  p1_base$bmi_p1_f_ms_std <- standardize_bmi(p1_base$bmi_p1_f)
} else if ("bmi_p1" %in% names(p1_base)) {
  p1_base$bmi_p1_f_ms_std <- make_bmi_numeric(p1_base$bmi_p1)
} else stop("BMI not found")
p1_base$bmi_p1_f_ms_std <- factor(p1_base$bmi_p1_f_ms_std, levels = c("normal_underweight", "overweight", "obese", "missing"))

if ("p1_year" %in% names(p1_base)) p1_base$p1_year <- as.integer(p1_base$p1_year) else if ("p1_date" %in% names(p1_base)) p1_base$p1_year <- extract_year_safe(p1_base$p1_date) else stop("p1 year/date not found")
median_year <- median(p1_base$p1_year, na.rm = TRUE)
p1_base$p1_year_5 <- (p1_base$p1_year - median_year) / 5

if ("side" %in% names(p1_base)) p1_base$side_std <- clean_side(p1_base$side) else if ("side_raw" %in% names(p1_base)) p1_base$side_std <- clean_side(p1_base$side_raw) else p1_base$side_std <- "missing_side"
p1_base$side_std <- factor(p1_base$side_std, levels = c("left", "right", "missing_side"))

if ("center_p1_label" %in% names(p1_base)) p1_base$center_p1_raw <- clean_center(p1_base$center_p1_label) else if ("center_p1" %in% names(p1_base)) p1_base$center_p1_raw <- clean_center(p1_base$center_p1) else p1_base$center_p1_raw <- "Missing_center"

center_lumping_table <- p1_base %>%
  count(center_p1_raw, name = "n_patient_sides") %>%
  arrange(desc(n_patient_sides)) %>%
  mutate(
    keep_as_separate_center = n_patient_sides >= min_center_n,
    center_p1_lumped = if_else(keep_as_separate_center, center_p1_raw, "Other_small_or_missing_centers")
  )
centers_to_keep <- center_lumping_table %>% filter(keep_as_separate_center) %>% pull(center_p1_raw)
center_levels <- unique(c(sort(centers_to_keep), "Other_small_or_missing_centers"))
p1_base$center_p1_lumped <- if_else(p1_base$center_p1_raw %in% centers_to_keep, p1_base$center_p1_raw, "Other_small_or_missing_centers")
p1_base$center_p1_lumped <- factor(p1_base$center_p1_lumped, levels = center_levels)

ps_dat <- p1_base %>%
  select(patient_side_id, CODPAT, hip_group, hip_binary, age_p1, age_p1_10, sex_ms_std, bmi_p1_f_ms_std, p1_year, p1_year_5, side_std, center_p1_lumped) %>%
  filter(if_all(everything(), ~ !is.na(.x))) %>%
  as.data.frame()

ps_formula <- hip_binary ~ age_p1_10 + sex_ms_std + bmi_p1_f_ms_std + p1_year_5 + side_std + center_p1_lumped
warn_msg <- character()
ps_fit <- tryCatch(
  withCallingHandlers(
    glm(ps_formula, data = ps_dat, family = binomial()),
    warning = function(w) { warn_msg <<- c(warn_msg, conditionMessage(w)); invokeRestart("muffleWarning") }
  ),
  error = function(e) e
)
if (inherits(ps_fit, "error")) stop(conditionMessage(ps_fit))

ps_dat$ps_raw <- as.numeric(predict(ps_fit, type = "response"))
ps_dat$ps <- pmin(pmax(ps_dat$ps_raw, eps), 1 - eps)
p_treated <- mean(ps_dat$hip_binary == 1)
ps_dat <- ps_dat %>%
  mutate(
    w_unweighted = 1,
    w_iptw = if_else(hip_binary == 1, p_treated / ps, (1 - p_treated) / (1 - ps)),
    w_overlap = if_else(hip_binary == 1, 1 - ps, ps),
    w_iptw_norm = w_iptw / mean(w_iptw),
    w_overlap_norm = w_overlap / mean(w_overlap)
  )

balance_matrix <- make_balance_matrix(ps_dat)
balance_long <- bind_rows(
  balance_table(ps_dat, balance_matrix, "w_unweighted", "unweighted"),
  balance_table(ps_dat, balance_matrix, "w_iptw", "iptw_stabilized"),
  balance_table(ps_dat, balance_matrix, "w_overlap", "overlap")
)
balance_summary <- balance_long %>%
  group_by(weighting) %>%
  summarise(
    n_variables = n(),
    max_abs_smd = max(abs_smd, na.rm = TRUE),
    mean_abs_smd = mean(abs_smd, na.rm = TRUE),
    n_abs_smd_gt_0_10 = sum(abs_smd > .10, na.rm = TRUE),
    n_abs_smd_gt_0_20 = sum(abs_smd > .20, na.rm = TRUE),
    .groups = "drop"
  )

weight_summary <- bind_rows(lapply(c("w_iptw", "w_overlap", "w_iptw_norm", "w_overlap_norm"), function(v) {
  ps_dat %>% group_by(hip_group) %>% summarise(
    weight_type = v, n = n(), mean = mean(.data[[v]]), sd = sd(.data[[v]]), min = min(.data[[v]]),
    p1 = quantile(.data[[v]], .01), p5 = quantile(.data[[v]], .05), median = median(.data[[v]]),
    p95 = quantile(.data[[v]], .95), p99 = quantile(.data[[v]], .99), max = max(.data[[v]]),
    ess = effective_sample_size(.data[[v]]), .groups = "drop"
  )
}))

ps_model_info <- tibble(
  quantity = c("n_patient_sides_used", "n_patients_used", "n_first_hip", "n_first_knee", "treated_probability_first_hip", "min_ps", "p1_ps", "p5_ps", "median_ps", "p95_ps", "p99_ps", "max_ps", "warnings"),
  value = c(
    nrow(ps_dat), n_distinct(ps_dat$CODPAT), sum(ps_dat$hip_binary == 1), sum(ps_dat$hip_binary == 0), round(p_treated, 4),
    round(min(ps_dat$ps), 4), round(quantile(ps_dat$ps, .01), 4), round(quantile(ps_dat$ps, .05), 4), round(median(ps_dat$ps), 4),
    round(quantile(ps_dat$ps, .95), 4), round(quantile(ps_dat$ps, .99), 4), round(max(ps_dat$ps), 4),
    ifelse(length(warn_msg) == 0, "none", paste(unique(warn_msg), collapse = " | "))
  )
)

ps_summary_by_group <- ps_dat %>% group_by(hip_group) %>% summarise(
  n = n(), mean_ps = mean(ps), sd_ps = sd(ps), min_ps = min(ps), p5_ps = quantile(ps, .05), median_ps = median(ps), p95_ps = quantile(ps, .95), max_ps = max(ps), .groups = "drop"
)

ps_overlap_flags <- tibble(
  check = c("overall_ps_below_0.05", "overall_ps_above_0.95", "overall_ps_between_0.10_and_0.90"),
  n = c(sum(ps_dat$ps < .05), sum(ps_dat$ps > .95), sum(ps_dat$ps >= .10 & ps_dat$ps <= .90))
) %>% mutate(pct = 100 * n / nrow(ps_dat))

max_smd_overlap <- balance_summary %>% filter(weighting == "overlap") %>% pull(max_abs_smd)
max_smd_iptw <- balance_summary %>% filter(weighting == "iptw_stabilized") %>% pull(max_abs_smd)
max_w_iptw <- max(ps_dat$w_iptw)
p99_w_iptw <- quantile(ps_dat$w_iptw, .99)
n_extreme_ps <- sum(ps_dat$ps < .05 | ps_dat$ps > .95)

feasibility_decision <- tibble(
  domain = c("propensity-score overlap", "IPTW weights", "overlap-weight balance", "IPTW balance", "recommended use"),
  diagnostic_value = c(
    paste0(n_extreme_ps, " patient-sides with PS <0.05 or >0.95"),
    paste0("max IPTW = ", round(max_w_iptw, 2), "; p99 = ", round(p99_w_iptw, 2)),
    paste0("max absolute SMD = ", round(max_smd_overlap, 3)),
    paste0("max absolute SMD = ", round(max_smd_iptw, 3)),
    ifelse(max_smd_overlap <= .10, "Use overlap weighting as a measured-confounding sensitivity for P1-origin outcomes.", "Residual imbalance remains; inspect before outcome weighting.")
  ),
  interpretation = c(
    "Describes empirical comparability of the two observed implant-order groups.",
    "Large IPTW values indicate limited positivity for some observations.",
    "Absolute SMD <=0.10 is used as a balance diagnostic.",
    "IPTW is retained only as a diagnostic comparison.",
    "This analysis does not by itself justify a causal interpretation."
  )
)

ps_coef <- summary(ps_fit)$coefficients %>% as.data.frame() %>% tibble::rownames_to_column("term") %>% as_tibble()

qc <- tibble(
  check = c("unique patient-side baseline rows", "two implant-order groups", "finite propensity scores", "finite positive overlap weights", "max overlap SMD <= 0.10"),
  value = c(
    n_distinct(ps_dat$patient_side_id) == nrow(ps_dat), n_distinct(ps_dat$hip_binary) == 2,
    all(is.finite(ps_dat$ps)), all(is.finite(ps_dat$w_overlap) & ps_dat$w_overlap > 0), max_smd_overlap <= .10
  ),
  status = if_else(value, "PASS", "CHECK")
)

p_ps <- ggplot(ps_dat, aes(x = ps, fill = hip_group, colour = hip_group)) + geom_density(alpha = .25, linewidth = .8) + labs(x = "Propensity score", y = "Density", fill = "Observed order", colour = "Observed order") + theme_minimal(base_size = 12)
ggsave(plot_ps_overlap, p_ps, width = 8, height = 5, dpi = 300)

wdat <- ps_dat %>% select(hip_group, w_iptw, w_overlap) %>% pivot_longer(c(w_iptw, w_overlap), names_to = "weight_type", values_to = "weight")
p_w <- ggplot(wdat, aes(x = weight, fill = hip_group)) + geom_histogram(bins = 60, alpha = .65, position = "identity") + facet_wrap(~weight_type, scales = "free_x") + labs(x = "Weight", y = "Count", fill = "Observed order") + theme_minimal(base_size = 12)
ggsave(plot_weight_distribution, p_w, width = 8, height = 5, dpi = 300)

love_vars <- balance_long %>% group_by(variable) %>% summarise(mx = max(abs_smd, na.rm = TRUE), .groups = "drop") %>% arrange(desc(mx)) %>% slice_head(n = 35) %>% pull(variable)
p_love <- balance_long %>% filter(variable %in% love_vars) %>% mutate(variable = factor(variable, levels = rev(love_vars))) %>% ggplot(aes(x = abs_smd, y = variable, shape = weighting)) + geom_vline(xintercept = .10, linetype = "dashed") + geom_point(size = 2.2, position = position_dodge(width = .4)) + labs(x = "Absolute standardized mean difference", y = NULL, shape = "Weighting") + theme_minimal(base_size = 11)
ggsave(plot_love, p_love, width = 9, height = 7, dpi = 300)

wb <- createWorkbook()
tabs <- list(qc = qc, feasibility_decision = feasibility_decision, ps_model_info = ps_model_info, ps_summary_by_group = ps_summary_by_group, ps_overlap_flags = ps_overlap_flags, weight_summary = weight_summary, balance_summary = balance_summary, balance_long = balance_long, center_lumping_table = center_lumping_table, ps_model_coefficients = ps_coef)
for (nm in names(tabs)) { addWorksheet(wb, substr(nm, 1, 31)); writeData(wb, substr(nm, 1, 31), tabs[[nm]]); setColWidths(wb, substr(nm, 1, 31), cols = 1:50, widths = "auto") }
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

save(ps_dat, ps_fit, ps_formula, centers_to_keep, center_levels, median_year, min_center_n, eps, balance_long, balance_summary, weight_summary, feasibility_decision, qc, file = out_rdata)

cat("\n12 complete\n")
print(qc)
print(feasibility_decision)
cat(out_xlsx, "\n", out_rdata, "\n")