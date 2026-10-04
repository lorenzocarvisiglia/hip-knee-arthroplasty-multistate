# ============================================================
# costruzione dataset PS
# ============================================================

source_search <- find_source_data()

source_data <- source_search$data
source_path <- source_search$path
source_object <- source_search$object

cols <- detect_source_columns(source_data)

column_mapping <- tibble::tibble(
  item = names(cols),
  column = unlist(cols, use.names = FALSE)
)

derive_order_from_dates <- function(hip_date, knee_date) {
  dplyr::case_when(
    !is.na(hip_date) & is.na(knee_date) ~ "first_hip",
    is.na(hip_date) & !is.na(knee_date) ~ "first_knee",
    !is.na(hip_date) & !is.na(knee_date) & hip_date < knee_date ~ "first_hip",
    !is.na(hip_date) & !is.na(knee_date) & knee_date < hip_date ~ "first_knee",
    TRUE ~ NA_character_
  )
}

source_n <- nrow(source_data)

hip_primary_date <- if (!is.na(cols$hip_primary_date)) {
  parse_date_vec(source_data[[cols$hip_primary_date]])
} else {
  rep(as.Date(NA), source_n)
}

knee_primary_date <- if (!is.na(cols$knee_primary_date)) {
  parse_date_vec(source_data[[cols$knee_primary_date]])
} else {
  rep(as.Date(NA), source_n)
}

p1_date <- if (!is.na(cols$p1_date)) {
  parse_date_vec(source_data[[cols$p1_date]])
} else {
  pmin_date(hip_primary_date, knee_primary_date)
}

order_group_from_col <- if (!is.na(cols$order_group)) {
  standardize_order_group(source_data[[cols$order_group]])
} else {
  rep(NA_character_, source_n)
}

order_group_from_dates <- derive_order_from_dates(
  hip_date = hip_primary_date,
  knee_date = knee_primary_date
)

order_group <- ifelse(
  !is.na(order_group_from_col),
  order_group_from_col,
  order_group_from_dates
)

CODPAT <- if (!is.na(cols$codpat)) {
  as.character(source_data[[cols$codpat]])
} else {
  as.character(seq_len(source_n))
}

patient_side_id <- if (!is.na(cols$patient_side_id)) {
  as.character(source_data[[cols$patient_side_id]])
} else if (!is.na(cols$side)) {
  paste(CODPAT, as.character(source_data[[cols$side]]), sep = "_")
} else {
  paste(CODPAT, seq_len(source_n), sep = "_")
}

side <- if (!is.na(cols$side)) {
  standardize_side(source_data[[cols$side]])
} else {
  rep("Missing", source_n)
}

age_p1 <- if (!is.na(cols$age_p1)) {
  to_num(source_data[[cols$age_p1]])
} else {
  rep(NA_real_, source_n)
}

sex <- if (!is.na(cols$sex)) {
  source_data[[cols$sex]]
} else {
  rep("Missing", source_n)
}

bmi_p1 <- if (!is.na(cols$bmi_p1)) {
  source_data[[cols$bmi_p1]]
} else {
  rep("Missing", source_n)
}

center_p1 <- if (!is.na(cols$center_p1)) {
  source_data[[cols$center_p1]]
} else {
  rep("Missing", source_n)
}

ps_data <- tibble::tibble(
  CODPAT = CODPAT,
  patient_side_id = patient_side_id,
  side = side,
  order_group = order_group,
  hip_first = as.integer(order_group == "first_hip"),
  p1_date = p1_date,
  p1_year = as.integer(format(p1_date, "%Y")),
  hip_primary_date = hip_primary_date,
  knee_primary_date = knee_primary_date,
  age_p1 = age_p1,
  age_p1_10 = age_p1 / 10,
  sex = safe_factor_missing(sex),
  bmi_p1 = safe_factor_missing(bmi_p1),
  side_f = safe_factor_missing(side),
  center_p1_raw = safe_factor_missing(center_p1),
  center_p1_lumped = lump_small_levels(center_p1, min_n = min_center_n)
) %>%
  filter(order_group %in% c("first_hip", "first_knee")) %>%
  mutate(
    order_group = factor(order_group, levels = c("first_knee", "first_hip")),
    order_group_label = dplyr::recode(
      as.character(order_group),
      "first_knee" = "First knee",
      "first_hip" = "First hip"
    ),
    p1_year_centered = p1_year - median(p1_year, na.rm = TRUE)
  )

# ============================================================
# variabili nel modello PS
# ============================================================

candidate_covariates <- tibble::tibble(
  variable = c(
    "age_p1_10",
    "sex",
    "bmi_p1",
    "p1_year_centered",
    "side_f",
    "center_p1_lumped"
  ),
  type = c(
    "continuous",
    "categorical",
    "categorical",
    "continuous",
    "categorical",
    "categorical"
  ),
  manuscript_label = c(
    "Age at P1, per 10 years",
    "Sex",
    "BMI category at P1",
    "Calendar year of P1",
    "Side",
    paste0("P1 centre, levels with n < ", min_center_n, " lumped")
  )
)

covariate_available <- function(dat, variable, type) {
  if (!variable %in% names(dat)) {
    return(FALSE)
  }
  
  x <- dat[[variable]]
  
  if (type == "continuous") {
    return(sum(!is.na(x)) > 20 && stats::sd(to_num(x), na.rm = TRUE) > 0)
  }
  
  if (type == "categorical") {
    return(dplyr::n_distinct(x[!is.na(x)]) >= 2)
  }
  
  FALSE
}

covariates_used <- candidate_covariates %>%
  rowwise() %>%
  mutate(
    used_in_ps_model = covariate_available(ps_data, variable, type)
  ) %>%
  ungroup()

ps_terms <- covariates_used %>%
  filter(used_in_ps_model) %>%
  pull(variable)

if (length(ps_terms) == 0) {
  stop("Nessuna covariata disponibile per il modello PS.")
}

ps_formula <- as.formula(
  paste0("hip_first ~ ", paste(ps_terms, collapse = " + "))
)

# ============================================================
# fit propensity score
# ============================================================

glm_warnings <- character()

ps_fit <- withCallingHandlers(
  glm(
    formula = ps_formula,
    data = ps_data,
    family = binomial(),
    control = glm.control(maxit = 100)
  ),
  warning = function(w) {
    glm_warnings <<- c(glm_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

ps_pred <- predict(ps_fit, type = "response")
ps_pred_clipped <- pmin(pmax(ps_pred, 1e-6), 1 - 1e-6)

ps_data <- ps_data %>%
  mutate(
    ps_first_hip = ps_pred,
    ps_first_hip_clipped = ps_pred_clipped,
    overlap_weight = ifelse(
      hip_first == 1L,
      1 - ps_first_hip_clipped,
      ps_first_hip_clipped
    ),
    overlap_weight_normalized = overlap_weight / mean(overlap_weight, na.rm = TRUE)
  )

ps_model_coef <- as.data.frame(summary(ps_fit)$coefficients) %>%
  tibble::rownames_to_column("term") %>%
  as_tibble_clean()

coef_names <- names(ps_model_coef)

estimate_col <- grep("estimate", coef_names, value = TRUE)[1]
se_col <- grep("std", coef_names, value = TRUE)[1]
p_col <- grep("pr", coef_names, value = TRUE)[1]

ps_model_coef <- ps_model_coef %>%
  mutate(
    odds_ratio = exp(.data[[estimate_col]]),
    ci_low = exp(.data[[estimate_col]] - 1.96 * .data[[se_col]]),
    ci_high = exp(.data[[estimate_col]] + 1.96 * .data[[se_col]]),
    or_ci = paste0(
      formatC(odds_ratio, format = "f", digits = 2),
      " [",
      formatC(ci_low, format = "f", digits = 2),
      ", ",
      formatC(ci_high, format = "f", digits = 2),
      "]"
    ),
    p_value = .data[[p_col]]
  )

# ============================================================
# overlap / weight summaries
# ============================================================

quantile_summary <- function(x) {
  tibble::tibble(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    p01 = as.numeric(stats::quantile(x, 0.01, na.rm = TRUE)),
    p05 = as.numeric(stats::quantile(x, 0.05, na.rm = TRUE)),
    p25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)),
    median = stats::median(x, na.rm = TRUE),
    p75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)),
    p95 = as.numeric(stats::quantile(x, 0.95, na.rm = TRUE)),
    p99 = as.numeric(stats::quantile(x, 0.99, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}

ps_summary_by_group <- ps_data %>%
  group_by(order_group, order_group_label) %>%
  summarise(
    quantile_summary(ps_first_hip),
    .groups = "drop"
  )

weight_summary_by_group <- ps_data %>%
  group_by(order_group, order_group_label) %>%
  summarise(
    quantile_summary(overlap_weight),
    .groups = "drop"
  )

effective_sample_size <- function(w) {
  if (sum(!is.na(w)) == 0 || sum(w, na.rm = TRUE) <= 0) {
    return(NA_real_)
  }
  
  sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
}

ess_summary <- bind_rows(
  ps_data %>%
    summarise(
      group = "Overall",
      n_unweighted = n(),
      sum_weights = sum(overlap_weight, na.rm = TRUE),
      ess_overlap = effective_sample_size(overlap_weight)
    ),
  
  ps_data %>%
    group_by(order_group_label) %>%
    summarise(
      group = unique(order_group_label),
      n_unweighted = n(),
      sum_weights = sum(overlap_weight, na.rm = TRUE),
      ess_overlap = effective_sample_size(overlap_weight),
      .groups = "drop"
    ) %>%
    select(group, n_unweighted, sum_weights, ess_overlap)
)

support_summary <- ps_data %>%
  group_by(order_group, order_group_label) %>%
  summarise(
    n = n(),
    ps_min = min(ps_first_hip, na.rm = TRUE),
    ps_max = max(ps_first_hip, na.rm = TRUE),
    ps_p01 = as.numeric(quantile(ps_first_hip, 0.01, na.rm = TRUE)),
    ps_p99 = as.numeric(quantile(ps_first_hip, 0.99, na.rm = TRUE)),
    .groups = "drop"
  )

common_support_min <- max(support_summary$ps_min, na.rm = TRUE)
common_support_max <- min(support_summary$ps_max, na.rm = TRUE)

common_support <- ps_data %>%
  mutate(
    in_common_support = ps_first_hip >= common_support_min &
      ps_first_hip <= common_support_max
  ) %>%
  group_by(order_group, order_group_label) %>%
  summarise(
    n = n(),
    n_in_common_support = sum(in_common_support, na.rm = TRUE),
    pct_in_common_support = n_in_common_support / n,
    weighted_pct_in_common_support = sum(overlap_weight[in_common_support], na.rm = TRUE) /
      sum(overlap_weight, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# standardized mean differences
# ============================================================

weighted_mean <- function(x, w) {
  x <- to_num(x)
  w <- to_num(w)
  
  ok <- !is.na(x) & !is.na(w)
  
  if (!any(ok) || sum(w[ok]) <= 0) {
    return(NA_real_)
  }
  
  sum(w[ok] * x[ok]) / sum(w[ok])
}

weighted_var <- function(x, w) {
  x <- to_num(x)
  w <- to_num(w)
  
  ok <- !is.na(x) & !is.na(w)
  
  if (!any(ok) || sum(w[ok]) <= 0) {
    return(NA_real_)
  }
  
  m <- weighted_mean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
}

smd_continuous <- function(dat, variable, w_col = NULL) {
  g <- dat$hip_first
  x <- dat[[variable]]
  w <- if (is.null(w_col)) rep(1, nrow(dat)) else dat[[w_col]]
  
  m1 <- weighted_mean(x[g == 1], w[g == 1])
  m0 <- weighted_mean(x[g == 0], w[g == 0])
  v1 <- weighted_var(x[g == 1], w[g == 1])
  v0 <- weighted_var(x[g == 0], w[g == 0])
  
  denom <- sqrt((v1 + v0) / 2)
  
  smd <- ifelse(is.na(denom) | denom == 0, NA_real_, (m1 - m0) / denom)
  
  tibble::tibble(
    variable = variable,
    level = "",
    variable_level = variable,
    type = "continuous",
    first_hip_value = m1,
    first_knee_value = m0,
    smd = smd,
    abs_smd = abs(smd)
  )
}

smd_categorical <- function(dat, variable, w_col = NULL) {
  g <- dat$hip_first
  x <- factor(dat[[variable]])
  w <- if (is.null(w_col)) rep(1, nrow(dat)) else dat[[w_col]]
  
  levs <- levels(x)
  
  purrr::map_dfr(
    levs,
    function(lv) {
      z <- as.integer(x == lv)
      
      p1 <- weighted_mean(z[g == 1], w[g == 1])
      p0 <- weighted_mean(z[g == 0], w[g == 0])
      
      denom <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
      
      smd <- ifelse(is.na(denom) | denom == 0, NA_real_, (p1 - p0) / denom)
      
      tibble::tibble(
        variable = variable,
        level = lv,
        variable_level = paste0(variable, ": ", lv),
        type = "categorical",
        first_hip_value = p1,
        first_knee_value = p0,
        smd = smd,
        abs_smd = abs(smd)
      )
    }
  )
}

compute_balance <- function(dat, covariates_used, weighted = FALSE) {
  w_col <- if (weighted) "overlap_weight" else NULL
  sample_label <- if (weighted) "Overlap weighted" else "Unweighted"
  
  purrr::pmap_dfr(
    covariates_used %>% filter(used_in_ps_model),
    function(variable, type, manuscript_label, used_in_ps_model) {
      if (type == "continuous") {
        out <- smd_continuous(dat, variable, w_col = w_col)
      } else {
        out <- smd_categorical(dat, variable, w_col = w_col)
      }
      
      out %>%
        mutate(
          manuscript_label = manuscript_label,
          sample = sample_label
        )
    }
  )
}

balance_unweighted <- compute_balance(ps_data, covariates_used, weighted = FALSE)
balance_weighted <- compute_balance(ps_data, covariates_used, weighted = TRUE)

balance_long <- bind_rows(
  balance_unweighted,
  balance_weighted
) %>%
  select(
    sample,
    variable,
    manuscript_label,
    level,
    variable_level,
    type,
    first_hip_value,
    first_knee_value,
    smd,
    abs_smd
  ) %>%
  arrange(variable, level, sample)

balance_wide <- balance_long %>%
  select(variable, manuscript_label, level, variable_level, type, sample, smd, abs_smd) %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = c(smd, abs_smd),
    names_sep = "_"
  ) %>%
  arrange(desc(abs_smd_Unweighted))

balance_summary_by_covariate <- balance_long %>%
  group_by(sample, variable, manuscript_label) %>%
  summarise(
    max_abs_smd = max(abs_smd, na.rm = TRUE),
    mean_abs_smd = mean(abs_smd, na.rm = TRUE),
    n_levels_or_terms = n(),
    .groups = "drop"
  ) %>%
  arrange(sample, desc(max_abs_smd))

balance_overall_summary <- balance_long %>%
  group_by(sample) %>%
  summarise(
    max_abs_smd = max(abs_smd, na.rm = TRUE),
    mean_abs_smd = mean(abs_smd, na.rm = TRUE),
    n_terms = n(),
    n_terms_abs_smd_gt_0_10 = sum(abs_smd > 0.10, na.rm = TRUE),
    n_terms_abs_smd_gt_0_20 = sum(abs_smd > 0.20, na.rm = TRUE),
    .groups = "drop"
  )

love_plot_data <- balance_long %>%
  group_by(variable_level) %>%
  mutate(max_abs_smd_any = max(abs_smd, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_abs_smd_any)) %>%
  filter(variable_level %in% unique(variable_level)[seq_len(min(love_plot_top_n, n_distinct(variable_level)))]) %>%
  mutate(
    variable_level = factor(variable_level, levels = rev(unique(variable_level[order(max_abs_smd_any)])))
  )
