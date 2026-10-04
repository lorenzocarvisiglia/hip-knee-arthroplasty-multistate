# ============================================================
# 17_p1_ps_diagnostics_for_paper.R
#
# Diagnostica propensity score / overlap weights per paper P1
# Rizzoli/RIPO
#
# Obiettivi:
# 1. Ricostruire il propensity score per first hip vs first knee
#    al baseline P1.
# 2. Calcolare overlap weights:
#       first hip  -> 1 - e(X)
#       first knee -> e(X)
# 3. Produrre diagnostica:
#       - distribuzione PS
#       - distribuzione overlap weights
#       - effective sample size
#       - standardized mean differences prima/dopo weighting
#       - love plot
#       - tabella supplementare pronta
#
# Nota:
# - Questo script NON fa bootstrap.
# - Serve come diagnostica supplementare per difendere
#   l'analisi overlap-weighted.
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c(
    "dplyr",
    "tidyr",
    "stringr",
    "purrr",
    "tibble",
    "openxlsx",
    "ggplot2",
    "scales"
  )
  
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop(
      "Installa questi pacchetti prima di eseguire lo script: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }
  
  invisible(lapply(pkgs, library, character.only = TRUE))
})

# ============================================================
# settings
# ============================================================

source("config.R")

out_dir <- file.path(output_dir, "diagnostics", "propensity_score")
fig_dir <- file.path(out_dir, "figures")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_xlsx <- file.path(out_dir, "p1_ps_diagnostics_for_paper.xlsx")
out_rdata <- file.path(out_dir, "p1_ps_diagnostics_for_paper.RData")
out_txt <- file.path(out_dir, "p1_ps_diagnostics_for_paper_notes.txt")

manual_source_rdata <- clean_data_file
manual_source_object <- "cohort_p1_clean"

min_center_n <- 100
love_plot_top_n <- 40

# ============================================================
# funzioni generali
# ============================================================

clean_names_simple <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\u00a0", " ")
  x <- stringr::str_squish(x)
  x <- stringr::str_to_lower(x)
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x
}

as_tibble_clean <- function(x) {
  out <- tibble::as_tibble(x, .name_repair = "unique")
  names(out) <- make.unique(clean_names_simple(names(out)), sep = "_dup")
  out
}

to_num <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  
  out <- suppressWarnings(as.numeric(x))
  out[is.infinite(out)] <- NA_real_
  out
}

parse_date_vec <- function(x) {
  if (inherits(x, "Date")) {
    return(as.Date(x))
  }
  
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) {
    return(as.Date(x))
  }
  
  if (is.factor(x)) {
    x <- as.character(x)
  }
  
  if (is.numeric(x)) {
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    
    if (is.na(med)) {
      return(as.Date(rep(NA_real_, length(x)), origin = "1970-01-01"))
    }
    
    if (med > 25000) {
      return(as.Date(x, origin = "1899-12-30"))
    } else {
      return(as.Date(x, origin = "1970-01-01"))
    }
  }
  
  x_chr <- as.character(x)
  x_chr <- stringr::str_squish(x_chr)
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_
  
  formats <- c(
    "%Y-%m-%d",
    "%d/%m/%Y",
    "%m/%d/%Y",
    "%d-%m-%Y",
    "%Y/%m/%d",
    "%d.%m.%Y"
  )
  
  out <- rep(as.Date(NA), length(x_chr))
  
  for (fmt in formats) {
    idx <- is.na(out) & !is.na(x_chr)
    
    if (!any(idx)) {
      break
    }
    
    parsed <- suppressWarnings(as.Date(x_chr[idx], format = fmt))
    out[idx] <- parsed
  }
  
  out
}

pmin_date <- function(...) {
  dots <- list(...)
  
  mat <- do.call(
    cbind,
    lapply(dots, function(x) as.numeric(as.Date(x)))
  )
  
  out_num <- apply(
    mat,
    1,
    function(z) {
      if (all(is.na(z))) {
        return(NA_real_)
      }
      min(z, na.rm = TRUE)
    }
  )
  
  as.Date(out_num, origin = "1970-01-01")
}

find_col <- function(dat, exact = character(), regex = character(), required = FALSE, label = "") {
  nm <- names(dat)
  
  exact <- clean_names_simple(exact)
  hit <- intersect(exact, nm)
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  for (pat in regex) {
    idx <- grep(pat, nm, ignore.case = TRUE, perl = TRUE, value = TRUE)
    
    if (length(idx) > 0) {
      return(idx[1])
    }
  }
  
  if (required) {
    stop(
      "Non trovo la colonna richiesta: ", label, "\n",
      "Colonne disponibili:\n",
      paste(nm, collapse = ", ")
    )
  }
  
  NA_character_
}

safe_factor_missing <- function(x, missing_label = "Missing") {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[is.na(x) | x == ""] <- missing_label
  factor(x)
}

lump_small_levels <- function(x, min_n = 50, other_label = "Other_small_centres") {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[is.na(x) | x == ""] <- "Missing"
  
  tab <- table(x, useNA = "no")
  small <- names(tab)[tab < min_n]
  
  x[x %in% small] <- other_label
  factor(x)
}

safe_create_style <- function(...) {
  tryCatch(
    openxlsx::createStyle(...),
    error = function(e) openxlsx::createStyle()
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    formatC(x, format = "f", digits = digits)
  )
}

fmt_pct <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    "",
    paste0(formatC(100 * x, format = "f", digits = digits), "%")
  )
}

standardize_order_group <- function(x) {
  y <- as.character(x)
  y <- stringr::str_to_lower(y)
  y <- stringr::str_replace_all(y, "\\s+", "_")
  y <- stringr::str_replace_all(y, "-", "_")
  y <- stringr::str_squish(y)
  
  dplyr::case_when(
    y %in% c(
      "first_hip", "firsthip", "hip_first", "hip",
      "p1_hip", "hipfirst", "first_hip_group"
    ) ~ "first_hip",
    
    y %in% c(
      "first_knee", "firstknee", "knee_first", "knee",
      "p1_knee", "kneefirst", "first_knee_group"
    ) ~ "first_knee",
    
    stringr::str_detect(y, "hip") & !stringr::str_detect(y, "knee") ~ "first_hip",
    stringr::str_detect(y, "knee") & !stringr::str_detect(y, "hip") ~ "first_knee",
    
    TRUE ~ NA_character_
  )
}

standardize_side <- function(x) {
  y <- as.character(x)
  y <- stringr::str_to_lower(y)
  y <- stringr::str_squish(y)
  
  dplyr::case_when(
    y %in% c("left", "l", "sx", "sinistra", "sinistro", "1") ~ "left",
    y %in% c("right", "r", "dx", "destra", "destro", "2") ~ "right",
    stringr::str_detect(y, "left|sin") ~ "left",
    stringr::str_detect(y, "right|des|dest|dx") ~ "right",
    TRUE ~ y
  )
}

# ============================================================
# trova dataset sorgente
# ============================================================

detect_source_columns <- function(dat) {
  list(
    codpat = find_col(
      dat,
      exact = c("codpat", "patient_id", "id_patient"),
      regex = c("^codpat$", "patient.*id"),
      required = FALSE,
      label = "CODPAT"
    ),
    
    patient_side_id = find_col(
      dat,
      exact = c("patient_side_id", "id_patient_side", "pat_side_id"),
      regex = c("patient.*side.*id", "side.*id"),
      required = FALSE,
      label = "patient_side_id"
    ),
    
    side = find_col(
      dat,
      exact = c("side", "lato", "side_ms"),
      regex = c("^side$", "^lato$"),
      required = FALSE,
      label = "side"
    ),
    
    order_group = find_col(
      dat,
      exact = c("order_group", "implant_order", "p1_order", "first_implant", "first_joint"),
      regex = c("order.*group", "implant.*order", "first.*implant", "first.*joint"),
      required = FALSE,
      label = "order group"
    ),
    
    hip_primary_date = find_col(
      dat,
      exact = c(
        "hip_primary_date",
        "primary_hip_date",
        "date_hip_primary",
        "hip_date_primary",
        "data_hip_primary",
        "dataint_hip",
        "dataint_anca",
        "anca_primary_date"
      ),
      regex = c(
        "hip.*primary.*date",
        "primary.*hip.*date",
        "hip.*prim.*date",
        "anca.*prim.*date",
        "data.*anca.*prim",
        "dataint.*anca",
        "dataint.*hip"
      ),
      required = FALSE,
      label = "hip primary date"
    ),
    
    knee_primary_date = find_col(
      dat,
      exact = c(
        "knee_primary_date",
        "primary_knee_date",
        "date_knee_primary",
        "knee_date_primary",
        "data_knee_primary",
        "dataint_knee",
        "dataint_ginocchio",
        "ginocchio_primary_date"
      ),
      regex = c(
        "knee.*primary.*date",
        "primary.*knee.*date",
        "knee.*prim.*date",
        "ginocchio.*prim.*date",
        "gin.*prim.*date",
        "data.*ginocchio.*prim",
        "dataint.*gin",
        "dataint.*knee"
      ),
      required = FALSE,
      label = "knee primary date"
    ),
    
    p1_date = find_col(
      dat,
      exact = c("p1_date", "date_p1", "first_primary_date"),
      regex = c("^p1.*date", "first.*primary.*date"),
      required = FALSE,
      label = "P1 date"
    ),
    
    age_p1 = find_col(
      dat,
      exact = c("age_p1", "age_at_p1", "p1_age", "age_primary", "age"),
      regex = c("age.*p1", "p1.*age", "^age$"),
      required = FALSE,
      label = "age at P1"
    ),
    
    sex = find_col(
      dat,
      exact = c("sex_ms", "sex", "sesso", "gender"),
      regex = c("^sex$", "sesso", "gender"),
      required = FALSE,
      label = "sex"
    ),
    
    bmi_p1 = find_col(
      dat,
      exact = c("bmi_p1_f_ms", "bmi_p1_f", "bmi_p1", "bmi_cat_p1", "cl_bmi", "bmi_category"),
      regex = c("bmi.*p1", "cl_bmi", "bmi.*cat", "^bmi$"),
      required = FALSE,
      label = "BMI at P1"
    ),
    
    center_p1 = find_col(
      dat,
      exact = c("center_p1", "centre_p1", "hospital_p1", "codosp_p1", "istituto_p1"),
      regex = c("center.*p1", "centre.*p1", "hospital.*p1", "codosp.*p1", "istituto.*p1"),
      required = FALSE,
      label = "center P1"
    )
  )
}

looks_like_ps_source <- function(dat) {
  dat <- as_tibble_clean(dat)
  cols <- detect_source_columns(dat)
  
  has_order <- !is.na(cols$order_group)
  has_dates_for_order <- !is.na(cols$hip_primary_date) && !is.na(cols$knee_primary_date)
  has_id <- !is.na(cols$codpat) || !is.na(cols$patient_side_id)
  has_side <- !is.na(cols$side)
  has_some_cov <- any(!is.na(c(cols$age_p1, cols$sex, cols$bmi_p1, cols$center_p1)))
  
  nrow(dat) > 1000 && has_id && has_side && (has_order || has_dates_for_order) && has_some_cov
}

find_source_data_in_file <- function(path, preferred_objects = c("cohort_p1_clean", "p1_ms_wide", "cohort_p1")) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  
  env <- new.env(parent = emptyenv())
  
  ok <- tryCatch(
    {
      load(path, envir = env)
      TRUE
    },
    error = function(e) FALSE
  )
  
  if (!ok) {
    return(NULL)
  }
  
  objs <- ls(env)
  
  preferred <- intersect(preferred_objects, objs)
  
  for (nm in preferred) {
    if (is.data.frame(env[[nm]])) {
      dat <- as_tibble_clean(env[[nm]])
      
      if (looks_like_ps_source(dat)) {
        return(list(data = dat, path = path, object = nm))
      }
    }
  }
  
  dfs <- objs[vapply(objs, function(nm) is.data.frame(env[[nm]]), logical(1))]
  
  for (nm in dfs) {
    dat <- as_tibble_clean(env[[nm]])
    
    if (looks_like_ps_source(dat)) {
      return(list(data = dat, path = path, object = nm))
    }
  }
  
  NULL
}

find_source_data <- function() {
  candidate_files <- unique(c(
    manual_source_rdata,
    clean_data_file,
    analysis_file
  ))
  candidate_files <- candidate_files[!is.na(candidate_files) & file.exists(candidate_files)]

  for (path in candidate_files) {
    hit <- find_source_data_in_file(
      path,
      preferred_objects = c(manual_source_object, "cohort_p1_clean", "p1_ms_wide", "cohort_p1")
    )
    if (!is.null(hit)) {
      return(hit)
    }
  }

  rdata_files <- list.files(
    output_dir,
    pattern = "\\.RData$|\\.rdata$|\\.Rda$|\\.rda$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  for (path in rdata_files) {
    hit <- find_source_data_in_file(path)
    if (!is.null(hit)) {
      return(hit)
    }
  }

  stop(
    "No compatible RData file was found. Run analysis/01_clean_data.R first, ",
    "or set clean_data_file in config.R."
  )
}

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