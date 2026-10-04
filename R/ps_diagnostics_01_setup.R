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
    if (!is.null(hit)) return(hit)
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
    if (!is.null(hit)) return(hit)
  }

  stop(
    "No compatible RData file was found. Run analysis/01_clean_data.R first, ",
    "or set clean_data_file in config.R."
  )
}