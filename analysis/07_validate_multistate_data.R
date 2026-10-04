source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

if (!file.exists(clean_data_file)) {
  stop("Missing cleaned data file: ", clean_data_file)
}
if (!file.exists(analysis_file)) {
  stop("Missing analysis file: ", analysis_file)
}

clean_env <- new.env(parent = emptyenv())
load(clean_data_file, envir = clean_env)

analysis_env <- new.env(parent = emptyenv())
load(analysis_file, envir = analysis_env)

checks <- list()

if (exists("cohort_p1_clean", envir = clean_env, inherits = FALSE)) {
  cohort <- clean_env$cohort_p1_clean
  
  checks[["nonempty_cohort"]] <- nrow(cohort) > 0
  checks[["unique_patient_side"]] <- if (all(c("CODPAT", "side") %in% names(cohort))) {
    nrow(cohort) == nrow(distinct(cohort, CODPAT, side))
  } else {
    NA
  }
  checks[["valid_order_indicator"]] <- if ("hip" %in% names(cohort)) {
    all(is.na(cohort$hip) | cohort$hip %in% c(0, 1))
  } else {
    NA
  }
  checks[["nonnegative_followup"]] <- if ("follow_up_p1" %in% names(cohort)) {
    all(is.na(cohort$follow_up_p1) | cohort$follow_up_p1 >= 0)
  } else {
    NA
  }
}

qc <- tibble(
  check = names(checks),
  value = unlist(checks),
  status = case_when(
    is.na(value) ~ "NOT_EVALUATED",
    value ~ "PASS",
    TRUE ~ "CHECK"
  )
)

print(qc)

write.xlsx(
  qc,
  file.path(output_dir, "p1_multistate_validation_checks.xlsx"),
  overwrite = TRUE
)

if (any(qc$status == "CHECK")) {
  warning("One or more validation checks require review.")
}
