source("config.R")

# ============================================================
# 10_p1_core_audit_and_alignment_updated.R
# consolidated audit after final files 03-09
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(openxlsx)
})

# ============================================================
# settings
# ============================================================

out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tol <- 1e-10

files <- tibble(
  step = c("03", "04", "05", "06", "08", "09"),
  object_file = c(
    "p1_multistate_analysis_objects.RData",
    "p1_multistate_final_models.RData",
    "p1_multistate_state_probabilities.RData",
    "p1_mstate_validation_check.RData",
    "p1_sensitivity_analyses_updated.RData",
    "p1_patient_bootstrap_probabilities_B500.RData"
  )
) %>%
  mutate(
    path = file.path(out_dir, object_file),
    exists = file.exists(path)
  )

if (!all(files$exists)) {
  stop(
    "mancano file necessari: ",
    paste(files$object_file[!files$exists], collapse = ", ")
  )
}

# ============================================================
# load objects in separate environments
# ============================================================

env03 <- new.env()
env04 <- new.env()
env05 <- new.env()
env06 <- new.env()
env08 <- new.env()
env09 <- new.env()

load(files$path[files$step == "03"], envir = env03)
load(files$path[files$step == "04"], envir = env04)
load(files$path[files$step == "05"], envir = env05)
load(files$path[files$step == "06"], envir = env06)
load(files$path[files$step == "08"], envir = env08)
load(files$path[files$step == "09"], envir = env09)

required_objects <- tibble(
  step = c(
    "03", "05", "05", "05", "05",
    "06", "08", "08", "09", "09", "09"
  ),
  object = c(
    "p1_ms_long",
    "coefficient_validation_summary",
    "origin_riskset_check",
    "probability_sum_check",
    "tv_specification_check",
    "overall_validation_summary",
    "bmi_main_validation",
    "joint_probability_sum_check",
    "bootstrap_validation",
    "origin_success",
    "contrast_ci"
  )
)

env_map <- list(
  `03` = env03,
  `04` = env04,
  `05` = env05,
  `06` = env06,
  `08` = env08,
  `09` = env09
)

required_objects <- required_objects %>%
  rowwise() %>%
  mutate(available = exists(object, envir = env_map[[step]])) %>%
  ungroup()

if (!all(required_objects$available)) {
  stop(
    "mancano oggetti richiesti nei file RData: ",
    paste(required_objects$object[!required_objects$available], collapse = ", ")
  )
}

p1_ms_long <- env03$p1_ms_long

# ============================================================
# 1) frozen transition structure and counts
# ============================================================

expected_transition_counts <- tibble(
  trans = c(
    "P1_to_P2",
    "P1_to_Rpre",
    "P1_to_D",
    "P2_to_Rpost",
    "P2_to_D",
    "Rpre_to_P2",
    "Rpre_to_D",
    "Rpost_to_D"
  ),
  expected_n_risk = c(
    9703, 9703, 9703,
    2667, 2667,
    401, 401,
    180
  ),
  expected_events = c(
    2668, 401, 1529,
    180, 567,
    101, 67,
    62
  )
)

observed_transition_counts <- p1_ms_long %>%
  mutate(trans = as.character(trans)) %>%
  group_by(trans) %>%
  summarise(
    n_risk = n(),
    n_patient_sides = n_distinct(patient_side_id),
    n_patients = n_distinct(CODPAT),
    events = sum(status == 1, na.rm = TRUE),
    .groups = "drop"
  )

transition_count_audit <- expected_transition_counts %>%
  left_join(observed_transition_counts, by = "trans") %>%
  mutate(
    diff_n_risk = n_risk - expected_n_risk,
    diff_events = events - expected_events,
    match = diff_n_risk == 0 & diff_events == 0
  )

# ============================================================
# 2) origin risk-set consistency
# ============================================================

risk_ids <- function(tr) {
  p1_ms_long %>%
    filter(as.character(trans) == tr) %>%
    distinct(patient_side_id) %>%
    pull(patient_side_id)
}

p1_ids_p2 <- risk_ids("P1_to_P2")
p1_ids_rpre <- risk_ids("P1_to_Rpre")
p1_ids_d <- risk_ids("P1_to_D")

p2_ids_rpost <- risk_ids("P2_to_Rpost")
p2_ids_d <- risk_ids("P2_to_D")

rpre_ids_p2 <- risk_ids("Rpre_to_P2")
rpre_ids_d <- risk_ids("Rpre_to_D")

rpost_ids_d <- risk_ids("Rpost_to_D")

p1_to_p2_events <- observed_transition_counts %>%
  filter(trans == "P1_to_P2") %>%
  pull(events)

p1_to_rpre_events <- observed_transition_counts %>%
  filter(trans == "P1_to_Rpre") %>%
  pull(events)

p2_to_rpost_events <- observed_transition_counts %>%
  filter(trans == "P2_to_Rpost") %>%
  pull(events)

p2_to_d_events <- observed_transition_counts %>%
  filter(trans == "P2_to_D") %>%
  pull(events)

origin_structure_audit <- tibble(
  check = c(
    "P1 outgoing transitions use identical patient-side risk set",
    "P2 direct outgoing transitions use identical patient-side risk set",
    "Rpre outgoing transitions use identical patient-side risk set",
    "P1_to_Rpre events equal Rpre risk set",
    "P2_to_Rpost events equal Rpost risk set",
    "P1_to_P2 events exceed direct P2 positive-follow-up risk set by exactly one"
  ),
  value = c(
    as.character(setequal(p1_ids_p2, p1_ids_rpre) && setequal(p1_ids_p2, p1_ids_d)),
    as.character(setequal(p2_ids_rpost, p2_ids_d)),
    as.character(setequal(rpre_ids_p2, rpre_ids_d)),
    paste0(p1_to_rpre_events, " vs ", length(rpre_ids_p2)),
    paste0(p2_to_rpost_events, " vs ", length(rpost_ids_d)),
    paste0(p1_to_p2_events, " vs ", length(p2_ids_rpost))
  ),
  pass = c(
    setequal(p1_ids_p2, p1_ids_rpre) && setequal(p1_ids_p2, p1_ids_d),
    setequal(p2_ids_rpost, p2_ids_d),
    setequal(rpre_ids_p2, rpre_ids_d),
    p1_to_rpre_events == length(rpre_ids_p2),
    p2_to_rpost_events == length(rpost_ids_d),
    p1_to_p2_events - length(p2_ids_rpost) == 1
  ),
  note = c(
    "P1 has competing transitions to P2, Rpre, and death.",
    "Only the direct P1->P2 branch contributes to post-P2 models.",
    "Rpre->P2 terminates in the distinct P2_after_Rpre state.",
    "Every P1->Rpre event enters the Rpre origin risk set.",
    "Every direct post-P2 revision event enters Rpost.",
    "The one-unit difference is expected from the P2 event occurring exactly at administrative censoring with zero positive follow-up after P2."
  )
)

# ============================================================
# 3) file 05 probability audit
# ============================================================

coef_validation_summary <- env05$coefficient_validation_summary
origin_riskset_check <- env05$origin_riskset_check
probability_sum_check <- env05$probability_sum_check
tv_specification_check <- env05$tv_specification_check

file05_audit <- tibble(
  check = c(
    "file 05 coefficients reproduce file 04",
    "file 05 origin risk sets are identical within origin",
    "file 05 probability sums are valid",
    "file 05 hip(t) specification matches final specification"
  ),
  value = c(
    as.character(all(coef_validation_summary$all_coefficients_match)),
    as.character(all(origin_riskset_check$same_subject_set)),
    paste0(sum(probability_sum_check$valid), "/", nrow(probability_sum_check)),
    paste0(sum(tv_specification_check$match), "/", nrow(tv_specification_check))
  ),
  pass = c(
    all(coef_validation_summary$all_coefficients_match),
    all(origin_riskset_check$same_subject_set),
    all(probability_sum_check$valid) && max(probability_sum_check$abs_error, na.rm = TRUE) < tol,
    all(tv_specification_check$match)
  )
)

# ============================================================
# 4) file 06 structural mstate validation
# ============================================================

mstate_summary <- env06$overall_validation_summary

file06_audit <- tibble(
  check = c(
    "all origin-specific msprep calls succeeded",
    "manual and mstate counts match",
    "no row-level mismatches",
    "clock-reset durations match"
  ),
  value = c(
    as.character(mstate_summary$all_msprep_ok[1]),
    as.character(mstate_summary$all_counts_match[1]),
    as.character(mstate_summary$total_row_mismatches[1]),
    as.character(mstate_summary$max_abs_duration_diff[1])
  ),
  pass = c(
    isTRUE(mstate_summary$all_msprep_ok[1]),
    isTRUE(mstate_summary$all_counts_match[1]),
    mstate_summary$total_row_mismatches[1] == 0,
    mstate_summary$max_abs_duration_diff[1] < tol
  )
)

# ============================================================
# 5) file 08 sensitivity audit
# ============================================================

bmi_validation <- env08$bmi_main_validation
joint_probability_sum_check <- env08$joint_probability_sum_check
joint_transition_counts <- env08$joint_transition_counts
joint_hip_effects <- env08$joint_hip_effects
bmi_counts <- env08$bmi_counts

n_joint_hip_revision <- joint_transition_counts %>%
  filter(cause == "hip_revision") %>%
  pull(events)

n_joint_knee_revision <- joint_transition_counts %>%
  filter(cause == "knee_revision") %>%
  pull(events)

n_joint_death <- joint_transition_counts %>%
  filter(cause == "death") %>%
  pull(events)

file08_audit <- tibble(
  check = c(
    "file 08 main refit reproduces file 04",
    "joint-specific revision events sum to P2_to_Rpost events",
    "joint-specific death events equal P2_to_D events",
    "joint-specific probability sums are valid",
    "joint-specific final hip effects have no fit errors"
  ),
  value = c(
    paste0(sum(bmi_validation$match), "/", nrow(bmi_validation)),
    paste0(n_joint_hip_revision, " + ", n_joint_knee_revision, " = ", n_joint_hip_revision + n_joint_knee_revision),
    as.character(n_joint_death),
    paste0(nrow(joint_probability_sum_check), " rows"),
    as.character(sum(!is.na(joint_hip_effects$error)))
  ),
  pass = c(
    all(bmi_validation$match),
    n_joint_hip_revision + n_joint_knee_revision == p2_to_rpost_events,
    n_joint_death == p2_to_d_events,
    max(joint_probability_sum_check$abs_error, na.rm = TRUE) < tol,
    all(is.na(joint_hip_effects$error))
  )
)

# ============================================================
# 6) file 09 patient-level bootstrap audit
# ============================================================

bootstrap_validation <- env09$bootstrap_validation
origin_success <- env09$origin_success
bootstrap_errors <- env09$bootstrap_errors
bootstrap_warnings <- env09$bootstrap_warnings

get_boot_value <- function(name) {
  bootstrap_validation %>%
    filter(check == name) %>%
    pull(value) %>%
    first()
}

file09_audit <- tibble(
  check = c(
    "500 bootstrap replicates requested",
    "500 bootstrap replicates completed",
    "all origins have 500 successful replicates",
    "no probability CI rows are missing",
    "no contrast CI rows are missing",
    "no bootstrap fitting errors"
  ),
  value = c(
    as.character(get_boot_value("requested_replicates")),
    as.character(get_boot_value("completed_replicates")),
    paste(origin_success$origin, origin_success$successful_replicates, sep = "=", collapse = "; "),
    as.character(get_boot_value("rows_probability_ci_with_missing_ci")),
    as.character(get_boot_value("rows_contrast_ci_with_missing_ci")),
    as.character(nrow(bootstrap_errors))
  ),
  pass = c(
    get_boot_value("requested_replicates") == 500,
    get_boot_value("completed_replicates") == 500,
    all(origin_success$successful_replicates == 500),
    get_boot_value("rows_probability_ci_with_missing_ci") == 0,
    get_boot_value("rows_contrast_ci_with_missing_ci") == 0,
    nrow(bootstrap_errors) == 0
  )
)

bootstrap_warning_summary <- tibble(
  n_warnings = nrow(bootstrap_warnings),
  warning_replicates = if (nrow(bootstrap_warnings) > 0) {
    paste(unique(bootstrap_warnings$replicate), collapse = ", ")
  } else {
    "none"
  },
  note = "Warnings are retained for audit and are not automatically treated as failed replicates."
)

# ============================================================
# 7) overall critical QC
# ============================================================

critical_checks <- bind_rows(
  transition_count_audit %>%
    transmute(section = "transition_counts", check = trans, pass = match, value = paste0(n_risk, "/", events)),
  origin_structure_audit %>%
    transmute(section = "origin_structure", check, pass, value),
  file05_audit %>%
    transmute(section = "file05", check, pass, value),
  file06_audit %>%
    transmute(section = "file06", check, pass, value),
  file08_audit %>%
    transmute(section = "file08", check, pass, value),
  file09_audit %>%
    transmute(section = "file09", check, pass, value)
)

overall_status <- tibble(
  n_critical_checks = nrow(critical_checks),
  n_passed = sum(critical_checks$pass),
  n_failed = sum(!critical_checks$pass),
  all_critical_checks_pass = all(critical_checks$pass),
  bootstrap_warnings = nrow(bootstrap_warnings),
  status = if_else(all(critical_checks$pass), "CORE ANALYSIS CLOSED", "CHECK REQUIRED")
)

if (!all(critical_checks$pass)) {
  warning("almeno un controllo critico non è passato. vedere critical_checks.")
}

# ============================================================
# save audit workbook
# ============================================================

out_xlsx <- file.path(out_dir, "p1_core_audit_and_alignment_updated.xlsx")
out_rdata <- file.path(out_dir, "p1_core_audit_and_alignment_updated.RData")

notes <- tibble(
  item = c(
    "scope",
    "P2 direct branch",
    "Rpre branch",
    "probabilities",
    "bootstrap"
  ),
  description = c(
    "This audit consolidates final files 03-09 before later sensitivity and reporting scripts are rerun.",
    "Post-P2 models use only direct P1->P2 subjects with positive post-P2 follow-up. One P1->P2 event occurs exactly at administrative censoring and therefore does not enter the P2 risk set.",
    "Rpre->P2 terminates in P2_after_Rpre and is not merged into the direct P2 downstream branch.",
    "File 05 estimates origin-specific clock-reset cumulative probabilities, not a single full-state global-time product integral.",
    "File 09 uses patient-level resampling so both sides of a patient remain in the same sampled cluster."
  )
)

wb <- createWorkbook()

tabs <- list(
  overall_status = overall_status,
  critical_checks = critical_checks,
  required_files = files,
  required_objects = required_objects,
  transition_counts = transition_count_audit,
  origin_structure = origin_structure_audit,
  file05_audit = file05_audit,
  file06_audit = file06_audit,
  file08_audit = file08_audit,
  file09_audit = file09_audit,
  bootstrap_warnings = bootstrap_warnings,
  bootstrap_warning_summary = bootstrap_warning_summary,
  bmi_counts = bmi_counts,
  origin_success = origin_success,
  notes = notes
)

for (nm in names(tabs)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, tabs[[nm]])
  setColWidths(wb, nm, cols = 1:100, widths = "auto")
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

save(
  overall_status,
  critical_checks,
  files,
  required_objects,
  transition_count_audit,
  origin_structure_audit,
  file05_audit,
  file06_audit,
  file08_audit,
  file09_audit,
  bootstrap_warning_summary,
  notes,
  file = out_rdata
)

cat("\n==============================\n")
cat("CORE AUDIT STATUS\n")
cat("==============================\n")
print(overall_status)

cat("\n==============================\n")
cat("FAILED CRITICAL CHECKS\n")
cat("==============================\n")
print(critical_checks %>% filter(!pass))

cat("\noutput salvati in:\n")
cat(out_xlsx, "\n")
cat(out_rdata, "\n")