source("config.R")

#final thesis tables from frozen main and overlap-weighted analyses

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
})

out_dir <- output_dir
paper_dir <- file.path(out_dir, "paper_outputs")
dir.create(paper_dir, showWarnings = FALSE, recursive = TRUE)

analysis_file <- file.path(out_dir, "p1_multistate_analysis_objects.RData")
main_boot_file <- file.path(out_dir, "p1_patient_bootstrap_probabilities_B500.RData")
overlap_file <- file.path(out_dir, "p1_overlap_weighted_P1_probabilities_B500.RData")
pathway_file <- file.path(out_dir, "p1_pathway_P1_P2_Rpost_point.RData")
implant_file <- file.path(out_dir, "p1_implant_specific_revision_secondary.RData")
out_xlsx <- file.path(paper_dir, "p1_final_thesis_outputs.xlsx")
out_rdata <- file.path(paper_dir, "p1_final_thesis_outputs.RData")

for (f in c(analysis_file, main_boot_file, overlap_file)) if (!file.exists(f)) stop("required file not found: ", f)

load_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}

analysis_env <- load_env(analysis_file)
main_env <- load_env(main_boot_file)
ow_env <- load_env(overlap_file)

required_main <- c("probability_ci", "contrast_ci", "bootstrap_validation", "origin_success")
if (!all(vapply(required_main, exists, logical(1), envir = main_env))) stop("main B500 RData is missing expected objects")
required_ow <- c("bootstrap_probability_ci", "bootstrap_contrast_ci", "success_by_method", "qc")
if (!all(vapply(required_ow, exists, logical(1), envir = ow_env))) stop("overlap B500 RData is missing expected objects")
if (!exists("p1_ms_long", envir = analysis_env)) stop("p1_ms_long not found in analysis RData")

probability_ci <- main_env$probability_ci
contrast_ci <- main_env$contrast_ci
bootstrap_validation <- main_env$bootstrap_validation
origin_success <- main_env$origin_success
overlap_probability_ci <- ow_env$bootstrap_probability_ci
overlap_contrast_ci <- ow_env$bootstrap_contrast_ci
overlap_success <- ow_env$success_by_method
overlap_qc <- ow_env$qc
p1_ms_long <- as_tibble(analysis_env$p1_ms_long)

transition_from_origin_state <- function(origin, state) {
  case_when(
    origin == "P1" & state == "P2" ~ "P1_to_P2",
    origin == "P1" & state == "R_pre" ~ "P1_to_Rpre",
    origin == "P1" & state == "D" ~ "P1_to_D",
    origin == "P2" & state == "R_post" ~ "P2_to_Rpost",
    origin == "P2" & state == "D" ~ "P2_to_D",
    origin == "R_pre" & state == "P2_after_Rpre" ~ "Rpre_to_P2",
    origin == "R_pre" & state == "D" ~ "Rpre_to_D",
    TRUE ~ paste(origin, state, sep = "_to_")
  )
}

transition_label <- function(x) {
  recode(x,
         P1_to_P2 = "P1 -> P2",
         P1_to_Rpre = "P1 -> Rpre",
         P1_to_D = "P1 -> death",
         P2_to_Rpost = "P2 -> Rpost",
         P2_to_D = "P2 -> death",
         Rpre_to_P2 = "Rpre -> P2 after Rpre",
         Rpre_to_D = "Rpre -> death",
         .default = x
  )
}

fmt_prob <- function(point, low, high) sprintf("%.1f%% [%.1f%%, %.1f%%]", 100 * point, 100 * low, 100 * high)
fmt_rd <- function(point, low, high) sprintf("%.1f pp [%.1f, %.1f]", 100 * point, 100 * low, 100 * high)

main_prob_long <- probability_ci %>%
  mutate(
    origin = as.character(origin), state = as.character(state), group = as.character(group),
    transition = transition_from_origin_state(origin, state), transition_label = transition_label(transition),
    probability_ci = fmt_prob(point_estimate, lower_95, upper_95)
  )

main_prob_table <- main_prob_long %>%
  select(origin, transition, transition_label, time_years, group, probability_ci) %>%
  pivot_wider(names_from = group, values_from = probability_ci) %>%
  left_join(
    contrast_ci %>% mutate(origin = as.character(origin), state = as.character(state), transition = transition_from_origin_state(origin, state)) %>%
      transmute(origin, transition, time_years, difference = fmt_rd(point_risk_difference, rd_lower_95, rd_upper_95)),
    by = c("origin", "transition", "time_years")
  ) %>%
  arrange(factor(origin, levels = c("P1", "P2", "R_pre")), transition, time_years)

main_prob_table_5_10_15 <- main_prob_table %>% filter(time_years %in% c(5, 10, 15))
main_rd_numeric <- contrast_ci %>%
  mutate(origin = as.character(origin), state = as.character(state), transition = transition_from_origin_state(origin, state), transition_label = transition_label(transition)) %>%
  transmute(origin, transition, transition_label, time_years, rd_percentage_points = point_rd_percentage_points, lower_95_pp = rd_lower_95_percentage_points, upper_95_pp = rd_upper_95_percentage_points, n_boot_success)

ow_prob_table <- overlap_probability_ci %>%
  filter(method == "overlap_weighted_AJ") %>%
  mutate(probability_ci = fmt_prob(point_estimate, ci_low, ci_high), transition_label = transition_label(transition)) %>%
  select(transition, transition_label, time_years, group, probability_ci) %>%
  pivot_wider(names_from = group, values_from = probability_ci) %>%
  left_join(
    overlap_contrast_ci %>% filter(method == "overlap_weighted_AJ") %>% transmute(transition, time_years, difference = fmt_rd(point_estimate, ci_low, ci_high)),
    by = c("transition", "time_years")
  ) %>% arrange(transition, time_years)

ow_rd_numeric <- overlap_contrast_ci %>%
  filter(method == "overlap_weighted_AJ") %>%
  mutate(transition_label = transition_label(transition)) %>%
  transmute(transition, transition_label, time_years, rd_percentage_points = point_pp, lower_95_pp = ci_low_pp, upper_95_pp = ci_high_pp, n_boot)

rep_transition <- c(P1 = "P1_to_P2", P2 = "P2_to_Rpost", R_pre = "Rpre_to_P2", R_post = "Rpost_to_D")
risk_times <- c(0, 1, 3, 5, 10, 15)

risk_base <- p1_ms_long %>%
  mutate(trans_chr = as.character(trans), group = as.character(hip)) %>%
  filter(trans_chr %in% unname(rep_transition)) %>%
  mutate(origin = names(rep_transition)[match(trans_chr, rep_transition)]) %>%
  distinct(origin, patient_side_id, .keep_all = TRUE) %>%
  filter(group %in% c("first_knee", "first_hip"), is.finite(Tstop), Tstop > 0)

number_at_risk <- risk_base %>%
  crossing(time_years = risk_times) %>%
  group_by(origin, group, time_years) %>%
  summarise(n_at_risk = sum(Tstop >= time_years), .groups = "drop") %>%
  arrange(factor(origin, levels = c("P1", "P2", "R_pre", "R_post")), group, time_years)

number_at_risk_wide <- number_at_risk %>%
  mutate(time = paste0("year_", time_years)) %>%
  select(origin, group, time, n_at_risk) %>%
  pivot_wider(names_from = time, values_from = n_at_risk)

implant_cif_report <- tibble()
implant_revision_difference <- tibble()
implant_cox_results <- tibble()
implant_ph_results <- tibble()
if (file.exists(implant_file)) {
  ie <- load_env(implant_file)
  if (exists("cif_report", envir = ie)) implant_cif_report <- ie$cif_report
  if (exists("revision_difference", envir = ie)) implant_revision_difference <- ie$revision_difference
  if (exists("cox_results", envir = ie)) implant_cox_results <- ie$cox_results
  if (exists("ph_results", envir = ie)) implant_ph_results <- ie$ph_results
}

pathway_report <- tibble()
pathway_contrast <- tibble()
if (file.exists(pathway_file)) {
  pe <- load_env(pathway_file)
  if (exists("pathway_report", envir = pe)) pathway_report <- pe$pathway_report
  if (exists("pathway_contrast", envir = pe)) pathway_contrast <- pe$pathway_contrast
}

main_complete <- suppressWarnings(as.numeric(bootstrap_validation$value[bootstrap_validation$check == "completed_replicates"]))
ow_complete <- overlap_success %>% filter(method == "overlap_weighted_AJ") %>% pull(successful_replicates)
qc <- tibble(
  check = c("main bootstrap completed 500", "all main origins have 500 successful replicates", "overlap bootstrap has 500 successful replicates", "overlap QC all PASS", "P1 risk set equals 9703", "implant-specific secondary output available", "pathway point output available"),
  value = c(main_complete == 500, all(origin_success$successful_replicates == 500), length(ow_complete) == 1 && ow_complete == 500, all(overlap_qc$status == "PASS"), sum(risk_base$origin == "P1") == 9703, nrow(implant_cif_report) > 0, nrow(pathway_report) > 0),
  status = c(if_else(main_complete == 500, "PASS", "CHECK"), if_else(all(origin_success$successful_replicates == 500), "PASS", "CHECK"), if_else(length(ow_complete) == 1 && ow_complete == 500, "PASS", "CHECK"), if_else(all(overlap_qc$status == "PASS"), "PASS", "CHECK"), if_else(sum(risk_base$origin == "P1") == 9703, "PASS", "CHECK"), if_else(nrow(implant_cif_report) > 0, "PASS", "INFO"), if_else(nrow(pathway_report) > 0, "PASS", "INFO"))
)

notes <- tibble(note = c(
  "Main uncertainty uses the frozen patient-level bootstrap with B=500.",
  "Overlap-weighted uncertainty also uses B=500 and re-estimates the propensity score in every patient-level bootstrap replicate.",
  "Overlap-weighted results are a measured-confounding sensitivity analysis and are not presented as causal estimates.",
  "Use 5 and 10 years as the main compact reporting horizons; 1 and 3 years are useful for early non-proportional effects and 15 years is a secondary long-term estimate.",
  "Rpre-origin long-term estimates should be interpreted cautiously because the risk set is sparse.",
  "The pathway output, when present, currently contains point estimates only and should not be shown with confidence intervals unless a dedicated bootstrap is later justified."
))

wb <- createWorkbook()
tabs <- list(notes = notes, qc = qc, main_probabilities = main_prob_table, main_prob_5_10_15 = main_prob_table_5_10_15, main_rd_numeric = main_rd_numeric, overlap_probabilities = ow_prob_table, overlap_rd_numeric = ow_rd_numeric, number_at_risk = number_at_risk, number_at_risk_wide = number_at_risk_wide, implant_cif = implant_cif_report, implant_difference = implant_revision_difference, implant_cox = implant_cox_results, implant_ph = implant_ph_results, pathway_point = pathway_report, pathway_contrast = pathway_contrast)
for (nm in names(tabs)) {
  addWorksheet(wb, substr(nm, 1, 31))
  writeData(wb, substr(nm, 1, 31), tabs[[nm]])
  setColWidths(wb, substr(nm, 1, 31), cols = 1:50, widths = "auto")
}
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
save(main_prob_table, main_prob_table_5_10_15, main_rd_numeric, ow_prob_table, ow_rd_numeric, number_at_risk, number_at_risk_wide, implant_cif_report, implant_revision_difference, implant_cox_results, implant_ph_results, pathway_report, pathway_contrast, qc, notes, file = out_rdata)

cat("\n15 complete\n")
print(qc)
cat(out_xlsx, "\n")