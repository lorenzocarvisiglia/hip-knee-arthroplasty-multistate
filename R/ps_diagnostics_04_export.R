# ============================================================
# salva Excel
# ============================================================

wb <- openxlsx::createWorkbook()

header_style <- safe_create_style(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "bottom",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

body_style <- safe_create_style(
  valign = "top",
  wrapText = TRUE
)

note_style <- safe_create_style(
  textDecoration = "italic",
  fontColour = "#555555",
  wrapText = TRUE
)

title_style <- safe_create_style(
  textDecoration = "bold",
  fontSize = 13
)

add_sheet <- function(wb, sheet_name, dat, title = NULL, note = NULL) {
  openxlsx::addWorksheet(wb, sheet_name)
  
  start_row <- 1
  
  if (!is.null(title)) {
    openxlsx::writeData(wb, sheet_name, title, startRow = start_row, startCol = 1)
    openxlsx::addStyle(
      wb,
      sheet_name,
      title_style,
      rows = start_row,
      cols = 1,
      gridExpand = TRUE
    )
    start_row <- start_row + 2
  }
  
  if (!is.null(note)) {
    openxlsx::writeData(wb, sheet_name, note, startRow = start_row, startCol = 1)
    openxlsx::addStyle(
      wb,
      sheet_name,
      note_style,
      rows = start_row,
      cols = 1,
      gridExpand = TRUE
    )
    start_row <- start_row + 2
  }
  
  openxlsx::writeData(
    wb,
    sheet_name,
    dat,
    startRow = start_row,
    startCol = 1,
    withFilter = TRUE
  )
  
  if (ncol(dat) > 0) {
    openxlsx::addStyle(
      wb,
      sheet_name,
      header_style,
      rows = start_row,
      cols = seq_len(ncol(dat)),
      gridExpand = TRUE
    )
  }
  
  if (nrow(dat) > 0 && ncol(dat) > 0) {
    openxlsx::addStyle(
      wb,
      sheet_name,
      body_style,
      rows = (start_row + 1):(start_row + nrow(dat)),
      cols = seq_len(ncol(dat)),
      gridExpand = TRUE
    )
  }
  
  openxlsx::freezePane(wb, sheet_name, firstActiveRow = start_row + 1)
  openxlsx::setColWidths(wb, sheet_name, cols = 1:max(1, ncol(dat)), widths = "auto")
}

add_sheet(
  wb,
  "00_QC",
  qc,
  title = "QC propensity-score diagnostics",
  note = "Check counts, PS model convergence, weight range, ESS, and weighted SMD."
)

add_sheet(
  wb,
  "01_methods_text",
  methods_text,
  title = "Methods text",
  note = "Text blocks for manuscript Methods/Supplement."
)

add_sheet(
  wb,
  "02_column_mapping",
  column_mapping,
  title = "Detected source columns",
  note = "Column names selected automatically from the source RData."
)

add_sheet(
  wb,
  "03_covariates_used",
  covariates_used,
  title = "Covariates considered and used in the PS model",
  note = "Variables used only if available and sufficiently variable."
)

add_sheet(
  wb,
  "04_ps_model_coefficients",
  ps_model_coef,
  title = "Propensity-score model coefficients",
  note = "Logistic regression for first hip vs first knee."
)

add_sheet(
  wb,
  "05_ps_summary_by_group",
  ps_summary_by_group,
  title = "Propensity-score summary by implant-order group",
  note = "Propensity score is probability of first hip."
)

add_sheet(
  wb,
  "06_weight_summary_by_group",
  weight_summary_by_group,
  title = "Overlap-weight summary by implant-order group",
  note = "First hip receives 1 - e(X); first knee receives e(X)."
)

add_sheet(
  wb,
  "07_ESS",
  ess_summary,
  title = "Effective sample size",
  note = "ESS = (sum weights)^2 / sum weights^2."
)

add_sheet(
  wb,
  "08_support_summary",
  support_summary,
  title = "Propensity-score support by group",
  note = "Used to assess empirical overlap."
)

add_sheet(
  wb,
  "09_common_support",
  common_support,
  title = "Common support summary",
  note = "Common support is based on overlapping min-max PS ranges."
)

add_sheet(
  wb,
  "10_balance_long",
  balance_long,
  title = "Standardized mean differences, long format",
  note = "Includes unweighted and overlap-weighted balance."
)

add_sheet(
  wb,
  "11_balance_wide",
  balance_wide,
  title = "Standardized mean differences, wide format",
  note = "Useful supplement table."
)

add_sheet(
  wb,
  "12_balance_by_covariate",
  balance_summary_by_covariate,
  title = "Balance summary by covariate",
  note = "Maximum and mean absolute SMD by covariate."
)

add_sheet(
  wb,
  "13_balance_overall",
  balance_overall_summary,
  title = "Overall balance summary",
  note = "Main QC table for balance before/after overlap weighting."
)

add_sheet(
  wb,
  "14_love_plot_data",
  love_plot_data,
  title = "Love plot data",
  note = "Top terms by imbalance used in the Love plot."
)

add_sheet(
  wb,
  "15_analysis_dataset",
  ps_data,
  title = "PS analysis dataset",
  note = "One row per defined-order patient-side unit."
)

add_sheet(
  wb,
  "16_figures_manifest",
  figures_manifest,
  title = "Figures created",
  note = "Paths to diagnostic figures."
)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# ============================================================
# note
# ============================================================

notes <- c(
  "Script 17 completed.",
  "",
  paste0("Source RData: ", normalizePath(source_path, winslash = "/", mustWork = FALSE)),
  paste0("Source object: ", source_object),
  paste0("Output Excel: ", normalizePath(out_xlsx, winslash = "/", mustWork = FALSE)),
  paste0("Output RData: ", normalizePath(out_rdata, winslash = "/", mustWork = FALSE)),
  paste0("Figures directory: ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE)),
  "",
  "Main checks:",
  paste0("n first hip = ", n_first_hip),
  paste0("n first knee = ", n_first_knee),
  paste0("overall ESS overlap = ", fmt_num(ess_summary$ess_overlap[ess_summary$group == 'Overall'], 1)),
  paste0("max abs SMD unweighted = ", fmt_num(max_smd_unweighted, 3)),
  paste0("max abs SMD overlap weighted = ", fmt_num(max_smd_weighted, 3)),
  paste0("weighted SMD terms > 0.10 = ", n_smd_gt_010_weighted),
  "",
  "Interpretation:",
  "This script provides diagnostic evidence for the overlap-weighted sensitivity analysis.",
  "It should be reported in the supplement rather than as a main clinical result."
)

writeLines(notes, con = out_txt, useBytes = TRUE)

# ============================================================
# console
# ============================================================

message("Script 17 completato.")
message("Excel salvato in: ", normalizePath(out_xlsx, winslash = "/", mustWork = FALSE))
message("RData salvato in: ", normalizePath(out_rdata, winslash = "/", mustWork = FALSE))
message("Figure salvate in: ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE))

print(qc)
print(balance_overall_summary)
print(ess_summary)
print(weight_summary_by_group)