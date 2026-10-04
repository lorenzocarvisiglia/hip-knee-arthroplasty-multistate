#example project configuration

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

data_dir <- file.path(project_dir, "data", "raw")
output_dir <- file.path(project_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

hip_file <- file.path(data_dir, "hip.xlsx")
knee_file <- file.path(data_dir, "knee.xlsx")
registry_file <- file.path(data_dir, "registry.xlsx")

hip_sheet <- "Esporta foglio di lavoro"
knee_sheet <- "Esporta foglio di lavoro"
interventions_sheet <- "INTERVENTI"
deaths_sheet <- "CON DECESSO COME EVENTO"

admin_censor_date <- as.Date("2021-12-31")

clean_data_file <- file.path(output_dir, "p1_multistate_option3_clean_data.RData")
analysis_file <- file.path(output_dir, "p1_multistate_analysis_objects.RData")
