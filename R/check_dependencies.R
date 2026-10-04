required_packages <- c(
  "dplyr", "tidyr", "stringr", "purrr", "readxl", "lubridate",
  "survival", "broom", "openxlsx", "ggplot2", "scales", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}
