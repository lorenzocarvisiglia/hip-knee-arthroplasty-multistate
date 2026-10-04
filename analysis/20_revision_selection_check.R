source("config.R")

library(dplyr)
library(stringr)
library(readxl)
library(lubridate)

administrative_censoring <- admin_censor_date

anca_file <- hip_file
ginocchio_file <- knee_file
dati_file <- registry_file

norm_txt <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\u00a0", " ") %>%
    str_squish() %>%
    str_to_lower()
}

norm_side <- function(x) {
  y <- norm_txt(x)
  case_when(
    y %in% c("destro", "dx", "right", "r") ~ "right",
    y %in% c("sinistro", "sx", "left", "l") ~ "left",
    TRUE ~ NA_character_
  )
}

to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  
  x <- str_squish(as.character(x))
  
  out <- suppressWarnings(as.Date(as.numeric(x), origin = "1899-12-30"))
  
  idx <- is.na(out) & !is.na(x)
  out[idx] <- suppressWarnings(ymd(x[idx]))
  
  idx <- is.na(out) & !is.na(x)
  out[idx] <- suppressWarnings(dmy(x[idx]))
  
  idx <- is.na(out) & !is.na(x)
  out[idx] <- suppressWarnings(mdy(x[idx]))
  
  out
}

anca <- read_excel(anca_file, sheet = hip_sheet)
ginocchio <- read_excel(ginocchio_file, sheet = knee_sheet)
interventi <- read_excel(dati_file, sheet = interventions_sheet)

hip_primary <- anca %>%
  transmute(
    CODPAT = as.character(CODPAT),
    side = norm_side(D_LATO),
    primary_date = to_date(DATAINT),
    residence = norm_txt(R_RESIDENZA),
    tipoint = norm_txt(D_TIPOINT),
    diagnosis_code = as.character(DIAGNOSI),
    diagnosis_label = norm_txt(D_DIAGNOSI)
  ) %>%
  filter(
    residence == "emilia romagna",
    str_detect(tipoint, "protesi primaria"),
    diagnosis_code == "1" | str_detect(diagnosis_label, "^1\\s*coxartrosi"),
    !is.na(CODPAT),
    !is.na(side),
    !is.na(primary_date)
  ) %>%
  arrange(CODPAT, side, primary_date) %>%
  group_by(CODPAT, side) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(joint = "hip")

knee_primary <- ginocchio %>%
  transmute(
    CODPAT = as.character(CODPAT),
    side = norm_side(D_LATO),
    primary_date = to_date(DATAINT),
    residence = norm_txt(R_RESIDENZA),
    tipoint = norm_txt(D_TIPOINT),
    diagnosis_code = as.character(DIAGNOSI),
    diagnosis_label = norm_txt(D_DIAGNOSI),
    prosthesis_code = as.character(TIPOPROT),
    prosthesis_label = norm_txt(D_TIPOPROT)
  ) %>%
  filter(
    residence == "emilia romagna",
    str_detect(tipoint, "protesi primaria"),
    diagnosis_code %in% c("1", "3", "10", "11") |
      str_detect(diagnosis_label, "artrosi primaria|deformit"),
    prosthesis_code %in% c("2", "3") |
      str_detect(prosthesis_label, "bicompartimentale|tricompartimentale"),
    !is.na(CODPAT),
    !is.na(side),
    !is.na(primary_date)
  ) %>%
  arrange(CODPAT, side, primary_date) %>%
  group_by(CODPAT, side) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(joint = "knee")

primaries <- bind_rows(hip_primary, knee_primary)

revisions <- interventi %>%
  transmute(
    CODPAT = as.character(CODPAT),
    joint = case_when(
      norm_txt(INTERVENTO_TIPO) == "anca" ~ "hip",
      norm_txt(INTERVENTO_TIPO) == "ginocchio" ~ "knee",
      TRUE ~ NA_character_
    ),
    side = norm_side(D_LATO),
    event_date = to_date(DATAINT),
    event_type = norm_txt(D_TIPOINT)
  ) %>%
  filter(
    str_detect(event_type, "reimpianto|espianto"),
    !is.na(CODPAT),
    !is.na(joint),
    !is.na(side),
    !is.na(event_date),
    event_date <= administrative_censoring
  )

check <- primaries %>%
  left_join(revisions, by = c("CODPAT", "side", "joint")) %>%
  group_by(CODPAT, side, joint, primary_date) %>%
  summarise(
    first_revision_any = if (
      all(is.na(event_date))
    ) as.Date(NA) else min(event_date, na.rm = TRUE),
    
    first_revision_after_primary = if (
      all(is.na(event_date[event_date > primary_date]))
    ) as.Date(NA) else min(event_date[event_date > primary_date], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    old_method = if_else(
      !is.na(first_revision_any) & first_revision_any > primary_date,
      first_revision_any,
      as.Date(NA)
    ),
    differs = old_method != first_revision_after_primary |
      xor(is.na(old_method), is.na(first_revision_after_primary))
  )

cat("\nnumero articolazioni eleggibili:", nrow(check), "\n")

cat(
  "revisioni post-primary con metodo corretto:",
  sum(!is.na(check$first_revision_after_primary)),
  "\n"
)

cat(
  "casi in cui il metodo attuale perde o cambia una revisione:",
  sum(check$differs, na.rm = TRUE),
  "\n\n"
)

print(
  check %>%
    filter(differs) %>%
    arrange(CODPAT, side, joint)
)