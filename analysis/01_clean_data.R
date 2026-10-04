#pacchetti
source("config.R")

pkgs <- c("dplyr", "tidyr", "stringr", "purrr", "readxl", "lubridate", "openxlsx")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("installa questi pacchetti prima di eseguire lo script: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(pkgs, library, character.only = TRUE))

#parameters
administrative_censoring <- admin_censor_date
input_dir <- data_dir
out_dir <- output_dir

required_inputs <- c(hip_file, knee_file, registry_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing input file(s): ", paste(missing_inputs, collapse = ", "))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

anca_file <- hip_file
ginocchio_file <- knee_file
dati_file <- registry_file

#funzioni
pick_col <- function(dat, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) {
    if (required) stop("colonna mancante. cercate: ", paste(candidates, collapse = ", "))
    return(NULL)
  }
  hit[1]
}

get_col <- function(dat, candidates, default = NA) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) return(rep(default, nrow(dat)))
  dat[[hit[1]]]
}

as_num <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[is.infinite(out)] <- NA_real_
  out
}

to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x_chr <- stringr::str_squish(as.character(x))
  x_chr[x_chr %in% c("", "-", "NA", "NaN", "NULL")] <- NA_character_
  out <- suppressWarnings(as.Date(as.numeric(x_chr), origin = "1899-12-30"))
  idx <- is.na(out) & !is.na(x_chr)
  if (any(idx)) out[idx] <- suppressWarnings(lubridate::ymd(x_chr[idx]))
  idx <- is.na(out) & !is.na(x_chr)
  if (any(idx)) out[idx] <- suppressWarnings(lubridate::dmy(x_chr[idx]))
  idx <- is.na(out) & !is.na(x_chr)
  if (any(idx)) out[idx] <- suppressWarnings(lubridate::mdy(x_chr[idx]))
  out
}

years_between <- function(end, start) {
  as.numeric(difftime(as.Date(end), as.Date(start), units = "days")) / 365.25
}

norm_txt <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\u00a0", " ")
  x <- stringr::str_squish(x)
  stringr::str_to_lower(x)
}

norm_side <- function(x) {
  y <- norm_txt(x)
  dplyr::case_when(
    y %in% c("destro", "dx", "right", "r") ~ "right",
    y %in% c("sinistro", "sx", "left", "l") ~ "left",
    TRUE ~ NA_character_
  )
}

norm_sex <- function(x) {
  y <- norm_txt(x)
  dplyr::case_when(
    y %in% c("f", "female", "femmina", "2") ~ "female",
    y %in% c("m", "male", "maschio", "1") ~ "male",
    TRUE ~ NA_character_
  )
}

bmi_class <- function(x) {
  z <- as_num(x)
  dplyr::case_when(
    is.na(z) ~ NA_character_,
    z < 25 ~ "normal_underweight",
    z < 30 ~ "overweight",
    z >= 30 ~ "obese"
  )
}

clean_bmi_cat <- function(raw_cat, raw_bmi) {
  y <- norm_txt(raw_cat)
  out <- dplyr::case_when(
    stringr::str_detect(y, "under|normal|normo|sottopeso") ~ "normal_underweight",
    stringr::str_detect(y, "over|sovrappeso") ~ "overweight",
    stringr::str_detect(y, "obes") ~ "obese",
    TRUE ~ NA_character_
  )
  dplyr::coalesce(out, bmi_class(raw_bmi), "missing")
}

min_date_na <- function(...) {
  x <- c(...)
  x <- as.Date(x)
  if (all(is.na(x))) return(as.Date(NA))
  min(x, na.rm = TRUE)
}

max_date_na <- function(...) {
  x <- c(...)
  x <- as.Date(x)
  if (all(is.na(x))) return(as.Date(NA))
  max(x, na.rm = TRUE)
}

first_date_after <- function(event_date, origin_date) {
  out <- as.Date(event_date)
  origin_date <- as.Date(origin_date)
  out[is.na(origin_date) | is.na(out) | out <= origin_date] <- as.Date(NA)
  out
}

date_min2 <- function(a, b) {
  aa <- as.numeric(as.Date(a))
  bb <- as.numeric(as.Date(b))
  out <- pmin(aa, bb, na.rm = TRUE)
  out[is.infinite(out)] <- NA_real_
  as.Date(out, origin = "1970-01-01")
}

date_max2_if_both <- function(a, b) {
  aa <- as.numeric(as.Date(a))
  bb <- as.numeric(as.Date(b))
  out <- pmax(aa, bb, na.rm = FALSE)
  out[is.infinite(out)] <- NA_real_
  as.Date(out, origin = "1970-01-01")
}

first_date_by_row <- function(...) {
  dots <- list(...)
  mat <- do.call(cbind, lapply(dots, as.Date))
  out <- apply(mat, 1, function(z) {
    z <- as.Date(z, origin = "1970-01-01")
    z <- z[!is.na(z)]
    if (length(z) == 0) as.Date(NA) else min(z)
  })
  as.Date(out, origin = "1970-01-01")
}

#lettura file
anca_raw <- readxl::read_excel(anca_file, sheet = hip_sheet)
ginocchio_raw <- readxl::read_excel(ginocchio_file, sheet = knee_sheet)
interventi_raw <- readxl::read_excel(dati_file, sheet = interventions_sheet)
decessi_raw <- readxl::read_excel(dati_file, sheet = deaths_sheet)

#primarie eleggibili
make_hip <- function(dat) {
  dat %>%
    transmute(
      CODPAT = as.character(.data[[pick_col(dat, c("CODPAT"))]]),
      side = norm_side(.data[[pick_col(dat, c("D_LATO"))]]),
      side_raw = as.character(.data[[pick_col(dat, c("D_LATO"))]]),
      hip_primary_date = to_date(.data[[pick_col(dat, c("DATAINT"))]]),
      sex_hip = norm_sex(get_col(dat, c("D_SESSO", "SESSO"))),
      age_hip = as_num(get_col(dat, c("AGE"))),
      bmi_hip = as_num(get_col(dat, c("BMI"))),
      bmi_cat_hip = clean_bmi_cat(get_col(dat, c("CL_BMI")), get_col(dat, c("BMI"))),
      residence_hip = as.character(get_col(dat, c("R_RESIDENZA"))),
      diagnosis_hip_code = as.character(get_col(dat, c("DIAGNOSI"))),
      diagnosis_hip_label = as.character(get_col(dat, c("D_DIAGNOSI"))),
      tipoint_hip = as.character(get_col(dat, c("D_TIPOINT"))),
      center_hip_code = as.character(get_col(dat, c("VREPARTI", "DVREPARTI"))),
      center_hip_label = as.character(get_col(dat, c("D_DVREPARTI", "DOISTAT", "D_OSPEDALE_FULL"))),
      center_hip_short = as.character(get_col(dat, c("DOISTAT", "D_DVREPARTI"))),
      visitnum_hip = as.character(get_col(dat, c("VISITNUM"))),
      visitnum_progr_hip = as.character(get_col(dat, c("VISITNUM_PROGR")))
    ) %>%
    mutate(
      resident_er_hip = norm_txt(residence_hip) == "emilia romagna",
      primary_hip = stringr::str_detect(norm_txt(tipoint_hip), "protesi primaria") | norm_txt(tipoint_hip) == "p protesi primaria",
      coxarthrosis_hip = diagnosis_hip_code == "1" | stringr::str_detect(norm_txt(diagnosis_hip_label), "^1\\s*coxartrosi")
    ) %>%
    filter(resident_er_hip, primary_hip, coxarthrosis_hip, !is.na(CODPAT), !is.na(side), !is.na(hip_primary_date)) %>%
    arrange(CODPAT, side, hip_primary_date) %>%
    group_by(CODPAT, side) %>%
    slice(1) %>%
    ungroup()
}

make_knee <- function(dat) {
  dat %>%
    transmute(
      CODPAT = as.character(.data[[pick_col(dat, c("CODPAT"))]]),
      side = norm_side(.data[[pick_col(dat, c("D_LATO"))]]),
      side_raw = as.character(.data[[pick_col(dat, c("D_LATO"))]]),
      knee_primary_date = to_date(.data[[pick_col(dat, c("DATAINT"))]]),
      sex_knee = norm_sex(get_col(dat, c("D_SESSO", "SESSO"))),
      age_knee = as_num(get_col(dat, c("AGE"))),
      bmi_knee = as_num(get_col(dat, c("BMI"))),
      bmi_cat_knee = clean_bmi_cat(get_col(dat, c("CL_BMI")), get_col(dat, c("BMI"))),
      residence_knee = as.character(get_col(dat, c("R_RESIDENZA"))),
      diagnosis_knee_code = as.character(get_col(dat, c("DIAGNOSI"))),
      diagnosis_knee_label = as.character(get_col(dat, c("D_DIAGNOSI"))),
      tipoint_knee = as.character(get_col(dat, c("D_TIPOINT"))),
      knee_prosthesis_code = as.character(get_col(dat, c("TIPOPROT"))),
      knee_prosthesis_label = as.character(get_col(dat, c("D_TIPOPROT"))),
      center_knee_code = as.character(get_col(dat, c("VREPARTI", "DVREPARTI"))),
      center_knee_label = as.character(get_col(dat, c("D_DVREPARTI", "DOISTAT"))),
      center_knee_short = as.character(get_col(dat, c("DOISTAT", "D_DVREPARTI"))),
      visitnum_knee = as.character(get_col(dat, c("VISITNUM"))),
      visitnum_progr_knee = as.character(get_col(dat, c("VISITNUM_PROGR")))
    ) %>%
    mutate(
      resident_er_knee = norm_txt(residence_knee) == "emilia romagna",
      primary_knee = stringr::str_detect(norm_txt(tipoint_knee), "protesi primaria") | norm_txt(tipoint_knee) == "p protesi primaria",
      diagnosis_ok_knee = diagnosis_knee_code %in% c("1", "3", "10", "11") |
        stringr::str_detect(norm_txt(diagnosis_knee_label), "artrosi primaria|deformit"),
      prosthesis_ok_knee = knee_prosthesis_code %in% c("2", "3") |
        stringr::str_detect(norm_txt(knee_prosthesis_label), "bicompartimentale|tricompartimentale")
    ) %>%
    filter(resident_er_knee, primary_knee, diagnosis_ok_knee, prosthesis_ok_knee, !is.na(CODPAT), !is.na(side), !is.na(knee_primary_date)) %>%
    arrange(CODPAT, side, knee_primary_date) %>%
    group_by(CODPAT, side) %>%
    slice(1) %>%
    ungroup()
}

hip_eligible <- make_hip(anca_raw)
knee_eligible <- make_knee(ginocchio_raw)

#eventi da dati.xlsx
interventions_clean <- interventi_raw %>%
  transmute(
    CODPAT = as.character(.data[[pick_col(interventi_raw, c("CODPAT"))]]),
    joint = dplyr::case_when(
      norm_txt(.data[[pick_col(interventi_raw, c("INTERVENTO_TIPO"))]]) == "anca" ~ "hip",
      norm_txt(.data[[pick_col(interventi_raw, c("INTERVENTO_TIPO"))]]) == "ginocchio" ~ "knee",
      TRUE ~ NA_character_
    ),
    side = norm_side(.data[[pick_col(interventi_raw, c("D_LATO"))]]),
    event_date = to_date(.data[[pick_col(interventi_raw, c("DATAINT"))]]),
    event_type_label = as.character(.data[[pick_col(interventi_raw, c("D_TIPOINT"))]]),
    visitnum = as.character(get_col(interventi_raw, c("VISITNUM"))),
    visitnum_progr = as.character(get_col(interventi_raw, c("VISITNUM_PROGR")))
  ) %>%
  mutate(
    event_type_norm = norm_txt(event_type_label),
    event_class = dplyr::case_when(
      stringr::str_detect(event_type_norm, "protesi primaria") ~ "primary",
      stringr::str_detect(event_type_norm, "reimpianto|espianto") ~ "revision_or_removal",
      TRUE ~ "other"
    )
  ) %>%
  filter(!is.na(CODPAT), !is.na(joint), !is.na(side), !is.na(event_date))

death_clean <- decessi_raw %>%
  transmute(
    CODPAT = as.character(.data[[pick_col(decessi_raw, c("CODPAT"))]]),
    event_label = as.character(.data[[pick_col(decessi_raw, c("Evento"))]]),
    death_date_raw = to_date(.data[[pick_col(decessi_raw, c("DATAINT o DATA DECESSO"))]])
  ) %>%
  filter(norm_txt(event_label) == "decesso", !is.na(death_date_raw)) %>%
  group_by(CODPAT) %>%
  summarise(death_date_raw = min(death_date_raw), .groups = "drop")

hip_revision_dates <- interventions_clean %>%
  filter(
    event_class == "revision_or_removal",
    joint == "hip",
    event_date <= administrative_censoring
  ) %>%
  inner_join(
    hip_eligible %>%
      select(CODPAT, side, hip_primary_date),
    by = c("CODPAT", "side")
  ) %>%
  filter(event_date > hip_primary_date) %>%
  group_by(CODPAT, side) %>%
  summarise(
    hip_revision_date_raw = min(event_date),
    .groups = "drop"
  )

knee_revision_dates <- interventions_clean %>%
  filter(
    event_class == "revision_or_removal",
    joint == "knee",
    event_date <= administrative_censoring
  ) %>%
  inner_join(
    knee_eligible %>%
      select(CODPAT, side, knee_primary_date),
    by = c("CODPAT", "side")
  ) %>%
  filter(event_date > knee_primary_date) %>%
  group_by(CODPAT, side) %>%
  summarise(
    knee_revision_date_raw = min(event_date),
    .groups = "drop"
  )

#coorte patient-side
cohort_p1_clean <- full_join(
  hip_eligible,
  knee_eligible,
  by = c("CODPAT", "side"),
  suffix = c("_hiprow", "_kneerow")
) %>%
  left_join(hip_revision_dates, by = c("CODPAT", "side")) %>%
  left_join(knee_revision_dates, by = c("CODPAT", "side")) %>%
  left_join(death_clean, by = "CODPAT") %>%
  mutate(
    patient_side_id = paste(CODPAT, side, sep = "_"),
    has_eligible_hip = !is.na(hip_primary_date),
    has_eligible_knee = !is.na(knee_primary_date),
    has_p2_raw = has_eligible_hip & has_eligible_knee,
    side_raw = dplyr::coalesce(side_raw_hiprow, side_raw_kneerow),
    p1_date = date_min2(hip_primary_date, knee_primary_date),
    p2_date_raw = date_max2_if_both(hip_primary_date, knee_primary_date),
    same_day_primaries = has_p2_raw & hip_primary_date == knee_primary_date,
    p1_joint = dplyr::case_when(
      has_eligible_hip & !has_eligible_knee ~ "hip",
      !has_eligible_hip & has_eligible_knee ~ "knee",
      has_eligible_hip & has_eligible_knee & hip_primary_date < knee_primary_date ~ "hip",
      has_eligible_hip & has_eligible_knee & knee_primary_date < hip_primary_date ~ "knee",
      has_eligible_hip & has_eligible_knee & hip_primary_date == knee_primary_date ~ "same_day",
      TRUE ~ NA_character_
    ),
    p2_joint = dplyr::case_when(
      !has_p2_raw ~ NA_character_,
      p1_joint == "hip" ~ "knee",
      p1_joint == "knee" ~ "hip",
      p1_joint == "same_day" ~ "same_day",
      TRUE ~ NA_character_
    ),
    hip = dplyr::case_when(
      p1_joint == "hip" ~ 1L,
      p1_joint == "knee" ~ 0L,
      TRUE ~ NA_integer_
    ),
    order_group = dplyr::case_when(
      p1_joint == "hip" ~ "first_hip",
      p1_joint == "knee" ~ "first_knee",
      p1_joint == "same_day" ~ "same_day_undefined",
      TRUE ~ NA_character_
    ),
    p2_date = if_else(has_p2_raw & p2_date_raw <= administrative_censoring, p2_date_raw, as.Date(NA)),
    has_p2_by_admin = !is.na(p2_date),
    t12 = if_else(has_p2_raw, years_between(p2_date_raw, p1_date), NA_real_),
    t12_by_admin = if_else(has_p2_by_admin, years_between(p2_date, p1_date), NA_real_),
    sex = dplyr::coalesce(sex_hip, sex_knee),
    age_p1 = dplyr::case_when(
      p1_joint == "hip" ~ age_hip,
      p1_joint == "knee" ~ age_knee,
      p1_joint == "same_day" ~ dplyr::coalesce(age_hip, age_knee),
      TRUE ~ NA_real_
    ),
    bmi_p1 = dplyr::case_when(
      p1_joint == "hip" ~ bmi_hip,
      p1_joint == "knee" ~ bmi_knee,
      p1_joint == "same_day" ~ dplyr::coalesce(bmi_hip, bmi_knee),
      TRUE ~ NA_real_
    ),
    bmi_p1_f = dplyr::case_when(
      p1_joint == "hip" ~ bmi_cat_hip,
      p1_joint == "knee" ~ bmi_cat_knee,
      p1_joint == "same_day" ~ dplyr::coalesce(bmi_cat_hip, bmi_cat_knee),
      TRUE ~ "missing"
    ),
    age_p2 = dplyr::case_when(
      p2_joint == "hip" ~ age_hip,
      p2_joint == "knee" ~ age_knee,
      p2_joint == "same_day" ~ dplyr::coalesce(age_hip, age_knee),
      TRUE ~ NA_real_
    ),
    bmi_p2 = dplyr::case_when(
      p2_joint == "hip" ~ bmi_hip,
      p2_joint == "knee" ~ bmi_knee,
      p2_joint == "same_day" ~ dplyr::coalesce(bmi_hip, bmi_knee),
      TRUE ~ NA_real_
    ),
    bmi_p2_f = dplyr::case_when(
      p2_joint == "hip" ~ bmi_cat_hip,
      p2_joint == "knee" ~ bmi_cat_knee,
      p2_joint == "same_day" ~ dplyr::coalesce(bmi_cat_hip, bmi_cat_knee),
      TRUE ~ NA_character_
    ),
    center_p1_code = dplyr::case_when(
      p1_joint == "hip" ~ center_hip_code,
      p1_joint == "knee" ~ center_knee_code,
      p1_joint == "same_day" ~ dplyr::coalesce(center_hip_code, center_knee_code),
      TRUE ~ NA_character_
    ),
    center_p1_label = dplyr::case_when(
      p1_joint == "hip" ~ center_hip_label,
      p1_joint == "knee" ~ center_knee_label,
      p1_joint == "same_day" ~ dplyr::coalesce(center_hip_label, center_knee_label),
      TRUE ~ NA_character_
    ),
    center_p1_short = dplyr::case_when(
      p1_joint == "hip" ~ center_hip_short,
      p1_joint == "knee" ~ center_knee_short,
      p1_joint == "same_day" ~ dplyr::coalesce(center_hip_short, center_knee_short),
      TRUE ~ NA_character_
    ),
    center_p2_code = dplyr::case_when(
      p2_joint == "hip" ~ center_hip_code,
      p2_joint == "knee" ~ center_knee_code,
      p2_joint == "same_day" ~ dplyr::coalesce(center_hip_code, center_knee_code),
      TRUE ~ NA_character_
    ),
    center_p2_label = dplyr::case_when(
      p2_joint == "hip" ~ center_hip_label,
      p2_joint == "knee" ~ center_knee_label,
      p2_joint == "same_day" ~ dplyr::coalesce(center_hip_label, center_knee_label),
      TRUE ~ NA_character_
    ),
    hip_revision_date = first_date_after(hip_revision_date_raw, hip_primary_date),
    knee_revision_date = first_date_after(knee_revision_date_raw, knee_primary_date),
    hip_revision_date = if_else(!is.na(hip_revision_date) & hip_revision_date <= administrative_censoring, hip_revision_date, as.Date(NA)),
    knee_revision_date = if_else(!is.na(knee_revision_date) & knee_revision_date <= administrative_censoring, knee_revision_date, as.Date(NA)),
    first_revision_after_p1_date = first_date_by_row(hip_revision_date, knee_revision_date),
    death_date = if_else(!is.na(death_date_raw) & death_date_raw <= administrative_censoring, death_date_raw, as.Date(NA)),
    censor_date = administrative_censoring,
    first_event_date_from_p1 = first_date_by_row(p2_date, first_revision_after_p1_date, death_date, censor_date),
    first_event_from_p1 = dplyr::case_when(
      !is.na(p2_date) & p2_date == first_event_date_from_p1 ~ "p2",
      !is.na(first_revision_after_p1_date) & first_revision_after_p1_date == first_event_date_from_p1 ~ "revision",
      !is.na(death_date) & death_date == first_event_date_from_p1 ~ "death",
      !is.na(censor_date) & censor_date == first_event_date_from_p1 ~ "censor",
      TRUE ~ NA_character_
    ),
    time_first_from_p1 = years_between(first_event_date_from_p1, p1_date),
    revision_before_p2 = !is.na(first_revision_after_p1_date) & has_p2_raw & first_revision_after_p1_date < p2_date_raw,
    p2_before_revision = has_p2_raw & (is.na(first_revision_after_p1_date) | p2_date_raw < first_revision_after_p1_date),
    death_before_p2 = !is.na(death_date) & has_p2_raw & death_date < p2_date_raw,
    censored_before_p2 = !has_p2_by_admin & is.na(first_revision_after_p1_date) & is.na(death_date)
  ) %>%
  filter(!is.na(p1_date), p1_date <= administrative_censoring) %>%
  mutate(
    bmi_p1_f = factor(bmi_p1_f, levels = c("normal_underweight", "overweight", "obese", "missing")),
    bmi_p2_f = factor(bmi_p2_f, levels = c("normal_underweight", "overweight", "obese", "missing")),
    sex = factor(sex, levels = c("male", "female")),
    order_group = factor(order_group, levels = c("first_knee", "first_hip", "same_day_undefined")),
    p1_joint = factor(p1_joint, levels = c("knee", "hip", "same_day")),
    p2_joint = factor(p2_joint, levels = c("hip", "knee", "same_day"))
  ) %>%
  select(
    patient_side_id, CODPAT, side, side_raw, sex,
    p1_date, p1_joint, p2_date, p2_date_raw, p2_joint, has_p2_raw, has_p2_by_admin,
    same_day_primaries, order_group, hip, t12, t12_by_admin,
    age_p1, bmi_p1, bmi_p1_f, age_p2, bmi_p2, bmi_p2_f,
    hip_primary_date, knee_primary_date,
    diagnosis_hip_code, diagnosis_hip_label,
    diagnosis_knee_code, diagnosis_knee_label,
    knee_prosthesis_code, knee_prosthesis_label,
    center_p1_code, center_p1_label, center_p1_short,
    center_p2_code, center_p2_label,
    center_hip_code, center_hip_label, center_hip_short,
    center_knee_code, center_knee_label, center_knee_short,
    hip_revision_date, knee_revision_date, first_revision_after_p1_date,
    death_date, death_date_raw, censor_date,
    first_event_date_from_p1, first_event_from_p1, time_first_from_p1,
    revision_before_p2, p2_before_revision, death_before_p2, censored_before_p2,
    everything()
  )

#storia eventi pulita
event_history_clean <- interventions_clean %>%
  semi_join(cohort_p1_clean %>% select(CODPAT, side), by = c("CODPAT", "side")) %>%
  left_join(cohort_p1_clean %>% select(patient_side_id, CODPAT, side, p1_date, p2_date_raw, has_p2_raw), by = c("CODPAT", "side")) %>%
  mutate(
    time_from_p1 = years_between(event_date, p1_date),
    period = dplyr::case_when(
      is.na(p1_date) ~ NA_character_,
      event_date < p1_date ~ "before_p1",
      has_p2_raw & !is.na(p2_date_raw) & event_date < p2_date_raw ~ "p1_to_p2",
      has_p2_raw & !is.na(p2_date_raw) & event_date >= p2_date_raw ~ "after_p2",
      !has_p2_raw & event_date >= p1_date ~ "after_p1_no_p2",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(patient_side_id, event_date, joint, event_class)

#controlli qualità
qc_overall <- tibble::tibble(
  quantity = c(
    "eligible hip primaries",
    "eligible knee primaries",
    "patient-side units in union cohort",
    "patients",
    "only eligible hip primary",
    "only eligible knee primary",
    "eligible hip and knee primaries",
    "same-day eligible primaries",
    "first hip",
    "first knee",
    "p2 observed by administrative censoring",
    "revision before p2",
    "death before p2",
    "censored before p2 among no p2 by admin"
  ),
  value = c(
    nrow(hip_eligible),
    nrow(knee_eligible),
    nrow(cohort_p1_clean),
    dplyr::n_distinct(cohort_p1_clean$CODPAT),
    sum(cohort_p1_clean$has_eligible_hip & !cohort_p1_clean$has_eligible_knee, na.rm = TRUE),
    sum(!cohort_p1_clean$has_eligible_hip & cohort_p1_clean$has_eligible_knee, na.rm = TRUE),
    sum(cohort_p1_clean$has_eligible_hip & cohort_p1_clean$has_eligible_knee, na.rm = TRUE),
    sum(cohort_p1_clean$same_day_primaries, na.rm = TRUE),
    sum(cohort_p1_clean$hip == 1, na.rm = TRUE),
    sum(cohort_p1_clean$hip == 0, na.rm = TRUE),
    sum(cohort_p1_clean$has_p2_by_admin, na.rm = TRUE),
    sum(cohort_p1_clean$revision_before_p2, na.rm = TRUE),
    sum(cohort_p1_clean$death_before_p2, na.rm = TRUE),
    sum(cohort_p1_clean$censored_before_p2, na.rm = TRUE)
  )
)

qc_first_event <- cohort_p1_clean %>%
  count(first_event_from_p1, order_group, name = "n") %>%
  arrange(order_group, first_event_from_p1)

qc_missing <- cohort_p1_clean %>%
  summarise(
    n = n(),
    missing_sex = sum(is.na(sex)),
    missing_age_p1 = sum(is.na(age_p1)),
    missing_bmi_p1 = sum(is.na(bmi_p1)),
    missing_center_p1 = sum(is.na(center_p1_label) | center_p1_label == ""),
    missing_hip = sum(is.na(hip))
  )

qc_by_order <- cohort_p1_clean %>%
  group_by(order_group) %>%
  summarise(
    patient_sides = n(),
    patients = n_distinct(CODPAT),
    median_age_p1 = median(age_p1, na.rm = TRUE),
    median_bmi_p1 = median(bmi_p1, na.rm = TRUE),
    median_t12_raw = median(t12, na.rm = TRUE),
    p2_by_admin = sum(has_p2_by_admin, na.rm = TRUE),
    first_revision_after_p1 = sum(first_event_from_p1 == "revision", na.rm = TRUE),
    death_first_from_p1 = sum(first_event_from_p1 == "death", na.rm = TRUE),
    .groups = "drop"
  )

#salvataggio
save(
  cohort_p1_clean,
  event_history_clean,
  hip_eligible,
  knee_eligible,
  interventions_clean,
  death_clean,
  qc_overall,
  qc_first_event,
  qc_missing,
  qc_by_order,
  file = clean_data_file
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "cohort_p1_clean")
openxlsx::writeData(wb, "cohort_p1_clean", cohort_p1_clean)
openxlsx::addWorksheet(wb, "event_history_clean")
openxlsx::writeData(wb, "event_history_clean", event_history_clean)
openxlsx::addWorksheet(wb, "qc_overall")
openxlsx::writeData(wb, "qc_overall", qc_overall)
openxlsx::addWorksheet(wb, "qc_first_event")
openxlsx::writeData(wb, "qc_first_event", qc_first_event)
openxlsx::addWorksheet(wb, "qc_missing")
openxlsx::writeData(wb, "qc_missing", qc_missing)
openxlsx::addWorksheet(wb, "qc_by_order")
openxlsx::writeData(wb, "qc_by_order", qc_by_order)
for (sh in names(wb)) openxlsx::setColWidths(wb, sh, cols = 1:200, widths = "auto")
openxlsx::saveWorkbook(wb, file.path(out_dir, "p1_multistate_option3_clean_data.xlsx"), overwrite = TRUE)

message("output salvati in: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
print(qc_overall)
print(qc_by_order)

rm(list = setdiff(ls(), "cohort_p1_clean"))
gc()
ls()