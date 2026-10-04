source("config.R")

# ============================================================
# 09_p1_patient_bootstrap_B500.R
# bootstrap patient-level per probabilità origin-specific
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(survival)
  library(openxlsx)
})

#impostazioni
out_dir <- output_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

B <- 500
master_seed <- 20260829
checkpoint_every <- 25
resume_from_checkpoint <- TRUE

times_report <- c(1, 3, 5, 10, 15)
max_time <- max(times_report)
tol <- 1e-10

checkpoint_file <- file.path(out_dir, "p1_patient_bootstrap_B500_checkpoint.rds")

load(file.path(out_dir, "p1_multistate_analysis_objects.RData"))
load(file.path(out_dir, "p1_multistate_state_probabilities.RData"))

#specifica finale del file 05
expected_transitions <- c(
  "P1_to_P2",
  "P1_to_Rpre",
  "P1_to_D",
  "P2_to_Rpost",
  "P2_to_D",
  "Rpre_to_P2",
  "Rpre_to_D"
)

expected_tv <- c(
  "P1_to_P2",
  "P2_to_Rpost",
  "Rpre_to_P2",
  "Rpre_to_D"
)

if (!all(expected_transitions %in% transition_map$trans)) {
  stop("transition_map del file 05 non contiene tutte le transizioni attese")
}

spec_check <- transition_map %>%
  filter(trans %in% expected_transitions) %>%
  mutate(
    expected_tv = trans %in% expected_tv,
    specification_ok = hip_time_varying == expected_tv
  )

if (!all(spec_check$specification_ok)) {
  print(spec_check)
  stop("specifica hip(t) del file 05 non coerente con il bootstrap")
}

origin_order <- c("P1", "P2", "R_pre")

#funzioni
safe_cox <- function(formula, data, keep_model = TRUE, ...) {
  tryCatch(
    coxph(
      formula,
      data = data,
      ties = "efron",
      x = TRUE,
      y = TRUE,
      model = keep_model,
      ...
    ),
    error = function(e) e
  )
}

tt_log <- function(x, t, ...) {
  x * log(t + 1)
}

prepare_prob_data <- function(long_data, tr) {
  long_data %>%
    filter(as.character(trans) == tr) %>%
    mutate(
      hip = factor(
        as.character(hip),
        levels = c("first_knee", "first_hip")
      ),
      hip_binary = if_else(hip == "first_hip", 1, 0),
      age_cat = cut(
        age_p1,
        breaks = c(-Inf, 60, 70, 80, Inf),
        right = FALSE,
        labels = c("<60", "60-69", "70-79", "80+")
      ),
      sex_ms = factor(
        as.character(sex_ms),
        levels = c("male", "female")
      ),
      bmi_p1_f_ms = factor(
        as.character(bmi_p1_f_ms),
        levels = c(
          "normal_underweight",
          "overweight",
          "obese",
          "missing"
        )
      )
    ) %>%
    filter(
      !is.na(Tstop),
      Tstop > 0,
      !is.na(status),
      !is.na(hip),
      !is.na(age_cat),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms),
      !is.na(CODPAT)
    )
}

fit_probability_model <- function(long_data, tr) {
  meta <- transition_map %>%
    filter(.data$trans == .env$tr)
  
  d <- prepare_prob_data(long_data, tr)
  hip_tv <- meta$hip_time_varying[1]
  
  if (nrow(d) == 0 || sum(d$status == 1, na.rm = TRUE) == 0) {
    return(list(
      trans = tr,
      origin = meta$origin[1],
      destination = meta$destination[1],
      hip_tv = hip_tv,
      data = d,
      fit = structure(
        list(message = "nessun evento nel campione bootstrap"),
        class = "error"
      )
    ))
  }
  
  if (hip_tv) {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip_binary +
        tt(hip_binary) +
        sex_ms +
        bmi_p1_f_ms +
        strata(age_cat) +
        cluster(CODPAT),
      data = d,
      keep_model = FALSE,
      tt = tt_log
    )
  } else {
    fit <- safe_cox(
      Surv(Tstop, status) ~
        hip +
        sex_ms +
        bmi_p1_f_ms +
        strata(age_cat) +
        cluster(CODPAT),
      data = d
    )
  }
  
  list(
    trans = tr,
    origin = meta$origin[1],
    destination = meta$destination[1],
    hip_tv = hip_tv,
    data = d,
    fit = fit
  )
}

calc_lp_times <- function(pattern, beta, times) {
  if (!is.data.frame(pattern)) {
    pattern <- as.data.frame(pattern)
  }
  
  n_pattern <- nrow(pattern)
  n_times <- length(times)
  n <- max(n_pattern, n_times)
  
  if (n_pattern == 1 && n > 1) {
    pattern <- pattern[rep(1, n), , drop = FALSE]
  } else if (n_pattern != n) {
    stop("dimensioni incompatibili tra pattern e times")
  }
  
  if (n_times == 1 && n > 1) {
    times <- rep(times, n)
  } else if (n_times != n) {
    stop("dimensioni incompatibili tra pattern e times")
  }
  
  hip_value <- as.numeric(pattern$hip_binary)
  sex_value <- as.character(pattern$sex_ms)
  bmi_value <- as.character(pattern$bmi_p1_f_ms)
  
  lp <- rep(0, n)
  
  if ("hipfirst_hip" %in% names(beta)) {
    lp <- lp + beta["hipfirst_hip"] * hip_value
  }
  
  if ("hip_binary" %in% names(beta)) {
    lp <- lp + beta["hip_binary"] * hip_value
  }
  
  if ("tt(hip_binary)" %in% names(beta)) {
    lp <- lp +
      beta["tt(hip_binary)"] *
      hip_value *
      log(times + 1)
  }
  
  if ("sex_msfemale" %in% names(beta)) {
    lp <- lp +
      beta["sex_msfemale"] *
      as.numeric(sex_value == "female")
  }
  
  if ("bmi_p1_f_msoverweight" %in% names(beta)) {
    lp <- lp +
      beta["bmi_p1_f_msoverweight"] *
      as.numeric(bmi_value == "overweight")
  }
  
  if ("bmi_p1_f_msobese" %in% names(beta)) {
    lp <- lp +
      beta["bmi_p1_f_msobese"] *
      as.numeric(bmi_value == "obese")
  }
  
  if ("bmi_p1_f_msmissing" %in% names(beta)) {
    lp <- lp +
      beta["bmi_p1_f_msmissing"] *
      as.numeric(bmi_value == "missing")
  }
  
  as.numeric(lp)
}

make_efron_baseline <- function(model_obj) {
  tr <- model_obj$trans
  d <- model_obj$data
  fit <- model_obj$fit
  
  if (inherits(fit, "error")) {
    stop(paste0(tr, ": ", fit$message))
  }
  
  beta <- coef(fit)
  
  if (length(beta) == 0 || any(!is.finite(beta))) {
    stop(paste0(tr, ": coefficienti non finiti"))
  }
  
  age_levels <- levels(d$age_cat)
  
  map_dfr(
    age_levels,
    function(age_level) {
      d_s <- d %>%
        filter(as.character(age_cat) == age_level)
      
      event_times <- sort(
        unique(
          d_s$Tstop[
            d_s$status == 1 &
              d_s$Tstop <= max_time
          ]
        )
      )
      
      if (length(event_times) == 0) {
        return(
          tibble(
            trans = character(),
            age_stratum = character(),
            time = numeric(),
            dLambda0 = numeric()
          )
        )
      }
      
      map_dfr(
        event_times,
        function(tt) {
          risk_idx <- d_s$Tstop >= tt - tol
          event_idx <-
            d_s$status == 1 &
            abs(d_s$Tstop - tt) < tol
          
          n_events <- sum(event_idx, na.rm = TRUE)
          
          risk_pattern <- d_s[risk_idx, , drop = FALSE]
          event_pattern <- d_s[event_idx, , drop = FALSE]
          
          lp_risk <- calc_lp_times(
            risk_pattern,
            beta,
            rep(tt, nrow(risk_pattern))
          )
          
          lp_event <- calc_lp_times(
            event_pattern,
            beta,
            rep(tt, nrow(event_pattern))
          )
          
          risk_weight <- sum(exp(lp_risk))
          event_weight <- sum(exp(lp_event))
          
          den <-
            risk_weight -
            (seq(0, n_events - 1) / n_events) *
            event_weight
          
          if (any(!is.finite(den)) || any(den <= 0)) {
            stop(
              paste0(
                tr,
                ": denominatore Efron non valido; age=",
                age_level,
                ", t=",
                signif(tt, 8)
              )
            )
          }
          
          tibble(
            trans = tr,
            age_stratum = age_level,
            time = tt,
            dLambda0 = sum(1 / den)
          )
        }
      )
    }
  )
}

make_pattern <- function(age_stratum, sex_value, bmi_value, hip_value) {
  tibble(
    age_cat = factor(
      age_stratum,
      levels = c("<60", "60-69", "70-79", "80+")
    ),
    sex_ms = factor(
      sex_value,
      levels = c("male", "female")
    ),
    bmi_p1_f_ms = factor(
      bmi_value,
      levels = c(
        "normal_underweight",
        "overweight",
        "obese",
        "missing"
      )
    ),
    hip = factor(
      if_else(hip_value == 1, "first_hip", "first_knee"),
      levels = c("first_knee", "first_hip")
    ),
    hip_binary = as.numeric(hip_value)
  )
}

get_origin_subjects <- function(long_data, origin_name) {
  long_data %>%
    filter(as.character(from) == origin_name) %>%
    mutate(
      hip = factor(
        as.character(hip),
        levels = c("first_knee", "first_hip")
      ),
      age_cat = cut(
        age_p1,
        breaks = c(-Inf, 60, 70, 80, Inf),
        right = FALSE,
        labels = c("<60", "60-69", "70-79", "80+")
      ),
      sex_ms = factor(
        as.character(sex_ms),
        levels = c("male", "female")
      ),
      bmi_p1_f_ms = factor(
        as.character(bmi_p1_f_ms),
        levels = c(
          "normal_underweight",
          "overweight",
          "obese",
          "missing"
        )
      )
    ) %>%
    distinct(patient_side_id, .keep_all = TRUE) %>%
    filter(
      !is.na(hip),
      !is.na(age_cat),
      !is.na(sex_ms),
      !is.na(bmi_p1_f_ms)
    )
}

make_pattern_weights <- function(long_data, origin_name) {
  d <- get_origin_subjects(long_data, origin_name)
  n_total <- nrow(d)
  
  if (n_total == 0) {
    stop(paste0(origin_name, ": popolazione di standardizzazione vuota"))
  }
  
  d %>%
    count(
      age_cat,
      sex_ms,
      bmi_p1_f_ms,
      name = "n_pattern",
      .drop = FALSE
    ) %>%
    filter(n_pattern > 0) %>%
    mutate(
      weight = n_pattern / n_total,
      n_standardization = n_total
    )
}

predict_pattern <- function(
    origin_name,
    age_stratum,
    sex_value,
    bmi_value,
    hip_value,
    model_objects,
    baseline_increments,
    times = times_report
) {
  meta_origin <- transition_map %>%
    filter(origin == origin_name)
  
  dest_states <- unique(meta_origin$destination)
  
  pattern <- make_pattern(
    age_stratum,
    sex_value,
    bmi_value,
    hip_value
  )
  
  hazard_parts <- map_dfr(
    meta_origin$trans,
    function(tr) {
      inc <- baseline_increments %>%
        filter(
          .data$trans == .env$tr,
          .data$age_stratum == .env$age_stratum,
          !is.na(.data$time),
          !is.na(.data$dLambda0),
          .data$time <= max(times)
        )
      
      if (nrow(inc) == 0) {
        return(
          tibble(
            time = numeric(),
            destination = character(),
            dh = numeric()
          )
        )
      }
      
      obj <- model_objects[[tr]]
      beta <- coef(obj$fit)
      
      lp <- calc_lp_times(
        pattern,
        beta,
        inc$time
      )
      
      destination <- meta_origin %>%
        filter(.data$trans == .env$tr) %>%
        pull(destination)
      
      tibble(
        time = inc$time,
        destination = destination,
        dh = inc$dLambda0 * exp(lp)
      )
    }
  )
  
  hazard_by_time <- hazard_parts %>%
    group_by(time, destination) %>%
    summarise(dh = sum(dh), .groups = "drop") %>%
    arrange(time)
  
  all_event_times <- sort(unique(hazard_by_time$time))
  
  survival_origin <- 1
  cif <- setNames(rep(0, length(dest_states)), dest_states)
  event_index <- 1
  out <- vector("list", length(times))
  
  for (ii in seq_along(times)) {
    target <- times[ii]
    
    while (
      event_index <= length(all_event_times) &&
      all_event_times[event_index] <= target + tol
    ) {
      tt <- all_event_times[event_index]
      
      d_tt <- hazard_by_time %>%
        filter(abs(time - tt) < tol)
      
      d_total <- sum(d_tt$dh)
      
      if (!is.finite(d_total) || d_total < -tol) {
        stop(
          paste0(
            "incremento di hazard non valido; origin=",
            origin_name,
            ", age=",
            age_stratum,
            ", t=",
            signif(tt, 8)
          )
        )
      }
      
      survival_before <- survival_origin
      
      if (d_total > tol) {
        event_probability <- 1 - exp(-d_total)
        
        for (jj in seq_len(nrow(d_tt))) {
          dest <- d_tt$destination[jj]
          hazard_share <- d_tt$dh[jj] / d_total
          
          cif[dest] <-
            cif[dest] +
            survival_before *
            event_probability *
            hazard_share
        }
        
        survival_origin <-
          survival_before *
          exp(-d_total)
      }
      
      event_index <- event_index + 1
    }
    
    out[[ii]] <- bind_rows(
      tibble(
        origin = origin_name,
        time_years = target,
        state = origin_name,
        probability = survival_origin
      ),
      tibble(
        origin = origin_name,
        time_years = target,
        state = names(cif),
        probability = as.numeric(cif)
      )
    )
  }
  
  bind_rows(out)
}

predict_standardized_origin <- function(
    long_data,
    origin_name,
    hip_value,
    model_objects,
    baseline_increments,
    times = times_report
) {
  patterns <- make_pattern_weights(long_data, origin_name)
  
  group_name <- if_else(
    hip_value == 1,
    "first_hip",
    "first_knee"
  )
  
  pred_patterns <- map_dfr(
    seq_len(nrow(patterns)),
    function(i) {
      pat <- patterns[i, , drop = FALSE]
      
      predict_pattern(
        origin_name = origin_name,
        age_stratum = as.character(pat$age_cat),
        sex_value = as.character(pat$sex_ms),
        bmi_value = as.character(pat$bmi_p1_f_ms),
        hip_value = hip_value,
        model_objects = model_objects,
        baseline_increments = baseline_increments,
        times = times
      ) %>%
        mutate(pattern_weight = pat$weight)
    }
  )
  
  pred_patterns %>%
    group_by(origin, time_years, state) %>%
    summarise(
      probability = sum(probability * pattern_weight),
      .groups = "drop"
    ) %>%
    mutate(group = group_name)
}

make_boot_sample <- function(long_data, seed) {
  set.seed(seed)
  
  source_long <- long_data %>%
    mutate(
      CODPAT_original = as.character(CODPAT),
      patient_side_original = as.character(patient_side_id)
    )
  
  patient_ids <- unique(source_long$CODPAT_original)
  n_patients <- length(patient_ids)
  
  draws <- tibble(
    draw_id = seq_len(n_patients),
    CODPAT_original = sample(
      patient_ids,
      size = n_patients,
      replace = TRUE
    )
  )
  
  suppressWarnings(
    draws %>%
      left_join(
        source_long,
        by = "CODPAT_original"
      )
  ) %>%
    mutate(
      CODPAT = paste0("boot_patient_", draw_id),
      patient_side_id = paste0(
        patient_side_original,
        "__boot_",
        draw_id
      )
    ) %>%
    select(
      -CODPAT_original,
      -patient_side_original
    )
}

run_origin_bootstrap <- function(boot_long, origin_name) {
  trs <- transition_map %>%
    filter(origin == origin_name) %>%
    pull(trans)
  
  model_objects <- map(
    trs,
    function(tr) fit_probability_model(boot_long, tr)
  )
  names(model_objects) <- trs
  
  fit_errors <- map_chr(
    model_objects,
    function(x) {
      if (inherits(x$fit, "error")) {
        x$fit$message
      } else {
        ""
      }
    }
  )
  
  if (any(nzchar(fit_errors))) {
    bad <- names(fit_errors)[nzchar(fit_errors)]
    stop(
      paste0(
        "fit fallito: ",
        paste(
          paste0(bad, " [", fit_errors[bad], "]"),
          collapse = "; "
        )
      )
    )
  }
  
  coef_ok <- map_lgl(
    model_objects,
    function(x) {
      b <- coef(x$fit)
      length(b) > 0 && all(is.finite(b))
    }
  )
  
  if (!all(coef_ok)) {
    stop(
      paste0(
        "coefficienti non finiti in: ",
        paste(names(coef_ok)[!coef_ok], collapse = ", ")
      )
    )
  }
  
  baseline_increments <- map_dfr(
    model_objects,
    make_efron_baseline
  )
  
  pred_all <- bind_rows(
    predict_standardized_origin(
      boot_long,
      origin_name,
      0,
      model_objects,
      baseline_increments,
      times_report
    ),
    predict_standardized_origin(
      boot_long,
      origin_name,
      1,
      model_objects,
      baseline_increments,
      times_report
    )
  )
  
  sum_check <- pred_all %>%
    group_by(group, time_years) %>%
    summarise(
      total_probability = sum(probability),
      min_probability = min(probability),
      max_probability = max(probability),
      .groups = "drop"
    ) %>%
    mutate(
      valid =
        abs(total_probability - 1) < 1e-8 &
        min_probability >= -1e-10 &
        max_probability <= 1 + 1e-10
    )
  
  if (!all(sum_check$valid)) {
    stop("controllo delle probabilità non superato")
  }
  
  probabilities <- pred_all %>%
    filter(as.character(state) != as.character(origin)) %>%
    select(
      origin,
      time_years,
      state,
      group,
      probability
    )
  
  contrasts <- probabilities %>%
    pivot_wider(
      names_from = group,
      values_from = probability
    ) %>%
    mutate(
      risk_difference = first_hip - first_knee,
      risk_ratio = if_else(
        first_knee > 0,
        first_hip / first_knee,
        NA_real_
      )
    )
  
  list(
    probabilities = probabilities,
    contrasts = contrasts
  )
}

run_bootstrap_rep <- function(b, seed) {
  warning_messages <- character()
  
  out <- withCallingHandlers(
    {
      boot_long <- make_boot_sample(p1_ms_long, seed)
      
      prob_list <- list()
      contrast_list <- list()
      error_list <- list()
      
      for (origin_name in origin_order) {
        origin_out <- tryCatch(
          run_origin_bootstrap(
            boot_long,
            origin_name
          ),
          error = function(e) e
        )
        
        if (inherits(origin_out, "error")) {
          error_list[[origin_name]] <- tibble(
            replicate = b,
            origin = origin_name,
            error = origin_out$message
          )
        } else {
          prob_list[[origin_name]] <- origin_out$probabilities %>%
            mutate(replicate = b, .before = 1)
          
          contrast_list[[origin_name]] <- origin_out$contrasts %>%
            mutate(replicate = b, .before = 1)
        }
      }
      
      list(
        probabilities = bind_rows(prob_list),
        contrasts = bind_rows(contrast_list),
        errors = bind_rows(error_list)
      )
    },
    warning = function(w) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )
  
  out$warnings <- if (length(warning_messages) == 0) {
    tibble(
      replicate = integer(),
      warning = character()
    )
  } else {
    tibble(
      replicate = b,
      warning = unique(warning_messages)
    )
  }
  
  out
}

#semi riproducibili per singola replica
set.seed(master_seed)
bootstrap_seeds <- sample.int(
  .Machine$integer.max,
  B,
  replace = FALSE
)

probability_results <- vector("list", B)
contrast_results <- vector("list", B)
error_results <- vector("list", B)
warning_results <- vector("list", B)
done <- rep(FALSE, B)
rep_elapsed_seconds <- rep(NA_real_, B)

#resume opzionale
if (resume_from_checkpoint && file.exists(checkpoint_file)) {
  ck <- readRDS(checkpoint_file)
  
  if (
    identical(ck$B, B) &&
    identical(ck$master_seed, master_seed) &&
    identical(ck$bootstrap_seeds, bootstrap_seeds)
  ) {
    probability_results <- ck$probability_results
    contrast_results <- ck$contrast_results
    error_results <- ck$error_results
    warning_results <- ck$warning_results
    done <- ck$done
    rep_elapsed_seconds <- ck$rep_elapsed_seconds
    message("checkpoint caricato: ", sum(done), "/", B, " repliche già completate")
  } else {
    stop("checkpoint non compatibile con B/seed correnti")
  }
}

cat("\nbootstrap patient-level: B = ", B, "\n", sep = "")
cat("la barra mostra la percentuale di repliche completate\n\n")

pb <- txtProgressBar(
  min = 0,
  max = B,
  initial = sum(done),
  style = 3
)

start_time <- Sys.time()

for (b in which(!done)) {
  rep_start <- Sys.time()
  
  out_b <- tryCatch(
    run_bootstrap_rep(
      b = b,
      seed = bootstrap_seeds[b]
    ),
    error = function(e) {
      list(
        probabilities = tibble(),
        contrasts = tibble(),
        errors = tibble(
          replicate = b,
          origin = "all",
          error = e$message
        ),
        warnings = tibble()
      )
    }
  )
  
  probability_results[[b]] <- out_b$probabilities
  contrast_results[[b]] <- out_b$contrasts
  error_results[[b]] <- out_b$errors
  warning_results[[b]] <- out_b$warnings
  
  rep_elapsed_seconds[b] <- as.numeric(
    difftime(
      Sys.time(),
      rep_start,
      units = "secs"
    )
  )
  
  done[b] <- TRUE
  setTxtProgressBar(pb, sum(done))
  
  if (
    b %% checkpoint_every == 0 ||
    all(done)
  ) {
    saveRDS(
      list(
        B = B,
        master_seed = master_seed,
        bootstrap_seeds = bootstrap_seeds,
        probability_results = probability_results,
        contrast_results = contrast_results,
        error_results = error_results,
        warning_results = warning_results,
        done = done,
        rep_elapsed_seconds = rep_elapsed_seconds
      ),
      checkpoint_file
    )
  }
}

close(pb)

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

cat(
  "\nbootstrap completato in ",
  round(elapsed_minutes, 1),
  " minuti\n",
  sep = ""
)

bootstrap_probabilities <- bind_rows(probability_results)
bootstrap_contrasts <- bind_rows(contrast_results)
bootstrap_errors <- bind_rows(error_results)
bootstrap_warnings <- bind_rows(warning_results)

#point estimates dal file 05
point_probabilities <- adjusted_probabilities_report %>%
  filter(
    time_years %in% times_report,
    as.character(state) != as.character(origin)
  ) %>%
  transmute(
    origin = as.character(origin),