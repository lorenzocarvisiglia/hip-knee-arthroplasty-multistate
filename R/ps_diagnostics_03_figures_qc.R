# ============================================================
# figure
# ============================================================

ps_density_path <- file.path(fig_dir, "p1_ps_density_first_hip_vs_first_knee.png")
ps_hist_path <- file.path(fig_dir, "p1_ps_histogram_first_hip_vs_first_knee.png")
weight_density_path <- file.path(fig_dir, "p1_overlap_weight_density.png")
weight_hist_path <- file.path(fig_dir, "p1_overlap_weight_histogram.png")
love_plot_path <- file.path(fig_dir, "p1_love_plot_overlap_weights.png")
smd_summary_path <- file.path(fig_dir, "p1_smd_summary_before_after_weighting.png")

p_ps_density <- ggplot(
  ps_data,
  aes(
    x = ps_first_hip,
    linetype = order_group_label
  )
) +
  geom_density(linewidth = 0.9) +
  labs(
    x = "Propensity score for first hip",
    y = "Density",
    linetype = "",
    title = "Propensity-score overlap by implant-order group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = ps_density_path,
  plot = p_ps_density,
  width = 7,
  height = 5,
  dpi = 300
)

p_ps_hist <- ggplot(
  ps_data,
  aes(
    x = ps_first_hip,
    fill = order_group_label
  )
) +
  geom_histogram(
    bins = 40,
    position = "identity",
    alpha = 0.45
  ) +
  labs(
    x = "Propensity score for first hip",
    y = "Number of patient-side units",
    fill = "",
    title = "Propensity-score distribution by implant-order group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = ps_hist_path,
  plot = p_ps_hist,
  width = 7,
  height = 5,
  dpi = 300
)

p_weight_density <- ggplot(
  ps_data,
  aes(
    x = overlap_weight,
    linetype = order_group_label
  )
) +
  geom_density(linewidth = 0.9) +
  labs(
    x = "Overlap weight",
    y = "Density",
    linetype = "",
    title = "Overlap-weight distribution by implant-order group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = weight_density_path,
  plot = p_weight_density,
  width = 7,
  height = 5,
  dpi = 300
)

p_weight_hist <- ggplot(
  ps_data,
  aes(
    x = overlap_weight,
    fill = order_group_label
  )
) +
  geom_histogram(
    bins = 40,
    position = "identity",
    alpha = 0.45
  ) +
  labs(
    x = "Overlap weight",
    y = "Number of patient-side units",
    fill = "",
    title = "Overlap-weight distribution by implant-order group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = weight_hist_path,
  plot = p_weight_hist,
  width = 7,
  height = 5,
  dpi = 300
)

p_love <- ggplot(
  love_plot_data,
  aes(
    x = abs_smd,
    y = variable_level,
    shape = sample
  )
) +
  geom_vline(xintercept = 0.10, linetype = "dashed") +
  geom_point(size = 2.2) +
  labs(
    x = "Absolute standardized mean difference",
    y = "",
    shape = "",
    title = "Covariate balance before and after overlap weighting"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = love_plot_path,
  plot = p_love,
  width = 8,
  height = max(5, 0.18 * n_distinct(love_plot_data$variable_level) + 2),
  dpi = 300
)

p_smd_summary <- ggplot(
  balance_overall_summary,
  aes(
    x = sample,
    y = max_abs_smd
  )
) +
  geom_col() +
  geom_hline(yintercept = 0.10, linetype = "dashed") +
  labs(
    x = "",
    y = "Maximum absolute SMD",
    title = "Maximum imbalance before and after overlap weighting"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = smd_summary_path,
  plot = p_smd_summary,
  width = 6,
  height = 4,
  dpi = 300
)

figures_manifest <- tibble::tibble(
  figure = c(
    "Propensity-score density",
    "Propensity-score histogram",
    "Overlap-weight density",
    "Overlap-weight histogram",
    "Love plot",
    "SMD summary"
  ),
  path = c(
    ps_density_path,
    ps_hist_path,
    weight_density_path,
    weight_hist_path,
    love_plot_path,
    smd_summary_path
  )
)

# ============================================================
# QC e testo metodi
# ============================================================

n_first_hip <- sum(ps_data$order_group == "first_hip", na.rm = TRUE)
n_first_knee <- sum(ps_data$order_group == "first_knee", na.rm = TRUE)

max_weight <- max(ps_data$overlap_weight, na.rm = TRUE)
p99_weight <- as.numeric(quantile(ps_data$overlap_weight, 0.99, na.rm = TRUE))

max_smd_unweighted <- balance_overall_summary %>%
  filter(sample == "Unweighted") %>%
  pull(max_abs_smd)

max_smd_weighted <- balance_overall_summary %>%
  filter(sample == "Overlap weighted") %>%
  pull(max_abs_smd)

n_smd_gt_010_weighted <- balance_overall_summary %>%
  filter(sample == "Overlap weighted") %>%
  pull(n_terms_abs_smd_gt_0_10)

qc <- tibble::tibble(
  check = c(
    "source RData path",
    "source object",
    "n patient-side units in PS dataset",
    "n first hip",
    "n first knee",
    "expected first hip count",
    "expected first knee count",
    "PS model formula",
    "glm convergence",
    "glm warnings",
    "minimum PS",
    "maximum PS",
    "maximum overlap weight",
    "p99 overlap weight",
    "overall ESS overlap",
    "max abs SMD unweighted",
    "max abs SMD overlap weighted",
    "weighted SMD terms > 0.10",
    "figures created"
  ),
  value = c(
    source_path,
    source_object,
    as.character(nrow(ps_data)),
    as.character(n_first_hip),
    as.character(n_first_knee),
    "4319",
    "5384",
    paste(deparse(ps_formula), collapse = " "),
    as.character(ps_fit$converged),
    ifelse(length(glm_warnings) == 0, "none", paste(unique(glm_warnings), collapse = " | ")),
    as.character(min(ps_data$ps_first_hip, na.rm = TRUE)),
    as.character(max(ps_data$ps_first_hip, na.rm = TRUE)),
    as.character(max_weight),
    as.character(p99_weight),
    as.character(ess_summary$ess_overlap[ess_summary$group == "Overall"]),
    as.character(max_smd_unweighted),
    as.character(max_smd_weighted),
    as.character(n_smd_gt_010_weighted),
    as.character(all(file.exists(figures_manifest$path)))
  ),
  status = case_when(
    check %in% c("source RData path", "source object", "PS model formula") ~ "OK",
    check == "n patient-side units in PS dataset" & value == "9703" ~ "OK",
    check == "n first hip" & value == "4319" ~ "OK",
    check == "n first knee" & value == "5384" ~ "OK",
    check == "expected first hip count" ~ "INFO",
    check == "expected first knee count" ~ "INFO",
    check == "glm convergence" & value == "TRUE" ~ "OK",
    check == "glm warnings" & value == "none" ~ "OK",
    check == "minimum PS" & suppressWarnings(as.numeric(value)) > 0 ~ "OK",
    check == "maximum PS" & suppressWarnings(as.numeric(value)) < 1 ~ "OK",
    check == "maximum overlap weight" & suppressWarnings(as.numeric(value)) <= 1 ~ "OK",
    check == "p99 overlap weight" & suppressWarnings(as.numeric(value)) <= 1 ~ "OK",
    check == "overall ESS overlap" & suppressWarnings(as.numeric(value)) > 0 ~ "OK",
    check == "max abs SMD overlap weighted" & suppressWarnings(as.numeric(value)) <= 0.10 ~ "OK",
    check == "weighted SMD terms > 0.10" & value == "0" ~ "OK",
    check == "figures created" & value == "TRUE" ~ "OK",
    TRUE ~ "CHECK"
  ),
  note = case_when(
    check == "n patient-side units in PS dataset" & value != "9703" ~
      "Expected 9703 defined-order patient-side units from the P1 cohort.",
    check == "n first hip" & value != "4319" ~
      "Expected 4319 first-hip units from the P1 cohort.",
    check == "n first knee" & value != "5384" ~
      "Expected 5384 first-knee units from the P1 cohort.",
    check == "glm warnings" & value != "none" ~
      "Inspect warnings; possible separation or sparse covariate levels.",
    check == "max abs SMD overlap weighted" & suppressWarnings(as.numeric(value)) > 0.10 ~
      "Some residual imbalance remains after overlap weighting.",
    check == "weighted SMD terms > 0.10" & value != "0" ~
      "Inspect balance_long and love plot.",
    TRUE ~ ""
  )
)

methods_text <- tibble::tibble(
  section = c(
    "Propensity score model",
    "Overlap weights",
    "Overlap population",
    "Balance diagnostics",
    "Effective sample size",
    "Interpretation"
  ),
  text = c(
    paste0(
      "A propensity score for receiving hip rather than knee as the first ipsilateral primary arthroplasty was estimated at P1 using logistic regression. The model included available baseline variables: ",
      paste(covariates_used$manuscript_label[covariates_used$used_in_ps_model], collapse = ", "),
      "."
    ),
    "Overlap weights were defined as 1 - e(X) for first-hip patient-side units and e(X) for first-knee patient-side units, where e(X) is the estimated propensity score for first hip.",
    "The overlap-weighted analysis targets the population of patient-side units with clinical and covariate profiles compatible with either implant-order pathway, conditional on measured baseline covariates.",
    "Covariate balance before and after weighting was assessed using standardized mean differences for continuous covariates and binary indicators for categorical covariate levels.",
    "The effective sample size was computed as the squared sum of overlap weights divided by the sum of squared overlap weights.",
    "The weighted analysis should be interpreted as a measured-confounder sensitivity analysis and not as proof of a causal effect of surgical order."
  )
)

# ============================================================
# salva RData
# ============================================================

save(
  source_data,
  source_path,
  source_object,
  column_mapping,
  ps_data,
  covariates_used,
  ps_formula,
  ps_fit,
  glm_warnings,
  ps_model_coef,
  ps_summary_by_group,
  weight_summary_by_group,
  ess_summary,
  support_summary,
  common_support,
  balance_long,
  balance_wide,
  balance_summary_by_covariate,
  balance_overall_summary,
  love_plot_data,
  figures_manifest,
  qc,
  methods_text,
  file = out_rdata
)
