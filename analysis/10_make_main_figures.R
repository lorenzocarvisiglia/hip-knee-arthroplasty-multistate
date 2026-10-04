source("config.R")

#final black-and-white figures for adjusted probabilities and bootstrap uncertainty

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(openxlsx)
})

out_dir <- output_dir
fig_dir <- file.path(out_dir, "figures_final")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

prob_file <- file.path(out_dir, "p1_multistate_state_probabilities.RData")
boot_file <- file.path(out_dir, "p1_patient_bootstrap_probabilities_B500.RData")
ow_file <- file.path(out_dir, "p1_overlap_weighted_P1_probabilities_B500.RData")

for (f in c(prob_file, boot_file, ow_file)) {
  if (!file.exists(f)) stop("required file not found: ", f)
}

load_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}

pe <- load_env(prob_file)
be <- load_env(boot_file)
oe <- load_env(ow_file)

required_prob <- c("adjusted_probabilities_curve", "adjusted_probabilities_report")
required_boot <- c("probability_ci", "contrast_ci", "origin_success")
required_ow <- c("bootstrap_probability_ci", "bootstrap_contrast_ci", "success_by_method", "qc")

if (!all(vapply(required_prob, exists, logical(1), envir = pe))) stop("05 RData missing expected objects")
if (!all(vapply(required_boot, exists, logical(1), envir = be))) stop("09 RData missing expected objects")
if (!all(vapply(required_ow, exists, logical(1), envir = oe))) stop("14 RData missing expected objects")

adjusted_curve <- as_tibble(pe$adjusted_probabilities_curve)
adjusted_report <- as_tibble(pe$adjusted_probabilities_report)
probability_ci <- as_tibble(be$probability_ci)
contrast_ci <- as_tibble(be$contrast_ci)
origin_success <- as_tibble(be$origin_success)
ow_probability_ci <- as_tibble(oe$bootstrap_probability_ci)
ow_contrast_ci <- as_tibble(oe$bootstrap_contrast_ci)
ow_success <- as_tibble(oe$success_by_method)
ow_qc <- as_tibble(oe$qc)

if (!all(origin_success$successful_replicates == 500)) stop("main bootstrap does not have 500 successful replicates for every origin")
if (!all(ow_qc$status == "PASS")) stop("overlap-weighted bootstrap QC is not fully PASS")
if (!any(ow_success$method == "overlap_weighted_AJ" & ow_success$successful_replicates == 500)) stop("overlap-weighted bootstrap does not have 500 successful replicates")

transition_from_origin_state <- function(origin, state) {
  case_when(
    origin == "P1" & state == "P2" ~ "P1_to_P2",
    origin == "P1" & state == "R_pre" ~ "P1_to_Rpre",
    origin == "P1" & state == "D" ~ "P1_to_D",
    origin == "P2" & state == "R_post" ~ "P2_to_Rpost",
    origin == "P2" & state == "D" ~ "P2_to_D",
    origin == "R_pre" & state == "P2_after_Rpre" ~ "Rpre_to_P2",
    origin == "R_pre" & state == "D" ~ "Rpre_to_D",
    TRUE ~ NA_character_
  )
}

transition_label <- function(x) {
  recode(
    x,
    P1_to_P2 = "P1 to P2",
    P1_to_Rpre = "P1 to pre-P2 revision",
    P1_to_D = "P1 to death",
    P2_to_Rpost = "Direct P2 to post-P2 revision",
    P2_to_D = "Direct P2 to death",
    Rpre_to_P2 = "Pre-P2 revision to P2",
    Rpre_to_D = "Pre-P2 revision to death",
    .default = x
  )
}

group_label <- function(x) {
  recode(x, first_knee = "First knee", first_hip = "First hip", .default = x)
}

transition_levels <- c(
  "P1_to_P2", "P1_to_Rpre", "P1_to_D",
  "P2_to_Rpost", "P2_to_D",
  "Rpre_to_P2", "Rpre_to_D"
)

#validate that bootstrap point estimates reproduce the frozen 05 report
point_check <- adjusted_report %>%
  mutate(
    origin = as.character(origin),
    state = as.character(state),
    group = as.character(group)
  ) %>%
  filter(time_years > 0, state != origin) %>%
  select(origin, state, group, time_years, point_05 = probability) %>%
  inner_join(
    probability_ci %>%
      transmute(
        origin = as.character(origin),
        state = as.character(state),
        group = as.character(group),
        time_years,
        point_09 = point_estimate
      ),
    by = c("origin", "state", "group", "time_years")
  ) %>%
  mutate(abs_diff = abs(point_05 - point_09))

if (nrow(point_check) == 0 || max(point_check$abs_diff, na.rm = TRUE) > 1e-10) {
  stop("05 and 09 point estimates are not aligned")
}

curve_data <- adjusted_curve %>%
  mutate(
    origin = as.character(origin),
    state = as.character(state),
    group = as.character(group),
    transition = transition_from_origin_state(origin, state),
    transition_label = transition_label(transition),
    group_label = group_label(group)
  ) %>%
  filter(!is.na(transition)) %>%
  mutate(
    transition = factor(transition, levels = transition_levels),
    group_label = factor(group_label, levels = c("First knee", "First hip"))
  )

ci_data <- probability_ci %>%
  mutate(
    origin = as.character(origin),
    state = as.character(state),
    group = as.character(group),
    transition = transition_from_origin_state(origin, state),
    transition_label = transition_label(transition),
    group_label = group_label(group)
  ) %>%
  filter(!is.na(transition)) %>%
  mutate(
    transition = factor(transition, levels = transition_levels),
    group_label = factor(group_label, levels = c("First knee", "First hip"))
  )

rd_data <- contrast_ci %>%
  mutate(
    origin = as.character(origin),
    state = as.character(state),
    transition = transition_from_origin_state(origin, state),
    transition_label = transition_label(transition),
    transition = factor(transition, levels = transition_levels)
  ) %>%
  filter(!is.na(transition))

#black-and-white publication settings
line_values <- c("First knee" = "solid", "First hip" = "longdash")
shape_values <- c("First knee" = 16, "First hip" = 17)
ci_grey_values <- c("First knee" = "grey65", "First hip" = "grey35")

common_prob_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )

save_plot <- function(plot, stub, width, height) {
  ggsave(file.path(fig_dir, paste0(stub, ".png")), plot, width = width, height = height, dpi = 400)
  ggsave(file.path(fig_dir, paste0(stub, ".pdf")), plot, width = width, height = height)
}

make_probability_plot <- function(transitions, x_label, width, height, filename_stub, ncol = NULL) {
  cd <- curve_data %>% filter(as.character(transition) %in% transitions)
  id <- ci_data %>% filter(as.character(transition) %in% transitions)
  
  p <- ggplot(cd, aes(x = time_years, y = probability, linetype = group_label, group = group_label)) +
    geom_line(linewidth = 0.9, colour = "black") +
    geom_errorbar(
      data = id,
      aes(
        x = time_years,
        ymin = lower_95,
        ymax = upper_95,
        colour = group_label,
        group = group_label
      ),
      width = 0.18,
      linewidth = 0.6,
      position = position_dodge(width = 0.22),
      inherit.aes = FALSE
    ) +
    geom_point(
      data = id,
      aes(
        x = time_years,
        y = point_estimate,
        shape = group_label,
        group = group_label
      ),
      size = 2.2,
      colour = "black",
      position = position_dodge(width = 0.22),
      inherit.aes = FALSE
    ) +
    facet_wrap(~ transition_label, scales = "free_y", ncol = ncol) +
    scale_linetype_manual(values = line_values) +
    scale_shape_manual(values = shape_values) +
    scale_colour_manual(values = ci_grey_values, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.08))) +
    scale_x_continuous(breaks = c(0, 1, 3, 5, 10, 15), limits = c(0, 15)) +
    labs(
      x = x_label,
      y = "Adjusted cumulative incidence"
    ) +
    common_prob_theme
  
  save_plot(p, filename_stub, width, height)
  p
}

#main P1 figure: exact filename already used in Overleaf
p_p1 <- make_probability_plot(
  c("P1_to_P2", "P1_to_Rpre", "P1_to_D"),
  x_label = "Years since P1",
  width = 10.5,
  height = 4.1,
  filename_stub = "p1_bootstrap_transition_probabilities_from_P1",
  ncol = 3
)

#P2-only figure: exact filename already used in Overleaf
p_p2 <- make_probability_plot(
  c("P2_to_Rpost", "P2_to_D"),
  x_label = "Years since direct P2",
  width = 8.2,
  height = 4.2,
  filename_stub = "p1_bootstrap_transition_probabilities_from_P2",
  ncol = 2
)

#Rpre figure kept separate rather than mixing clocks/origins with P2
p_rpre <- make_probability_plot(
  c("Rpre_to_P2", "Rpre_to_D"),
  x_label = "Years since revision before P2",
  width = 8.2,
  height = 4.2,
  filename_stub = "p1_bootstrap_transition_probabilities_from_Rpre",
  ncol = 2
)

#main P1 risk differences only
p1_rd_data <- rd_data %>%
  filter(as.character(transition) %in% c("P1_to_P2", "P1_to_Rpre", "P1_to_D"))

p1_rd_limits <- range(
  c(p1_rd_data$rd_lower_95_percentage_points, p1_rd_data$rd_upper_95_percentage_points),
  na.rm = TRUE
)
p1_rd_pad <- max(0.25, diff(p1_rd_limits) * 0.08)

p_rd_p1 <- ggplot(p1_rd_data, aes(x = time_years, y = point_rd_percentage_points)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = rd_lower_95_percentage_points, ymax = rd_upper_95_percentage_points),
    width = 0.22,
    linewidth = 0.55
  ) +
  geom_point(size = 2.2) +
  facet_wrap(~ transition_label, ncol = 3) +
  scale_x_continuous(breaks = c(1, 3, 5, 10, 15)) +
  coord_cartesian(ylim = c(p1_rd_limits[1] - p1_rd_pad, p1_rd_limits[2] + p1_rd_pad)) +
  labs(
    x = "Years since P1",
    y = "Risk difference, first hip minus first knee (percentage points)"
  ) +
  common_prob_theme +
  theme(legend.position = "none")

save_plot(
  p_rd_p1,
  "p1_bootstrap_risk_differences_from_P1",
  width = 10.5,
  height = 4.1
)

#secondary risk differences from P2 and Rpre
secondary_rd_data <- rd_data %>%
  filter(as.character(transition) %in% c("P2_to_Rpost", "P2_to_D", "Rpre_to_P2", "Rpre_to_D"))

sec_rd_limits <- range(
  c(secondary_rd_data$rd_lower_95_percentage_points, secondary_rd_data$rd_upper_95_percentage_points),
  na.rm = TRUE
)
sec_rd_pad <- max(0.5, diff(sec_rd_limits) * 0.06)

p_rd_secondary <- ggplot(secondary_rd_data, aes(x = time_years, y = point_rd_percentage_points)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = rd_lower_95_percentage_points, ymax = rd_upper_95_percentage_points),
    width = 0.22,
    linewidth = 0.55
  ) +
  geom_point(size = 2.2) +
  facet_wrap(~ transition_label, ncol = 2) +
  scale_x_continuous(breaks = c(1, 3, 5, 10, 15)) +
  coord_cartesian(ylim = c(sec_rd_limits[1] - sec_rd_pad, sec_rd_limits[2] + sec_rd_pad)) +
  labs(
    x = "Years since entry into origin state",
    y = "Risk difference, first hip minus first knee (percentage points)"
  ) +
  common_prob_theme +
  theme(legend.position = "none")

save_plot(
  p_rd_secondary,
  "p1_bootstrap_risk_differences_secondary_origins",
  width = 8.5,
  height = 6.2
)

#overlap-weighted P1 probabilities
ow_prob <- ow_probability_ci %>%
  filter(method == "overlap_weighted_AJ") %>%
  mutate(
    transition_label = transition_label(transition),
    group_label = group_label(as.character(group)),
    transition = factor(transition, levels = c("P1_to_P2", "P1_to_Rpre", "P1_to_D")),
    group_label = factor(group_label, levels = c("First knee", "First hip"))
  )

p_ow <- ggplot(
  ow_prob,
  aes(x = time_years, y = point_estimate, linetype = group_label, shape = group_label, group = group_label)
) +
  geom_line(linewidth = 0.9, colour = "black") +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high, colour = group_label),
    width = 0.18,
    linewidth = 0.6,
    position = position_dodge(width = 0.22)
  ) +
  geom_point(size = 2.2, colour = "black", position = position_dodge(width = 0.22)) +
  facet_wrap(~ transition_label, scales = "free_y", ncol = 3) +
  scale_linetype_manual(values = line_values) +
  scale_shape_manual(values = shape_values) +
  scale_colour_manual(values = ci_grey_values, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.08))) +
  scale_x_continuous(breaks = c(1, 3, 5, 10, 15)) +
  labs(
    x = "Years since P1",
    y = "Overlap-weighted cumulative incidence"
  ) +
  common_prob_theme

save_plot(
  p_ow,
  "p1_overlap_weighted_transition_probabilities_from_P1",
  width = 10.5,
  height = 4.1
)

#main vs overlap risk-difference comparison
main_p1_rd <- p1_rd_data %>%
  transmute(
    transition = as.character(transition),
    transition_label,
    time_years,
    method = "Main standardized Cox",
    estimate = point_rd_percentage_points,
    lower = rd_lower_95_percentage_points,
    upper = rd_upper_95_percentage_points
  )

ow_p1_rd <- ow_contrast_ci %>%
  filter(method == "overlap_weighted_AJ") %>%
  transmute(
    transition,
    transition_label = transition_label(transition),
    time_years,
    method = "Overlap-weighted AJ",
    estimate = point_pp,
    lower = ci_low_pp,
    upper = ci_high_pp
  )

compare_rd <- bind_rows(main_p1_rd, ow_p1_rd) %>%
  mutate(
    transition = factor(transition, levels = c("P1_to_P2", "P1_to_Rpre", "P1_to_D")),
    method = factor(method, levels = c("Main standardized Cox", "Overlap-weighted AJ"))
  )

p_compare <- ggplot(
  compare_rd,
  aes(x = time_years, y = estimate, linetype = method, shape = method, group = method)
) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5) +
  geom_line(linewidth = 0.85, colour = "black", position = position_dodge(width = 0.25)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, colour = method),
    width = 0.18,
    linewidth = 0.55,
    position = position_dodge(width = 0.25)
  ) +
  geom_point(size = 2.1, colour = "black", position = position_dodge(width = 0.25)) +
  facet_wrap(~ transition_label, scales = "free_y", ncol = 3) +
  scale_linetype_manual(values = c("Main standardized Cox" = "solid", "Overlap-weighted AJ" = "longdash")) +
  scale_shape_manual(values = c("Main standardized Cox" = 16, "Overlap-weighted AJ" = 17)) +
  scale_colour_manual(values = c("Main standardized Cox" = "grey35", "Overlap-weighted AJ" = "grey65"), guide = "none") +
  scale_x_continuous(breaks = c(1, 3, 5, 10, 15)) +
  labs(
    x = "Years since P1",
    y = "Risk difference, first hip minus first knee (percentage points)"
  ) +
  common_prob_theme

save_plot(
  p_compare,
  "p1_main_vs_overlap_risk_differences_from_P1",
  width = 10.5,
  height = 4.1
)

expected_png <- c(
  "p1_bootstrap_transition_probabilities_from_P1.png",
  "p1_bootstrap_transition_probabilities_from_P2.png",
  "p1_bootstrap_transition_probabilities_from_Rpre.png",
  "p1_bootstrap_risk_differences_from_P1.png",
  "p1_bootstrap_risk_differences_secondary_origins.png",
  "p1_overlap_weighted_transition_probabilities_from_P1.png",
  "p1_main_vs_overlap_risk_differences_from_P1.png"
)
expected_pdf <- sub("\\.png$", ".pdf", expected_png)

qc <- tibble(
  check = c(
    "05 and 09 point estimates match",
    "main bootstrap 500 for every origin",
    "main probability CIs complete",
    "main risk-difference CIs complete",
    "overlap bootstrap QC all PASS",
    "overlap bootstrap 500 successful",
    "all final PNG figures written",
    "all final PDF figures written"
  ),
  value = c(
    max(point_check$abs_diff, na.rm = TRUE),
    all(origin_success$successful_replicates == 500),
    !any(is.na(ci_data$lower_95) | is.na(ci_data$upper_95)),
    !any(is.na(rd_data$rd_lower_95_percentage_points) | is.na(rd_data$rd_upper_95_percentage_points)),
    all(ow_qc$status == "PASS"),
    any(ow_success$method == "overlap_weighted_AJ" & ow_success$successful_replicates == 500),
    all(file.exists(file.path(fig_dir, expected_png))),
    all(file.exists(file.path(fig_dir, expected_pdf)))
  ),
  status = c(
    if_else(max(point_check$abs_diff, na.rm = TRUE) <= 1e-10, "PASS", "CHECK"),
    if_else(all(origin_success$successful_replicates == 500), "PASS", "CHECK"),
    if_else(!any(is.na(ci_data$lower_95) | is.na(ci_data$upper_95)), "PASS", "CHECK"),
    if_else(!any(is.na(rd_data$rd_lower_95_percentage_points) | is.na(rd_data$rd_upper_95_percentage_points)), "PASS", "CHECK"),
    if_else(all(ow_qc$status == "PASS"), "PASS", "CHECK"),
    if_else(any(ow_success$method == "overlap_weighted_AJ" & ow_success$successful_replicates == 500), "PASS", "CHECK"),
    if_else(all(file.exists(file.path(fig_dir, expected_png))), "PASS", "CHECK"),
    if_else(all(file.exists(file.path(fig_dir, expected_pdf))), "PASS", "CHECK")
  )
)

plot_data_xlsx <- file.path(fig_dir, "p1_final_figure_data_bw.xlsx")
wb <- createWorkbook()
plot_tabs <- list(
  qc = qc,
  main_probability_ci = ci_data,
  main_risk_difference = rd_data,
  overlap_probability_ci = ow_prob,
  main_vs_overlap_rd = compare_rd
)
for (nm in names(plot_tabs)) {
  addWorksheet(wb, substr(nm, 1, 31))
  writeData(wb, substr(nm, 1, 31), plot_tabs[[nm]])
  setColWidths(wb, substr(nm, 1, 31), cols = 1:50, widths = "auto")
}
saveWorkbook(wb, plot_data_xlsx, overwrite = TRUE)

cat("\n09b complete\n")
print(qc)
cat("\nfigures written to: ", fig_dir, "\n", sep = "")
cat("figure data: ", plot_data_xlsx, "\n", sep = "")