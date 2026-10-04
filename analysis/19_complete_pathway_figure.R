source("config.R")

# ============================================================
# figure 2.4: complete P1 -> P2_direct -> Rpost pathway
# black-and-white thesis version
# point estimates only
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)
library(openxlsx)

out_dir <- output_dir
fig_dir <- file.path(out_dir, "figures_final")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

#find final pathway workbook
candidate_files <- c(
  file.path(out_dir, "p1_pathway_P1_P2_Rpost_point.xlsx"),
  file.path(out_dir, "p1_pathway_P1_P2_Rpost_point(1).xlsx"),
  file.path(out_dir, "p1_pathway_P1_P2_Rpost_point(2).xlsx")
)

xlsx_file <- candidate_files[file.exists(candidate_files)][1]

if (is.na(xlsx_file)) {
  stop("final pathway workbook not found")
}

message("using: ", xlsx_file)

#read final pathway estimates
report <- read.xlsx(
  xlsx_file,
  sheet = "pathway_report"
)

curve <- read.xlsx(
  xlsx_file,
  sheet = "pathway_curve"
)

#standardize labels
clean_group <- function(x) {
  
  x <- as.character(x)
  
  case_when(
    x %in% c("first_knee", "First knee", "knee", "0") ~ "First knee",
    x %in% c("first_hip", "First hip", "hip", "1") ~ "First hip",
    TRUE ~ x
  )
}

report <- report %>%
  mutate(
    group_label = factor(
      clean_group(group),
      levels = c("First knee", "First hip")
    )
  ) %>%
  filter(time_years %in% c(1, 3, 5, 10, 15))

curve <- curve %>%
  mutate(
    group_label = factor(
      clean_group(group),
      levels = c("First knee", "First hip")
    )
  ) %>%
  filter(
    time_years >= 0,
    time_years <= 15
  )

#black-and-white figure
p <- ggplot(
  curve,
  aes(
    x = time_years,
    y = probability,
    linetype = group_label,
    group = group_label
  )
) +
  
  geom_line(
    linewidth = 0.9,
    colour = "black"
  ) +
  
  geom_point(
    data = report,
    aes(
      x = time_years,
      y = probability,
      shape = group_label
    ),
    inherit.aes = FALSE,
    size = 2.5,
    colour = "black"
  ) +
  
  scale_linetype_manual(
    values = c(
      "First knee" = "solid",
      "First hip" = "longdash"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "First knee" = 16,
      "First hip" = 17
    )
  ) +
  
  scale_x_continuous(
    breaks = c(0, 1, 3, 5, 10, 15),
    limits = c(0, 15),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 0.5),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  
  labs(
    x = "Years since P1",
    y = "Complete pathway probability",
    linetype = NULL,
    shape = NULL
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.25
    ),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.text = element_text(colour = "black")
  )

#same filename already used in Overleaf
png_file <- file.path(
  fig_dir,
  "p1_pathway_P1_P2_Rpost.png"
)

pdf_file <- file.path(
  fig_dir,
  "p1_pathway_P1_P2_Rpost.pdf"
)

ggsave(
  png_file,
  p,
  width = 7.2,
  height = 4.6,
  dpi = 400
)

ggsave(
  pdf_file,
  p,
  width = 7.2,
  height = 4.6
)

cat("\nfigure written to:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")

cat("\npoint estimates used:\n")

print(
  report %>%
    select(
      time_years,
      group_label,
      probability
    ) %>%
    arrange(
      time_years,
      group_label
    )
)