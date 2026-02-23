#!/usr/bin/env Rscript

# =============================
#   Load libraries
# =============================
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================
#   Command-line arguments
# =============================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_scfr.R <input_file.csv> <output_file.png/pdf>")
}

input_file  <- args[1]
output_file <- args[2]

# =============================
#   Read input
# =============================
df <- read.table(input_file, sep = ",", header = TRUE)

df$Window <- factor(df$Window, levels = df$Window)

# log10 transform
df$log10_total <- log10(df$Total_No_SCFR)

# Axis scaling
left_min <- 0
left_max <- max(df$log10_total, na.rm = TRUE)
scale_factor <- left_max / 100   # map % to left-axis height

# Scale % metrics for plotting
df$scaled_unfiltered <- df$Percent_unfiltered_by_filtered_SCFR * scale_factor
df$scaled_genome     <- df$Percent_genome_by_SCFR * scale_factor
df$scaled_by_coding  <- df$Percent_SCFR_by_coding * scale_factor

# Bar shading fractions
df$stack_coding    <- df$log10_total * df$Percent_coding_SCFR_count / 100
df$stack_noncoding <- df$log10_total * df$Percent_noncoding_SCFR_count / 100

df_long <- df %>%
  pivot_longer(cols = c(stack_coding, stack_noncoding),
               names_to = "type", values_to = "value")

# For percentage-line legend
df_line <- bind_rows(
  df %>% mutate(metric = "unfiltered_vs_filtered",
                scaled_value = scaled_unfiltered),
  df %>% mutate(metric = "genome_by_SCFR",
                scaled_value = scaled_genome),
  df %>% mutate(metric = "SCFR_by_coding",
                scaled_value = scaled_by_coding)
)

metric_labels <- c(
  "unfiltered_vs_filtered" = "% unfiltered SCFR covered by filtered SCFR",
  "genome_by_SCFR"         = "% genome covered by SCFR",
  "SCFR_by_coding"         = "% SCFR covered by coding exons"
)

metric_colors <- c(
  "unfiltered_vs_filtered" = "darkblue",
  "genome_by_SCFR"         = "darkgreen",
  "SCFR_by_coding"         = "darkred"
)

# Species from DF for title
species_name <- unique(df$Species)

# =============================
#   Build plot
# =============================
p <- ggplot() +
  
  # Stacked bars (log10 total SCFR count)
  geom_bar(data = df_long,
           aes(x = Window, y = value, fill = type),
           stat = "identity", width = 0.8) +
  
  scale_fill_manual(values = c("stack_coding" = "orange",
                               "stack_noncoding" = "grey70"),
                    name = "SCFR fraction",
                    labels = c("Coding SCFR %", "Non-coding SCFR %")) +
  
  # Percentage lines
  geom_line(data = df_line,
            aes(x = Window, y = scaled_value, color = metric, group = metric),
            linewidth = 1) +
  
  geom_point(data = df_line,
             aes(x = Window, y = scaled_value, color = metric),
             size = 2) +
  
  scale_color_manual(
    name = "",
    values = metric_colors,
    labels = metric_labels
  ) +
  
  scale_y_continuous(
    name = "log10(Total SCFR count)",
    limits = c(left_min, left_max),
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Percentage (%)")
  ) +
  
  labs(x = "Window size",
       title = paste0(species_name)
  ) +

  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.box = "vertical"
  ) +

guides(
    fill = guide_legend(
      title = "SCFR fraction",
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.3, "cm"),  # Increase height of legend boxes
      keywidth  = unit(0.3, "cm")   # Increase width of legend boxes

    ),
    color = guide_legend(
      title = "Percentage (%)",
      ncol = 1,
      byrow = TRUE
    )
  ) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.direction = "horizontal",

    # left-align legend text inside keys
    legend.text = element_text(size = 8, hjust = 0),
    legend.title = element_text(size = 8, hjust = 0),
    legend.spacing.y = unit(0.002, "cm"),
    legend.key.size = unit(0.01, "cm")
  )

# =============================
#   Save output
# =============================
if (grepl("\\.png$", output_file, ignore.case = TRUE)) {
  ggsave(output_file, p, width = 12, height = 7, dpi = 300)
} else if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
  ggsave(output_file, p, width = 7, height = 7)
} else {
  stop("Output file must end in .png or .pdf")
}


