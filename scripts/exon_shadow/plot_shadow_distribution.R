#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(forcats)

# -----------------------------
# Input
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript plot_shadow_distribution.R <input_file.tsv> <output_pdf>")
}

input_file <- args[1]
output_pdf <- args[2]

# -----------------------------
# Helper: format bp / kb / mb
# -----------------------------
# Function to format length values into bp, kb, or mb
format_length <- function(x) {
  x <- as.numeric(x)

  if (x >= 1000000) {
    return(paste0(round(x / 1000000, 1), " mb"))
  } else if (x >= 1000) {
    return(paste0(round(x / 1000, 1), " kb"))
  } else {
    return(paste0(x, " bp"))
  }
}

species_order <- c(
  "gibbon",
  "gorilla",
  "human",
  "bonobo",
  "chimpanzee",
  "bornean orangutan",
  "sumatran orangutan"
)

# -----------------------------
# Read data
# -----------------------------
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

df <- df %>%
  mutate(
    Species = factor(Species, levels = species_order),    
    SCFR_type = factor(
      SCFR_type,
      levels = c("single_exon", "multi_exon", "composite_exon", "all")
    ),
    
    min_lab  = sapply(min, format_length),
    max_lab  = sapply(max, format_length),
#    mean_lab = sapply(mean, format_length),
    N_lab    = paste0("N=", format(N, big.mark = ",")),
    
    zero_offset = 0.8,
    min_plot    = ifelse(min <= 0, zero_offset, min),
    max_plot    = ifelse(max <= 0, zero_offset, max),
    mean_plot   = ifelse(mean <= 0, zero_offset, mean),
    median_plot = ifelse(median <= 0, zero_offset, median)
  )

global_top_y <- max(df$max_plot) * 2

df_N <- df %>%
  group_by(Species) %>%
  summarise(
    N_label = paste0("N=", format(sum(N), big.mark = ","))
  )

# -----------------------------
# Geometry constants
# -----------------------------
box_width    <- 0.65
line_width   <- box_width * 0.75
cap_width    <- box_width
label_offset <- 0.2
med_mean_offset <- -0.2

dodge <- position_dodge(width = 0.8)

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(df, aes(x = Species, fill = SCFR_type )) +
  
  # Base box
  geom_boxplot(
    aes(
      ymin   = min_plot,
      lower  = Q1,
      middle = median_plot,
      upper  = Q3,
      ymax   = max_plot
    ),
    stat = "identity",
    width = box_width,
    color = "black",
    linewidth = 0.35,
    position = dodge
  ) +
  
  # Median line
  geom_errorbar(
    aes(
      y = median_plot,
      ymin = median_plot,
      ymax = median_plot,
      color = "Median"
    ),
    width = line_width,
    linewidth = 1.0,
    position = dodge
  ) +
  
  # Median label
  geom_text(
    aes(
      y = median_plot * (1 + med_mean_offset),
      label = paste0(median)
    ),
    size = 2.5,
    color = "darkgreen",
    position = dodge
  ) +
  
  # Mean line
  geom_errorbar(
    aes(
      y = mean_plot,
      ymin = mean_plot,
      ymax = mean_plot,
      color = "Mean"
    ),
    width = line_width,
    linewidth = 0.8,
    linetype = "dotted",
    position = dodge
  ) +
  
  # Mean point (no legend)
  geom_point(
    aes(y = mean_plot),
    shape = 16,
    size = 2,
    color = "blue",
    position = dodge,
    show.legend = FALSE
  ) +
  
  # Mean label
  geom_text(
    aes(
      y = mean_plot * (1 + med_mean_offset),
      label = paste0(mean)
    ),
    size = 2.5,
    color = "blue",
    position = dodge
  ) +
  
  # Min / Max caps
  geom_errorbar(
    aes(
      y = min_plot,
      ymin = min_plot,
      ymax = min_plot,
      color = "Min / Max"
    ),
    width = cap_width,
    linewidth = 0.9,
    position = dodge
  ) +
  geom_errorbar(
    aes(
      y = max_plot,
      ymin = max_plot,
      ymax = max_plot,
      color = "Min / Max"
    ),
    width = cap_width,
    linewidth = 0.9,
    position = dodge
  ) +
  
  # Min / Max labels
  geom_text(
    aes(y = max_plot * (1 + label_offset * 1.5), label = max_lab),
    size = 2.5,
    position = dodge
  ) +
  geom_text(
    aes(y = min_plot * (1 - label_offset), label = min_lab),
    size = 2.5,
    position = dodge
  ) +
  
  # N above species
  geom_text(
    data = df_N,
    aes(
      x = Species,
      y = global_top_y,
      label = N_label
    ),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  
  scale_y_log10() +
  
  scale_fill_brewer(palette = "Set2", name = "SCFR type") +
  
  scale_color_manual(
    name = "Summary statistics",
    values = c(
      "Mean"      = "blue",
      "Median"    = "darkgreen",
      "Min / Max" = "red"
    ),
    breaks = c("Mean", "Median", "Min / Max"),
    guide = guide_legend(
      override.aes = list(
        linetype = c("dotted", "solid", "solid"),
        linewidth = c(0.8, 1.0, 0.9)
      )
    )
  ) +
  
  labs(
    x = "",
    y = "Shadow length (bp)",
    title = "SCFR shadow length distribution across primates"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 30, 10)
  )

# -----------------------------
# Draw
# -----------------------------
#print(p)

ggsave(
  output_pdf,
  p,
  width = 13,
  height = 7,
  dpi = 800,
  device = cairo_pdf
)

cat("Saved:", output_pdf, "\n")
