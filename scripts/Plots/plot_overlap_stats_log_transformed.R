#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# ---------------------------------------------------------
# Input file
# ---------------------------------------------------------
if (length(args) < 2) {
  stop("Usage: Rscript plot.R <input.txt> <output.pdf>")
}

input_file  <- args[1]
output_pdf  <- args[2]
species <- args[3]

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

# -----------------------------------------------------------
# Read input file and prepare data
# -----------------------------------------------------------
df <- read.table(input_file, header = TRUE, sep = "\t")

df$Length_threshold <- factor(df$Length_threshold, levels = df$Length_threshold[order(as.numeric(as.character(df$Length_threshold)))])

# Apply the formatting function
df$max_kb <- sapply(df$max, format_length)
df$min_kb <- sapply(df$min, format_length)
df$mean_kb <- sapply(df$mean, format_length)

num_categories <- length(levels(df$Length_threshold))

# ------------------------------------
# Build the plot
# ------------------------------------
p <- ggplot(df, aes(x = Length_threshold)) +
  
  # 1. Base Boxplot 
  geom_boxplot(
    aes(ymin = min, lower = q1, middle = median, upper = q3, ymax = max),
    stat="identity",
    fill="grey85",
    color="black",
    linewidth=0.4,
    width=0.5,
    outlier.shape = NA
  ) +
  
  # 2. Q1 / Q3 Edges (Thicker Black)
  geom_segment(aes(x=as.numeric(Length_threshold)-0.25, xend=as.numeric(Length_threshold)+0.25, y=q1, yend=q1, color="Q1 / Q3"), linewidth=1.2) +
  geom_segment(aes(x=as.numeric(Length_threshold)-0.25, xend=as.numeric(Length_threshold)+0.25, y=q3, yend=q3, color="Q1 / Q3"), linewidth=1.2) +
  
  # 3. Median Line (Dark Green)
  geom_segment(aes(x=as.numeric(Length_threshold)-0.22, xend=as.numeric(Length_threshold)+0.22, y=median, yend=median, color="Median"), linewidth=1.2) +
  
  # 4. Min/Max Caps (Red)
  geom_segment(aes(x=as.numeric(Length_threshold)-0.12, xend=as.numeric(Length_threshold)+0.12, y=min, color="Min / Max"), linewidth=0.8) +
  geom_segment(aes(x=as.numeric(Length_threshold)-0.12, xend=as.numeric(Length_threshold)+0.12, y=max, color="Min / Max"), linewidth=0.8) +
  
  # 5. Mean Line (Dashed Light Curve) - Background
  geom_line(aes(y = mean, group = 1), color = "darkgrey", linetype = "dashed", linewidth = 0.5) +
  
  # 6. Mean Dotted Line + Point (Blue)
  geom_segment(aes(x=as.numeric(Length_threshold)-0.22, xend=as.numeric(Length_threshold)+0.22, y=mean, yend=mean, color="Mean"), linetype="dotted", linewidth=0.7) +
  geom_point(aes(y=mean, color="Mean"), size=2) +
  
  # 7. Color Mapping for Legend Keys
  scale_color_manual(
    name = "", 
    values = c("Mean"="blue", "Median"="darkgreen", "Q1 / Q3"="black", "Min / Max"="red"),
    breaks = c("Mean", "Median", "Q1 / Q3", "Min / Max"),
    guide = guide_legend(
      override.aes = list(
        linetype = c("dotted", "solid", "solid", "solid"),
        shape = c(19, NA, NA, NA),
        linewidth = c(0.7, 1.2, 1.2, 1.2),
        color = c("blue", "darkgreen", "black", "red")
      )
    )
  ) +
  
  # 8. Labels (Fixing the Mean Label Position)
  
  # Mean label: Left of the box, significantly above the mean point
  geom_text(aes(
    # Position X: Left of the box
    x = (as.numeric(Length_threshold) * 0.999 ) - 0.1 , 
    # Position Y: Above the mean point (e.g., mean * 1.6 - a high offset on log scale)
    y = (mean * 1.3) + 50, 
    label = paste0(mean_kb)
  ), color="blue", size=3.5, hjust=1) + 
  
  # Max label: Significantly higher than the new mean label position
  geom_text(aes(
    y = max * 2, # Increased offset to ensure it's higher than the mean label
    label = paste0(max_kb)
  ), size=3.5) +
  
  # Min label: Slightly lower
  geom_text(aes(
    y = min * 0.75,
    label = paste0(min_kb)
  ), size=3.5) +
  
  scale_y_log10() +
  
  labs(
    x = "Length Threshold (bp)",
    y = "SCFR - Gene Desert Overlap (bp)",
    title = species
  ) +
  
  coord_cartesian(xlim = c(0.1, num_categories + 0.5), expand = TRUE) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1, size=11),
    axis.text.y = element_text(hjust=1, size=11),
    legend.position = "right", 
    legend.title = element_blank()
  )

#print(p)

# --------------------------------------------------------------
# SAVE as high-resolution PDF (No changes needed here)
# --------------------------------------------------------------

ggsave(
  filename = output_pdf,
  plot = p,
  width = 9,
  height = 8,
  units = "in",
  dpi = 800,
  device = cairo_pdf
)

cat("Saved PDF:", output_pdf, "\n")

