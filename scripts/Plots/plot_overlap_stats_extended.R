#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(cowplot)

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
  # Convert to numeric in case it's not already
  x <- as.numeric(x)
  
  # Check for MB (1,000,000 bp or more)
  if (x >= 1000000) {
    return(paste0(round(x / 1000000, 1), " mb"))
  } 
  # Check for KB (1,000 bp or more, but less than 1,000,000 bp)
  else if (x >= 1000) {
    return(paste0(round(x / 1000, 1), " kb"))
  } 
  # Less than 1,000 bp
  else {
    return(paste0(x, " bp"))
  }
}

# -----------------------------------------------------------
# Read input file
# -----------------------------------------------------------
df <- read.table(input_file, header = TRUE, sep = "\t")

# --- MODIFICATION START ---
# New data uses 'Length_threshold' instead of 'species'.
# We need to ensure 'Length_threshold' is treated as a factor for the x-axis
# but maintain the numerical order.

# The new data has a 'Length_threshold' column, which should be the X-axis.
# Convert to a factor for discrete plotting, ordering numerically.
df$Length_threshold <- factor(df$Length_threshold, levels = df$Length_threshold[order(as.numeric(as.character(df$Length_threshold)))])

# Apply the formatting function
df$max_kb <- sapply(df$max, format_length)
df$min_kb <- sapply(df$min, format_length)
df$mean_kb <- sapply(df$mean, format_length)

# --- MODIFICATION END ---
num_categories <- length(levels(df$Length_threshold))
# ------------------------------------
# Build main plot WITHOUT legend
# ------------------------------------
p <- ggplot(df, aes(
  x = Length_threshold,
  ymin = min, lower = q1, middle = median, upper = q3, ymax = max
)) +
  
  # Grey box
  geom_boxplot(
    stat="identity",
    fill="grey85",
    color="black",
    linewidth=0.4,
    width=0.5,
    outlier.shape = NA
  ) +

  # --- NEW LAYER: Dashed Light Curve through the means ---
  geom_line(
    aes(y = mean, group = 1), # 'group = 1' is essential for connecting points across categorical X-axis
    color = "darkgrey",       # Light curve color
    linetype = "dashed",      # Dashed line style
    linewidth = 0.5           # Thin line weight
  ) +
  
  # Thicker Q1/Q3 edges
  geom_segment(aes(x=as.numeric(Length_threshold)-0.25,
                   xend=as.numeric(Length_threshold)+0.25,
                   y=q1, yend=q1),
               linewidth=1.2) +
  geom_segment(aes(x=as.numeric(Length_threshold)-0.25,
                   xend=as.numeric(Length_threshold)+0.25,
                   y=q3, yend=q3),
               linewidth=1.2) +
  
  # Thin whiskers + caps
  geom_segment(aes(x=as.numeric(Length_threshold), xend=as.numeric(Length_threshold),
                   y=min, yend=q1), color="black", linewidth=0.4) +
  geom_segment(aes(x=as.numeric(Length_threshold)-0.12, xend=as.numeric(Length_threshold)+0.12,
                   y=min), color="red", linewidth=0.8) +
  geom_segment(aes(x=as.numeric(Length_threshold), xend=as.numeric(Length_threshold),
                   y=q3, yend=max), color="black", linewidth=0.4) +
  geom_segment(aes(x=as.numeric(Length_threshold)-0.12, xend=as.numeric(Length_threshold)+0.12,
                   y=max), color="red", linewidth=0.8) +
  
  # Median (dark green)
  geom_segment(aes(
    x = as.numeric(Length_threshold)-0.22,
    xend = as.numeric(Length_threshold)+0.22,
    y = median, yend = median
  ), color="darkgreen", linewidth=1.2) +
  
  # Mean (blue dotted) + point
  geom_point(aes(y=mean, color="Mean"), size=2) +
  geom_segment(aes(
    x=as.numeric(Length_threshold)-0.22, xend=as.numeric(Length_threshold)+0.22,
    y=mean, yend=mean, color="Mean"
  ), linewidth=0.5, linetype="dotted") +
  
  scale_color_manual(values=c("Mean"="blue")) +
  
  # Labels --------------------------------

# Median (lower half)

  # Mean label
  geom_text(aes(
    x = as.numeric(Length_threshold) + 0.35,
    y = q3 + 200, # Adjusted Y offset
    label = paste0(mean_kb)
  ), color="blue", size=3.5) +
  
  # Min label
  # --- MODIFICATION: Changed 'species' to 'Length_threshold' and adjusted Y offset for new data scale ---
  geom_text(aes(
    y = min -250, # Adjusted Y offset
    label = paste0(min_kb)
  ), size=3.5) +
  
  # Max label (bp now, not kb)
  geom_text(aes(
    y = max +550, # Adjusted Y offset
    label = paste0(max_kb)
  ), size=3.5) +
  
    # --- MODIFICATION: Updated titles and axis labels for the new data context ---
  labs(
    x = "Length Threshold (bp)",
    y = "SCFR - Gene Desert Overlap (bp)",
    title = species,
    color = ""  # Legend title empty
  ) +
  
  coord_cartesian(xlim = c(0.1, num_categories + 0.8),
    ylim = c(min(df$min) - 1000, max(df$max) + 1000), # Assuming min(df$min) is 1 and max(df$max) is 6801
    expand = FALSE
  ) +

  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1, size=11),
    axis.text.y = element_text(hjust=1, size=11),
    legend.position = "none"  # We'll place legend externally
  )

# ----------------------------------------------------------------
# Build external legend (No changes needed here as it's static)
# ----------------------------------------------------------------

legend_plot <- ggplot() +
  geom_segment(aes(x=1, xend=1.4, y=6.2, yend=6.2), color="red", linewidth=1.2) +
  geom_text(aes(x=1.45, y=6.2, label="Min / Max"), hjust=0, size=3.5)

legend_plot <- legend_plot +
  geom_segment(aes(x=1, xend=1.4, y=6.4, yend=6.4), linewidth=1.2) +
  geom_text(aes(x=1.45, y=6.4, label="Q1 / Q3"), hjust=0, size=3.5)

legend_plot <- legend_plot +
  geom_segment(aes(x=1, xend=1.4, y=6.6, yend=6.6), color="darkgreen", linewidth=1.2) +
  geom_text(aes(x=1.45, y=6.6, label="Median"), hjust=0, size=3.5)

legend_plot <- legend_plot +
  geom_segment(aes(x=1, xend=1.4, y=6.8, yend=6.8), color="blue", linetype="dotted", linewidth=0.7) +
  geom_text(aes(x=1.45, y=6.8, label="Mean"), hjust=0, size=3.5)

legend_plot <- legend_plot +
  xlim(1, 3) + ylim(1, 7) +
  theme_void()

# --------------------------------------------------------------
# Combine main plot + external legend side by side (No changes needed here)
# --------------------------------------------------------------

final_plot <- plot_grid(
  p, legend_plot,
  rel_widths = c(6, 1),
  nrow = 1
)

#print(final_plot)

# --------------------------------------------------------------
# SAVE as high-resolution PDF (No changes needed here)
# --------------------------------------------------------------

ggsave(
  filename = output_pdf,
  plot = final_plot,
  width = 9,
  height = 8,
  units = "in",
  dpi = 800,
  device = cairo_pdf
)

cat("Saved PDF:", output_pdf, "\n")
