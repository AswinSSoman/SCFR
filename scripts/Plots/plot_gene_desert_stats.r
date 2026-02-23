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

# Set the desired order for the species factor on the X-axis
species_order <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean orangutan",
  "Sumatran orangutan"
)
df$species <- factor(df$species, levels = species_order)
# Apply the formatting function
df$max_kb <- sapply(df$max, format_length)
df$min_kb <- sapply(df$min, format_length)
df$mean_kb <- sapply(df$mean, format_length)
df$median_kb <- sapply(df$median, format_length)
df$q1_kb <- sapply(df$q1, format_length)
df$q3_kb <- sapply(df$q3, format_length)

# ------------------------------------
# Build main plot WITHOUT legend
# ------------------------------------

p <- ggplot(df, aes(
  x = species,
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
    linewidth = 0.4           # Thin line weight
  ) +
  
  # Thicker Q1/Q3 edges
  geom_segment(aes(x=as.numeric(species)-0.25,
                   xend=as.numeric(species)+0.25,
                   y=q1, yend=q1),
               linewidth=1.2) +
  geom_segment(aes(x=as.numeric(species)-0.25,
                   xend=as.numeric(species)+0.25,
                   y=q3, yend=q3),
               linewidth=1.2) +
  
  # Thin whiskers + caps
  geom_segment(aes(x=as.numeric(species), xend=as.numeric(species),
                   y=min, yend=q1), color="black", linewidth=0.4) +
  geom_segment(aes(x=as.numeric(species)-0.12, xend=as.numeric(species)+0.12,
                   y=min), color="red", linewidth=0.8) +
  geom_segment(aes(x=as.numeric(species), xend=as.numeric(species),
                   y=q3, yend=max), color="black", linewidth=0.4) +
  geom_segment(aes(x=as.numeric(species)-0.12, xend=as.numeric(species)+0.12,
                   y=max), color="red", linewidth=0.8) +
  
  # Median (dark green)
  geom_segment(aes(
    x = as.numeric(species)-0.22,
    xend = as.numeric(species)+0.22,
    y = median, yend = median
  ), color="darkgreen", linewidth=1.2) +
  
  # Mean (blue dotted) + point
  geom_point(aes(y=mean, color="Mean"), size=2) +
  geom_segment(aes(
    x=as.numeric(species)-0.22, xend=as.numeric(species)+0.22,
    y=mean, yend=mean, color="Mean"
  ), linewidth=0.5, linetype="dotted") +
  
  scale_color_manual(values=c("Mean"="blue")) +
  
  # Labels --------------------------------

# Median (lower half)
geom_text(aes(
  y = q1 + (median-q1)*0.5,
  label = paste0(median_kb)
), size=3.5, color="darkgreen") +
  
  # Q1 (medium spacing)
  geom_text(aes(
    x = as.numeric(species)-0.28,
    y = q1 - q1 * 0.1,
    label = paste0(q1_kb)
  ), size=3.5) +
  
  # Q3 (medium spacing)
  geom_text(aes(
    x = as.numeric(species)-0.25,
    y = q3 + q3 * 0.20,
    label = paste0(q3_kb)
  ), size=3.5) +
  
  # Mean label — HIGHER OFFSET as requested
  geom_text(aes(
    x = as.numeric(species)+0.25,
    y =  q3 + q3 * 0.20,   # <-- Increased offset
    label = paste0(mean_kb)
  ), color="blue", size=3.5) +
  
  # Min label
  geom_text(aes(
    y = min - min*0.1,
    label = paste0(min_kb)
  ), size=3.5) +
  
  # Max label (kb)
  geom_text(aes(
    y = max + max*0.1,
    label = paste0(max_kb)
  ), size=3.5) +
  
  scale_y_log10() +
  
  labs(
    x = "Species",
    y = "Length (log10 bp)",
    title = "Gene Deserts Length Distribution Across Primate Genomes",
    color = ""   # Legend title empty
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=35, hjust=1, size=12),
    axis.text.y = element_text(hjust=1, size=11),
    legend.position = "none"   # We'll place legend externally
  )

# ----------------------------------------------------------------
# Build external legend
# ----------------------------------------------------------------

legend_plot <- ggplot() +
  geom_segment(aes(x=1, xend=1.4, y=6.2, yend=6.2), color="red", linewidth=1.2) +
  geom_text(aes(x=1.45, y=6.2, label="Min / Max"), hjust=0, size=4)

legend_plot <- legend_plot +
  geom_segment(aes(x=1, xend=1.4, y=6.4, yend=6.4), linewidth=1.2) +
  geom_text(aes(x=1.45, y=6.4, label="Q1 / Q3"), hjust=0, size=4)

legend_plot <- legend_plot +
  geom_segment(aes(x=1, xend=1.4, y=6.6, yend=6.6), color="darkgreen", linewidth=1.2) +
  geom_text(aes(x=1.45, y=6.6, label="Median"), hjust=0, size=4)

legend_plot <- legend_plot +
  geom_segment(aes(x=1, xend=1.4, y=6.8, yend=6.8), color="blue", linetype="dotted", linewidth=0.7) +
  geom_text(aes(x=1.45, y=6.8, label="Mean"), hjust=0, size=4)

legend_plot <- legend_plot +
  xlim(1, 3) + ylim(1, 7) +
  theme_void()

# --------------------------------------------------------------
# Combine main plot + external legend side by side
# --------------------------------------------------------------

final_plot <- plot_grid(
  p, legend_plot,
  rel_widths = c(6, 1),
  nrow = 1
)

#print(final_plot)

# --------------------------------------------------------------
# SAVE as high-resolution PDF
# --------------------------------------------------------------

ggsave(
  filename = output_pdf,
  plot = final_plot,
  width =12,
  height = 9,
  #   dpi = 600,       # High resolution
  #   device = cairo_pdf    # Best for vector export
  units = "in",
  dpi = 800,
  device = cairo_pdf
)

cat("Saved PDF:", output_pdf, "\n")

