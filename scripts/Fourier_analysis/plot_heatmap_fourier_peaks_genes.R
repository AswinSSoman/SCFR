#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# ---------------------------------------------------------
# Input file
# ---------------------------------------------------------
if (length(args) < 2) {
  stop("Usage: Rscript plot.R <input.txt> <output.pdf> <X-bin_size> <Y-bin_size>")
}


# --- 1. Load Libraries ---
# Make sure you have the tidyverse package installed: install.packages("tidyverse")
library(tidyverse)

# --- 2. Define File and Bin Parameters ---
file_name <- args[1]
output <- args[2]
x_bin_size <- 0.01  # Bin size for Frequency (X-axis)
y_bin_size <- 100   # Bin size for Magnitude (Y-axis) - set smaller to highlight other peaks

# --- 3. Load Data ---
# Ensure the 'gene_fft_data.csv' file is in your working directory
data <- read.table(file_name, header = T, sep = "\t")

# --- 4. Bin Data and Count Genes ---
# --- 4. Bin Data and Count Genes ---
heatmap_data <- data %>%
  # Create binned variables. 'right = FALSE' ensures bins are [start, end)
  mutate(
    # Bins for Frequency (X)
    Freq_Binned = cut(
      Frequency,
      # Generate breaks from min(Frequency) down to 0, and up to max(Frequency)
      breaks = seq(floor(min(Frequency, 0)), ceiling(max(Frequency) + x_bin_size), by = x_bin_size),
      include.lowest = TRUE,
      right = FALSE
    ),
    # Bins for Magnitude (Y)
    Mag_Binned = cut(
      Magnitude,
      # Generate breaks from a reasonable minimum (e.g., 0) up to max(Magnitude)
      breaks = seq(floor(min(Magnitude, 0)), ceiling(max(Magnitude) + y_bin_size), by = y_bin_size),
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  # Group by the created bins and count the number of unique genes in each bin
  group_by(Freq_Binned, Mag_Binned) %>%
  summarise(
    N_Genes = n_distinct(Gene),
    .groups = 'drop'
  )

# --- 5. Plot the Heatmap ---
p <- ggplot(heatmap_data, aes(x = Freq_Binned, y = Mag_Binned, fill = N_Genes)) +
  # Add the tiles for the heatmap
  geom_tile(color = "gray", size = 0.5) + # Add gray borders for better separation
  # Set the color gradient from blue (low) to red (high)
  scale_fill_gradient(
    low = "blue",
    high = "red",
    name = "Number of Peaks" # Legend title for number of genes
  ) +
  # Customize axis labels and plot title
  labs(
    title = "Heatmap of Gene Fourier Peaks: Frequency vs. Magnitude",
    subtitle = paste(
      "Frequency Bin Size (X): ", x_bin_size, " | Magnitude Bin Size (Y): ", y_bin_size,
      sep = ""
    ),
    x = "Frequency Bin",
    y = "Magnitude Bin"
  ) +
  # Optional: Improve appearance by adjusting theme and rotating X-axis labels for readability
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

# --- 6. Display/Save the Plot ---
#print(p)

# To save the plot to a file (e.g., PNG)
#ggsave("gene_heatmap.png", plot = p, width = 10, height = 8)
ggsave(
  filename = output,
  plot = p,
  width =12,
  height = 9,
  #   dpi = 600,       # High resolution
  #   device = cairo_pdf    # Best for vector export
  units = "in",
  dpi = 800,
  device = cairo_pdf
)

cat("Saved PDF:", output, "\n")

# --- 7. Print Bin Sizes for user confirmation ---
#cat("\n--- Binning Information ---\n")
#cat("Frequency (X-axis) Bin Size: ", x_bin_size, "\n")
#cat("Magnitude (Y-axis) Bin Size: ", y_bin_size, "\n")
#cat("Note: The smaller Magnitude bin size (", y_bin_size, ") is intended to resolve distinct magnitude peaks.\n")

