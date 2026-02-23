#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot_downstream_exon_feature_heatmap.R <input.tsv> <BIN_WIDTH> <MAX_LEN> <output.pdf>")
}

# ============================
# PARAMETERS (ONLY CHANGE THIS)
# ============================
BIN_WIDTH <- as.numeric(args[2])
MAX_LEN   <- as.numeric(args[3])

# ============================
# Read data
# ============================

# ---------------------------------------------------------
# Input file
# ---------------------------------------------------------

input_file  <- args[1]
output_pdf <- args[4]

df <- read_tsv(input_file, col_types = cols())

# ============================
# Filter
# ============================
df <- df %>%
  filter(
    downstream_len_in_scfr >= -2,
    downstream_len_in_scfr <= MAX_LEN
  )

# ============================
# Create downstream bins
# ============================

# breaks: -2 to 0, then 1 onward
bin_breaks <- c(
  -2, 1,
  seq(1 + BIN_WIDTH, MAX_LEN + BIN_WIDTH, by = BIN_WIDTH)
)

bin_labels <- c(
  "-2–0",
  paste0(
    seq(1, MAX_LEN, by = BIN_WIDTH),
    "–",
    seq(BIN_WIDTH, MAX_LEN + BIN_WIDTH - 1, by = BIN_WIDTH)
  )
)

df <- df %>%
  mutate(
    downstream_bin = cut(
      downstream_len_in_scfr,
      breaks = bin_breaks,
      include.lowest = TRUE,
      right = FALSE,
      labels = bin_labels
    ),
    downstream_bin = factor(downstream_bin, levels = bin_labels)
  )

# ============================
# Categories for heatmap
# ============================

df <- df %>%
  mutate(
    exon_count_class = ifelse(exon_count > 1, "multi", "single"),
    exon_order = last_exon_order,
    splicing = last_exon_splicing
  )

# ============================
# Heatmap data
# ============================

heatmap_df <- df %>%
  select(downstream_bin, exon_count_class, exon_order, splicing) %>%
  pivot_longer(
    cols = c(exon_count_class, exon_order, splicing),
    names_to = "group",
    values_to = "category"
  ) %>%
  count(downstream_bin, category) %>%
  complete(
    downstream_bin,
    category = factor(
      category,
      levels = c(
        "alternative", "constitutive", "unique",
        "single", "multi",
        "first", "middle", "last"
      )
    ),
    fill = list(n = 0)
  )

# ============================
# Bin counts + percentages
# ============================

bin_stats <- df %>%
  count(downstream_bin, name = "bin_count") %>%
  mutate(
    pct = bin_count / sum(bin_count) * 100,
    label = paste0(bin_count, " (", sprintf("%.2f", pct), "%)")
  )

# ============================
# Heatmap plot
# ============================
heatmap_df <- heatmap_df %>%
  mutate(
    category = factor(
      category,
      levels = c(
        "alternative", "constitutive", "unique",
        "single", "multi",
        "first", "middle", "last"
      )
    )
  )

p_heatmap <- ggplot(heatmap_df, aes(x = category, y = downstream_bin, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Count") +
  labs(x = NULL, y = "Downstream length bins (bp)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 15)),
    panel.grid = element_blank()
  )

# ============================
# Count + percentage text plot 
# ============================

p_counts <- ggplot(bin_stats, aes(y = downstream_bin, x = 0, label = label)) +
  geom_text(hjust = 1, size = 3) +
  scale_x_continuous(limits = c(-0.1, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void()

# ============================
# Combine
# ============================

final_plot <- p_counts + p_heatmap +
  plot_layout(widths = c(2, 6))

#print(final_plot)

# --------------------------------------------------------------
# SAVE as high-resolution PDF (No changes needed here)
# --------------------------------------------------------------

ggsave(
  filename = output_pdf,
  plot = final_plot,
  width = 9,
  height = 10,
  units = "in",
  dpi = 800,
  device = cairo_pdf
)

cat("Saved PDF:", output_pdf, "\n")
