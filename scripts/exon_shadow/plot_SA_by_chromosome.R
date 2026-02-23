#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(zoo)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_SA_by_chromosome.R <input.tsv> <output_prefix>")
}

input_file <- args[1]
out_prefix <- args[2]

# ----------------------------
# Parameters
# ----------------------------
moderate_cutoff <- 0.4
extreme_cutoff  <- 0.8
min_counts      <- 5
smooth_k        <- 5

options(scipen = 999)

# ----------------------------
# Read data
# ----------------------------
df <- read_table(
  input_file,
  col_names = c("chr","start","end","plus","minus","SA"),
  show_col_types = FALSE
) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end),
    total = plus + minus
  ) %>%
  filter(total >= min_counts)

if (nrow(df) == 0) stop("No windows left after filtering")

# ----------------------------
# Bias classification
# ----------------------------
df <- df %>%
  mutate(
    bias_class = case_when(
      SA >=  extreme_cutoff ~ "Extreme +",
      SA <= -extreme_cutoff ~ "Extreme -",
      SA >=  moderate_cutoff ~ "Moderate +",
      SA <= -moderate_cutoff ~ "Moderate -",
      TRUE                  ~ "Neutral"
    ),
    is_hotspot = abs(SA) >= extreme_cutoff
  )

df$bias_class <- factor(
  df$bias_class,
  levels = c("Extreme +","Moderate +","Neutral","Moderate -","Extreme -")
)

# ----------------------------
# Chromosome summary
# ----------------------------
chrom_summary <- df %>%
  group_by(chr) %>%
  summarise(
    chr_len = max(end),
    mean_SA = mean(SA),
    median_SA = median(SA),
    n_windows = n(),
    n_extreme = sum(is_hotspot),
    dominant_bias = ifelse(mean_SA > 0, "-", "+"),
    unit = ifelse(chr_len >= 1e7, "Mb", "kb"),
    scale = ifelse(unit == "Mb", 1e6, 1e3),
    .groups = "drop"
  )

# ----------------------------
# Prepare plotting data
# ----------------------------
df_plot <- df %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(
    mid = (start + end) / 2,
    SA_smooth = rollmean(SA, k = smooth_k, fill = NA, align = "center")
  ) %>%
  ungroup() %>%
  left_join(chrom_summary, by = "chr") %>%
  mutate(x_scaled = mid / scale)

# ----------------------------
# Custom facet labels (single line, with units)
# ----------------------------
facet_labels <- chrom_summary %>%
  mutate(
    label = paste0(
      chr,
      " | mean=", round(mean_SA, 3),
      " median=", round(median_SA, 3),
      " windows=", n_windows,
      " extreme=", n_extreme,
      " bias=", dominant_bias,
      " | x: ", unit
    )
  ) %>%
  select(chr, label)

facet_labeller <- setNames(facet_labels$label, facet_labels$chr)

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(df_plot, aes(x = x_scaled, y = SA)) +

  geom_point(aes(color = bias_class), size = 0.6, alpha = 0.7) +

  # Highlight hotspot windows
#  geom_point(
#    data = subset(df_plot, is_hotspot),
#    shape = 21,
#    size = 1.4,
#    stroke = 0,
#    fill = NA,
#    color = "black"
 # ) +

  geom_line(
    aes(y = SA_smooth),
    linewidth = 0.25,
    linetype = "dashed",
    color = "black"
  ) +

  facet_wrap(
    ~ chr,
    scales = "free_x",
    ncol = 1,
    labeller = labeller(chr = facet_labeller)
  ) +

  scale_x_continuous(breaks = pretty_breaks(n = 6)) +

  scale_color_manual(
    values = c(
      "Extreme +"  = "#b2182b",
      "Moderate +" = "#ef8a62",
      "Neutral"    = "grey70",
      "Moderate -" = "#67a9cf",
      "Extreme -"  = "#2166ac"
    )
  ) +

  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-moderate_cutoff, moderate_cutoff),
             linetype = "dotted") +
  geom_hline(yintercept = c(-extreme_cutoff, extreme_cutoff),
             linetype = "solid") +

  labs(
    x = "Genomic position",
    y = "Strand Asymmetry (SA)",
    color = "Bias class"
  ) +

  theme_bw() +
  theme(
    strip.text = element_text(size = 8, hjust = 0, face="bold"),
    legend.position = "top"
  )

#Change chro size units
#format_genomic_axis <- function(x) {
#  max_x <- max(x, na.rm = TRUE)
#  if (max_x >= 1e6) {
#    scale_x_continuous(
#      labels = function(v) paste0(v / 1e6, " Mb"),
#      breaks = scales::pretty_breaks(n = 6)
#    )
#  } else {
#    scale_x_continuous(
#      labels = function(v) paste0(v / 1e3, " kb"),
#      breaks = scales::pretty_breaks(n = 6)
#    )
#  }
#}

#p1 <- p + format_genomic_axis(df$end)

ggsave(
  paste0(out_prefix, "_SA_by_chromosome.pdf"),
  p,
  width = 14,
  height = max(4, length(unique(df$chr)) * 1.5)
)

