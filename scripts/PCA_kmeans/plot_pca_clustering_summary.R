#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
})

# --------------------------------------------------
# ARGUMENTS
# --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript multi_panel_clustering_plot.R <input_tsv> [output_png]")
}

input_file  <- args[1]
output_file <- ifelse(length(args) >= 2, args[2], "multi_panel_clustering.png")

if (!file.exists(input_file)) {
  stop("Input file not found.")
}

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

df <- read_tsv(input_file, show_col_types = FALSE)

# --------------------------------------------------
# AUTO-DETECT LENGTH COLUMN
# --------------------------------------------------

if ("Length_threshold" %in% colnames(df)) {
  length_col <- "Length_threshold"
} else if ("Length_bin" %in% colnames(df)) {
  length_col <- "Length_bin"
} else {
  stop("Input must contain either 'Length_threshold' or 'Length_bin'.")
}

required_cols <- c("Species", length_col,
                   "PC1_PC2", "k", "Silhouette", "DBI")

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# --------------------------------------------------
# CLEAN + ORDER FACTORS
# --------------------------------------------------

if (length_col == "Length_threshold") {
  df[[length_col]] <- factor(
    df[[length_col]],
    levels = c("gt2500", "gt5000", "gt7500", "gt10000"),
    ordered = TRUE
  )
} else {
  df[[length_col]] <- factor(
    df[[length_col]],
    levels = sort(unique(df[[length_col]])),
    ordered = TRUE
  )
}

df <- df %>%
  mutate(Species = str_to_title(Species))

# --------------------------------------------------
# RESHAPE FOR FACETING
# --------------------------------------------------

plot_df <- df %>%
  select(Species, all_of(length_col),
         PC1_PC2, Silhouette, DBI, k) %>%
  pivot_longer(
    cols = c(PC1_PC2, Silhouette, DBI, k),
    names_to = "Metric",
    values_to = "Value"
  )

plot_df$Metric <- recode(plot_df$Metric,
  PC1_PC2   = "PC1 + PC2 Variance",
  Silhouette = "Silhouette Score",
  DBI        = "Davies–Bouldin Index",
  k          = "Optimal k"
)

# Convert length factor to numeric index for smoothing
plot_df$Length_index <- as.numeric(plot_df[[length_col]])

# --------------------------------------------------
# COLOR PALETTE
# --------------------------------------------------

species_colors <- c(
  "Human" = "#D55E00",
  "Bonobo" = "#E69F00",
  "Borangutan" = "#009E73",
  "Sorangutan" = "#0072B2",
  "Chimpanzee" = "#CC79A7",
  "Gorilla" = "#56B4E9",
  "Gibbon" = "#000000"
)

# --------------------------------------------------
# PLOT
# --------------------------------------------------

p <- ggplot(plot_df,
            aes(x = Length_index,
                y = Value,
                group = Species,
                color = Species)) +

  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 2.5) +

  # --- NEW: smoothing curve ---
  geom_smooth(
    method = "loess",
    se = FALSE,
    span = 0.8,
    size = 1.2,
    alpha = 0.9
  ) +

  facet_wrap(~ Metric,
             scales = "free_y",
             ncol = 2) +

  scale_color_manual(values = species_colors) +

  scale_x_continuous(
    breaks = unique(plot_df$Length_index),
    labels = levels(df[[length_col]])
  ) +

  labs(
    title = paste("Comparative Clustering Metrics Across", length_col),
    x = gsub("_", " ", length_col),
    y = NULL,
    color = "Species"
  ) +

  theme_classic(base_size = 14) +

  theme(
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --------------------------------------------------
# SAVE OUTPUT
# --------------------------------------------------

ggsave(output_file,
       plot = p,
       width = 12,
       height = 9,
       dpi = 300)

cat("Multi-panel figure saved to:", output_file, "\n")

