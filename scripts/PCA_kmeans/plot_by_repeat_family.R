#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(grid)
  library(ggrepel)
})

# -----------------------------------
# ARGUMENTS
# -----------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  stop("Usage: Rscript plot_pca_repeat_biplot.R <input_folder> <output_folder> <repeat.tsv> <species_name> <plot_width> <plot_height> <point_size> <length>")
}

input_folder  <- args[1]
output_folder <- args[2]
repeat_file   <- args[3]
species_name <- args[4]
pw <- as.numeric(args[5])
ph <- as.numeric(args[6])
ps <- as.numeric(args[7])
ln <- args[8]

scores_file    <- file.path(input_folder, "pca_scores.tsv")
cluster_file   <- file.path(output_folder, "cluster_assignments.tsv")
loadings_file  <- file.path(input_folder, "pca_loadings.tsv")
variance_file  <- file.path(input_folder, "explained_variance.tsv")

# -----------------------------------
# CHECK FILES
# -----------------------------------

for (f in c(scores_file, cluster_file, loadings_file, variance_file, repeat_file)) {
  if (!file.exists(f)) stop(paste("Missing file:", f))
}

# -----------------------------------
# LOAD DATA
# -----------------------------------

scores   <- read_tsv(scores_file, show_col_types = FALSE)
clusters <- read_tsv(cluster_file, show_col_types = FALSE)
loadings <- read_tsv(loadings_file, show_col_types = FALSE)
variance <- read_tsv(variance_file, show_col_types = FALSE)
repeats  <- read_tsv(repeat_file, show_col_types = FALSE)

# -----------------------------------
# CLEAN REPEAT DATA
# -----------------------------------

repeats <- repeats %>%
  select(query_sequence, repeat_class_family)

# -----------------------------------
# MERGE DATA
# -----------------------------------

plot_df <- scores %>%
  inner_join(clusters, by = "sequence") %>%
  left_join(repeats, by = c("sequence" = "query_sequence"))

plot_df$repeat_class_family[is.na(plot_df$repeat_class_family)] <- "No_repeat"

# -----------------------------------
# PCs
# -----------------------------------

pc_cols  <- paste0("PC", 1:5)
pc_pairs <- combn(pc_cols, 2, simplify = FALSE)

arrow_scale  <- 5
label_offset <- 2

# -----------------------------------
# LOOP OVER PC PAIRS
# -----------------------------------

for (pair in pc_pairs) {

  x_pc <- pair[1]
  y_pc <- pair[2]

  var_x <- variance$variance_explained[variance$PC == x_pc]
  var_y <- variance$variance_explained[variance$PC == y_pc]

  top_loadings <- loadings %>%
    select(codon, all_of(x_pc), all_of(y_pc)) %>%
    mutate(magnitude = abs(.data[[x_pc]]) + abs(.data[[y_pc]])) %>%
    arrange(desc(magnitude)) %>%
    slice(1:5)

  pdf_file <- file.path(
    output_folder,
    paste0("pca_repeat_", x_pc, "_", y_pc, ".pdf")
  )

  pdf(pdf_file, width = pw, height = ph)

  p <- ggplot(plot_df,
              aes(x = !!sym(x_pc),
                  y = !!sym(y_pc),
                  color = repeat_class_family)) +

    geom_point(size = ps, alpha = 0.8) +

    geom_segment(
      data = top_loadings,
      aes(x = 0, y = 0,
          xend = !!sym(x_pc) * arrow_scale,
          yend = !!sym(y_pc) * arrow_scale),
      arrow = arrow(length = unit(0.2, "cm")),
      inherit.aes = FALSE
    ) +

geom_text_repel(
  data = top_loadings,
  aes(x = !!sym(x_pc) * arrow_scale,
      y = !!sym(y_pc) * arrow_scale,
      label = codon),
  size = 4,
  fontface = "bold",
  inherit.aes = FALSE,
  box.padding = 0.5,
  point.padding = 0.3,
  segment.color = "grey40",
  max.overlaps = Inf
) +

    labs(
      title = paste(species_name, "-", ln, "-", x_pc, "vs", y_pc),
      x = paste0(x_pc, " (", round(var_x, 4)*100, "%)"),
      y = paste0(y_pc, " (", round(var_y, 4)*100, "%)"),
      color = "Repeat Class"
    ) +

    theme_minimal() +
theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16)
)

  print(p)
  dev.off()
}

cat("All PCA plots saved colored by repeat_class_family.\n")
