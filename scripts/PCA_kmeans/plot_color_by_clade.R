#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(rlang)
  library(ggrepel)
})

# -----------------------------------
# ARGUMENTS
# -----------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript plot_pca_biplot.R <input_folder> <output_folder> <species_name> <plot_width> <plot_height> <point_size> <length>")
}

input_folder  <- args[1]
output_folder <- args[2]
species_name <- args[3]
pw <- as.numeric(args[4])
ph <- as.numeric(args[5])
ps <- as.numeric(args[6])
ln <- args[7]

scores_file    <- file.path(input_folder, "pca_scores.tsv")
loadings_file  <- file.path(input_folder, "pca_loadings.tsv")
variance_file  <- file.path(input_folder, "explained_variance.tsv")
cluster_file   <- file.path(output_folder, "cluster_assignments.tsv")

# -----------------------------------
# CHECK FILES
# -----------------------------------

for (f in c(scores_file, loadings_file, variance_file, cluster_file)) {
  if (!file.exists(f)) stop(paste("Missing file:", f))
}

# -----------------------------------
# LOAD DATA
# -----------------------------------

scores   <- read_tsv(scores_file, show_col_types = FALSE)
loadings <- read_tsv(loadings_file, show_col_types = FALSE)
variance <- read_tsv(variance_file, show_col_types = FALSE)
clusters <- read_tsv(cluster_file, show_col_types = FALSE)

# -----------------------------------
# MERGE + CLADE
# -----------------------------------

plot_df <- scores %>%
  inner_join(clusters, by = "sequence") %>%
  mutate(Clade = sub(".*::([^:]+):.*", "\\1", sequence))

# -----------------------------------
# PCs
# -----------------------------------

pc_cols  <- paste0("PC", 1:5)
pc_pairs <- combn(pc_cols, 2, simplify = FALSE)
#arrow_scale <- 1
#label_offset <- 0.5

# -----------------------------------
# LOOP OVER PC PAIRS - SEPARATE PDFs
# -----------------------------------

for (pair in pc_pairs) {

  x_pc <- pair[1]
  y_pc <- pair[2]

  # ---- Variance explained for axis labels ----
  var_x <- variance$variance_explained[variance$PC == x_pc]
  var_y <- variance$variance_explained[variance$PC == y_pc]

  # ---- Top 5 loadings ----
  top_loadings <- loadings %>%
    select(codon, all_of(x_pc), all_of(y_pc)) %>%
    mutate(magnitude = abs(.data[[x_pc]]) + abs(.data[[y_pc]])) %>%
    arrange(desc(magnitude)) %>%
    slice(1:5)

#Dynamic labelling
max_score <- max(abs(plot_df[[x_pc]]), abs(plot_df[[y_pc]]))
max_loading <- max(abs(top_loadings[[x_pc]]),
                     abs(top_loadings[[y_pc]]))

arrow_scale <- 0.1 * max_score / max_loading
label_offset <- 0.1 * arrow_scale

  # ---- PDF FILE ----
  pdf_file <- file.path(output_folder,
                        paste0("pca_clade_", x_pc, "_", y_pc, ".pdf"))

  pdf(pdf_file, width = pw, height = ph)

  # ---- PLOT ----
  p <- ggplot(plot_df, aes(x = !!sym(x_pc), y = !!sym(y_pc), color = Clade)) +
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
      color = "Chromsome"
    ) +

theme_minimal() +

guides(
  color = guide_legend(ncol = 1)   # ← single column
) +

theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 16),
  axis.text  = element_text(size = 16),

  legend.position = "right",

  legend.title = element_text(size = 7, face = "bold"),
  legend.text  = element_text(size = 7),

  # 🔥 Compact legend spacing
  legend.key.height = unit(0.35, "cm"),
  legend.key.width  = unit(0.35, "cm"),
  legend.spacing.y  = unit(0.05, "cm"),
  legend.margin     = margin(2, 2, 2, 2)
)

  print(p)
  dev.off()
}

cat("PCA biplots with top 5 codon loadings and variance annotated saved as separate PDFs (pca_clade_...).\n")
