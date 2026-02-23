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
  stop("Usage: Rscript plot_pca_gc_biplot.R <input_folder> <output_folder> <gc_file> <species_name> <plot_width> <plot_height> <point_size> <length>")
}

input_folder  <- args[1]
output_folder <- args[2]
gc_file       <- args[3]
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

for (f in c(scores_file, cluster_file, loadings_file, variance_file, gc_file)) {
  if (!file.exists(f)) stop(paste("Missing file:", f))
}

# -----------------------------------
# LOAD DATA
# -----------------------------------

scores   <- read_tsv(scores_file, show_col_types = FALSE)
clusters <- read_tsv(cluster_file, show_col_types = FALSE)
loadings <- read_tsv(loadings_file, show_col_types = FALSE)
variance <- read_tsv(variance_file, show_col_types = FALSE)
gc_data  <- read_tsv(gc_file, show_col_types = FALSE)

# -----------------------------------
# CLEAN GC DATA
# -----------------------------------

gc_data <- gc_data %>%
  select(seqname, GC)

# Auto-detect scale
if (max(gc_data$GC, na.rm = TRUE) <= 1) {
  gc_data <- gc_data %>%
    mutate(GC = GC * 100)
}

# -----------------------------------
# MERGE DATA
# -----------------------------------

plot_df <- scores %>%
  inner_join(clusters, by = "sequence") %>%
  left_join(gc_data, by = c("sequence" = "seqname"))

# Remove rows without GC
plot_df <- plot_df %>%
  filter(!is.na(GC))

# -----------------------------------
# Compute cluster centroids + mean GC
# -----------------------------------

cluster_summary <- plot_df %>%
  group_by(Cluster) %>%
  summarise(
    across(starts_with("PC"), mean, .names = "mean_{.col}"),
    mean_GC = mean(GC),
    .groups = "drop"
  )

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

x_range <- range(plot_df[[x_pc]])
y_range <- range(plot_df[[y_pc]])

x_min <- x_range[1]
x_max <- x_range[2]
y_min <- y_range[1]
y_max <- y_range[2]

# small padding
pad_x <- 0.1 * diff(x_range)
pad_y <- 0.1 * diff(y_range)

  var_x <- variance$variance_explained[variance$PC == x_pc]
  var_y <- variance$variance_explained[variance$PC == y_pc]

  top_loadings <- loadings %>%
    select(codon, all_of(x_pc), all_of(y_pc)) %>%
    mutate(magnitude = abs(.data[[x_pc]]) + abs(.data[[y_pc]])) %>%
    arrange(desc(magnitude)) %>%
    slice(1:5)

cluster_labels <- cluster_summary %>%
  mutate(
    x_centroid = .data[[paste0("mean_", x_pc)]],
    y_centroid = .data[[paste0("mean_", y_pc)]],
    label = paste0("C",Cluster, ": ",
               round(mean_GC, 2), "%")
  ) %>%
  mutate(
    # assign quadrant-based corner
    x_label = ifelse(x_centroid >= 0, x_max + pad_x, x_min - pad_x),
    y_label = ifelse(y_centroid >= 0, y_max + pad_y, y_min - pad_y)
  )

  pdf_file <- file.path(
    output_folder,
    paste0("pca_gc_", x_pc, "_", y_pc, ".pdf")
  )

  pdf(pdf_file, width = pw, height = ph)

  p <- ggplot(plot_df,
              aes(x = !!sym(x_pc),
                  y = !!sym(y_pc),
                  color = GC)) +

    geom_point(size = ps, alpha = 0.9) +

    # Continuous high-resolution gradient
    scale_color_gradientn(
      colours = c("navy", "blue", "cyan", "yellow", "orange", "red"),
      values = scales::rescale(c(
        min(plot_df$GC),
        45,
        50,
        55,
        60,
        max(plot_df$GC)
      )),
      name = "GC (%)"
    ) +

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

geom_segment(
  data = cluster_labels,
  aes(x = x_centroid,
      y = y_centroid,
      xend = x_label,
      yend = y_label),
  inherit.aes = FALSE,
  linewidth = 0.5,
  linetype = "dashed",
  arrow = arrow(length = unit(0.25, "cm"))
) +

geom_label(
  data = cluster_labels,
  aes(x = x_label, y = y_label, label = label),
  inherit.aes = FALSE,
  size = 4.5,
  fontface = "bold",
  fill = "white",
  label.size = 0.4
) +

    labs(
      title = paste(species_name, "-", ln, "-", x_pc, "vs", y_pc),
      x = paste0(x_pc, " (", round(var_x, 4)*100, "%)"),
      y = paste0(y_pc, " (", round(var_y, 4)*100, "%)")
    ) +

    theme_minimal() +

    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )

  print(p)
  dev.off()
}

cat("All PCA plots saved colored by GC content.\n")

