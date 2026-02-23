# Load libraries
library(ggplot2)
library(ggrepel)
library(patchwork)

# Species list and display names
species_list <- c("human", "bonobo", "chimpanzee", "borangutan", "sorangutan", "gibbon")
labels <- LETTERS[1:6]
species_display <- c(
  "human" = "Human",
  "bonobo" = "Bonobo",
  "chimpanzee" = "Chimpanzee",
  "borangutan" = "Bornean Orangutan",
  "sorangutan" = "Sumatran Orangutan",
  "gibbon" = "Gibbon"
)

# Helper to extract lower bound of bin
extract_low <- function(x) as.numeric(sub("-.*", "", x))

# STEP 1: Build gene presence map
gene_sets <- list()
for (sp in species_list) {
  ann_file <- paste0(sp, "_genes_of_interest.txt")
  ann_raw <- read.table(ann_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("raw_bin2", "raw_bin1", "gene"))
  gene_sets[[sp]] <- unique(ann_raw$gene)
}
# Frequency of gene across species
gene_freq <- table(unlist(gene_sets))
gene_colors <- sapply(names(gene_freq), function(g) {
  if (gene_freq[g] == length(species_list)) {
    "green"
  } else if (gene_freq[g] >= 2) {
    "magenta"
  } else {
    "orange"
  }
})

# STEP 2: Global max count for shared color scale
max_count <- 0
for (sp in species_list) {
  awk_file <- paste0(sp, "_bins.out")
  awk_data <- read.table(awk_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         col.names = c("bin1", "bin2", "count"))
  max_count <- max(max_count, max(awk_data$count, na.rm = TRUE))
}
fill_scale <- scale_fill_gradient(low = "blue", high = "red", name = "Count", limits = c(0, max_count))

# STEP 3: Plot function for one species
create_plot <- function(species, label) {
  awk_file <- paste0(species, "_bins.out")
  ann_file <- paste0(species, "_genes_of_interest.txt")
  
  awk_data <- read.table(awk_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         col.names = c("bin1", "bin2", "count"))
  ann_raw <- read.table(ann_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("raw_bin2", "raw_bin1", "gene"))

  # Binning raw_bin1 and raw_bin2
  ann_raw$bin1 <- sapply(ann_raw$raw_bin1, function(x) {
    if (x < 1000) {
      low <- floor(x / 100) * 100
      high <- low + 99
    } else {
      low <- floor(x / 1000) * 1000
      high <- low + 999
    }
    sprintf("%d-%d", low, high)
  })
  ann_raw$bin2 <- sapply(ann_raw$raw_bin2, function(x) {
    low <- floor(x * 10) / 10
    high <- low + 0.1
    sprintf("%.1f-%.1f", low, high)
  })

  # Merge
  merged <- merge(awk_data, ann_raw[, c("bin1", "bin2", "gene")], by = c("bin1", "bin2"), all.x = TRUE)
  merged$gene[is.na(merged$gene)] <- ""

  # Assign label color
  merged$label_color <- gene_colors[merged$gene]
  merged$label_color[is.na(merged$label_color)] <- NA

  # Sort bins
  bin1_levels <- unique(merged$bin1)
  bin1_levels <- bin1_levels[order(extract_low(bin1_levels))]
  merged$bin1 <- factor(merged$bin1, levels = bin1_levels)

  bin2_levels <- unique(merged$bin2)
  bin2_levels <- bin2_levels[order(extract_low(bin2_levels))]
  merged$bin2 <- factor(merged$bin2, levels = bin2_levels)

  # Label data: allow multiple genes per bin, prefer non-LOC, no duplicates
  label_df <- subset(merged, !is.na(label_color) & gene != "")
  label_df$priority <- ifelse(grepl("^LOC", label_df$gene), 1, 0)
  label_df <- label_df[order(label_df$bin1, label_df$bin2, label_df$priority), ]
  label_df <- label_df[!duplicated(label_df$gene), ]

  # Plot
  species_label <- species_display[[species]]
  ggplot(merged, aes(x = bin2, y = bin1, fill = count)) +
    geom_tile(color = "white") +
    fill_scale +
    geom_text_repel(
      data = label_df,
      aes(label = gene, color = label_color),
      size = 4.5,
      box.padding = 0.6,
      point.padding = 0.4,
      segment.color = "grey50",
      max.overlaps = 500
    ) +
    scale_color_identity() +
    annotate("text", x = 0.5, y = length(bin1_levels) + 1, label = label,
             hjust = 0, vjust = 1, size = 6, fontface = "bold") +
    labs(
      title = species_label,
      x = "% AT Content",
      y = "SCFR length"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
}

# STEP 4: Generate plots
plots <- mapply(create_plot, species_list, labels, SIMPLIFY = FALSE)

# STEP 5: Combine plots with patchwork
combined_plot <- (plots[[1]] | plots[[2]] | plots[[3]]) /
                 (plots[[4]] | plots[[5]] | plots[[6]]) +
                 plot_layout(guides = "collect") & theme(legend.position = "right")

# STEP 6: Save to PNG
ggsave("combined_species_heatmaps.png", plot = combined_plot, width = 24, height = 16, dpi = 300)
