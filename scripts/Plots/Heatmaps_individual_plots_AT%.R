
# Load libraries
library(ggplot2)
library(ggrepel)
library(patchwork) 

# Define the fixed maximum SCFR bin to display for Gorilla
MAX_SCFR_BIN <- "33000-33999"

# Set working directory (Assuming this is required for file access)
setwd("/home/abhishek/SCFR/SCFR_all")

# Species list and display names (7 species)
species_list <- c("human", "bonobo", "chimpanzee", "borangutan", "sorangutan", "gibbon", "gorilla")
# Define 7 labels (A-G) for the 7 species
labels <- LETTERS[1:7] 
species_display <- c(
  "human" = "Human",
  "bonobo" = "Bonobo",
  "chimpanzee" = "Chimpanzee",
  "borangutan" = "Bornean Orangutan",
  "sorangutan" = "Sumatran Orangutan",
  "gibbon" = "Gibbon",
  "gorilla" = "Gorilla"
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

# Frequency of gene across species (now based on 7 species)
gene_freq <- table(unlist(gene_sets))
gene_colors <- sapply(names(gene_freq), function(g) {
  if (gene_freq[g] == length(species_list)) {
    "green" # Present in all 7 species
  } else if (gene_freq[g] >= 2) {
    "magenta" # Present in 2 to 6 species
  } else {
    "orange" # Present in 1 species
  }
})

# STEP 2: Global max count for shared color scale (Fill Color)
max_count <- 0
for (sp in species_list) {
  awk_file <- paste0(sp, "_bins.out")
  awk_data <- read.table(awk_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         col.names = c("bin1", "bin2", "count"))
  max_count <- max(max_count, max(awk_data$count, na.rm = TRUE))
}
fill_scale <- scale_fill_gradient(low = "blue", high = "red", name = "Count", limits = c(0, max_count))

# STEP 3: Plot function for one species - NOW INCLUDES SAVING
create_plot <- function(species, label) {
  awk_file <- paste0(species, "_bins.out")
  ann_file <- paste0(species, "_genes_of_interest.txt")
  
  awk_data <- read.table(awk_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                         col.names = c("bin1", "bin2", "count"))
  ann_raw <- read.table(ann_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                        col.names = c("raw_bin2", "raw_bin1", "gene"))
  
  # Binning raw_bin1 (SCFR length) and raw_bin2 (GC Content or AT Content, depending on data)
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
  
  # Merge data and annotations
  merged <- merge(awk_data, ann_raw[, c("bin1", "bin2", "gene")], by = c("bin1", "bin2"), all.x = TRUE)
  merged$gene[is.na(merged$gene)] <- ""
  
  # Assign label color based on gene frequency
  merged$label_color <- gene_colors[merged$gene]
  merged$label_color[is.na(merged$label_color)] <- NA
  
  
  # CONDITIONAL SCALING LOGIC: Limits Y-axis (bin1) ONLY for Gorilla
  if (species == "gorilla") {
    current_bin1_levels <- unique(merged$bin1)
    current_bin1_levels <- current_bin1_levels[order(extract_low(current_bin1_levels))]
    
    # Find the index of the max bin we want to keep
    max_index <- which(current_bin1_levels == MAX_SCFR_BIN)
    
    if (length(max_index) > 0) {
      # Trim levels to the specified max bin
      bin1_levels_to_use <- current_bin1_levels[1:max_index]
      
      # Filter the merged data to exclude bins beyond the limit
      merged <- subset(merged, merged$bin1 %in% bin1_levels_to_use)
      
    } else {
      # If the max bin isn't present, use all current levels
      bin1_levels_to_use <- current_bin1_levels
    }
    
  } else {
    # For all other species, use only the levels present in their data (original logic)
    bin1_levels_to_use <- unique(merged$bin1)
    bin1_levels_to_use <- bin1_levels_to_use[order(extract_low(bin1_levels_to_use))]
  }
  
  # Apply the determined factor levels to the data
  merged$bin1 <- factor(merged$bin1, levels = bin1_levels_to_use)
  
  # Sort bin2 (GC or AT Content) levels (always local)
  bin2_levels <- unique(merged$bin2)
  bin2_levels <- bin2_levels[order(extract_low(bin2_levels))]
  merged$bin2 <- factor(merged$bin2, levels = bin2_levels)
  
  # Prepare label data (prioritize non-LOC genes, remove duplicates)
  label_df <- subset(merged, !is.na(label_color) & gene != "")
  label_df$priority <- ifelse(grepl("^LOC", label_df$gene), 1, 0)
  label_df <- label_df[order(label_df$bin1, label_df$bin2, label_df$priority), ]
  label_df <- label_df[!duplicated(label_df$gene), ]
  
  # Determine annotation position based on the levels used
  y_annotation_pos <- length(bin1_levels_to_use) + 1
  
  # Plot
  species_label <- species_display[[species]]
  p <- ggplot(merged, aes(x = bin2, y = bin1, fill = count)) +
    geom_tile(color = "white") +
    fill_scale +
    geom_text_repel(
      # Filter label data one last time to match the filtered merged data
      data = subset(label_df, label_df$bin1 %in% bin1_levels_to_use), 
      aes(label = gene, color = label_color),
      size = 4.5,
      box.padding = 0.6,
      point.padding = 0.4,
      segment.color = "grey50",
      max.overlaps = 500
    ) +
    scale_color_identity() +
    # Use the determined y_annotation_pos
    annotate("text", x = 0.5, y = y_annotation_pos, label = label, 
             hjust = 0, vjust = 1, size = 6, fontface = "bold") +
    labs(
      title = species_label,
      # CORRECTED X-AXIS LABEL
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
  
  # SAVE PLOT HERE AS A SEPARATE PDF FILE
  filename <- paste0(species, "_heatmap.pdf")
  ggsave(filename, plot = p, width = 12, height = 8)
  
  # Return the plot object (optional, but good practice)
  return(p)
}

# STEP 4: Generate plots and SAVE INDIVIDUALLY
# mapply calls create_plot for each species, which saves the PDF inside the function.
plots <- mapply(create_plot, species_list, labels, SIMPLIFY = FALSE)

