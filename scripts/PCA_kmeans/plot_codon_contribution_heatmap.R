library(pheatmap)
library(dplyr)
library(grid)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 7){
  stop("Usage: Rscript script.R <loadings> <output> <species> <len> <n_pcs> <top_n> <variance_file>")
}

input_file    <- args[1]
output_file   <- args[2]
sp            <- args[3]
ln            <- args[4]
n_pcs         <- as.numeric(args[5])
top_n         <- as.numeric(args[6])
variance_file <- args[7]

############################################################
# Load Data
############################################################

loadings_table <- read.table(input_file, header=TRUE)
variance_table <- read.table(variance_file, header=TRUE)

codon_names <- loadings_table$codon
loadings <- as.matrix(loadings_table[,-1])

pc_names <- colnames(loadings)[1:n_pcs]
loadings <- loadings[, pc_names, drop=FALSE]

# Keep raw loadings separately for magnitude ranking
raw_loadings <- loadings

############################################################
# Get Variance Explained
############################################################

variance_lookup <- variance_table %>%
  filter(PC %in% pc_names)

variance_percent <- round(variance_lookup$variance_explained * 100, 2)
names(variance_percent) <- variance_lookup$PC

cumulative_selected <- round(
  variance_lookup$cumulative_variance[n_pcs] * 100,
  2
)

############################################################
# Signed % Contribution (for heatmap display)
############################################################

sq_loadings <- loadings^2
contrib <- sweep(sq_loadings, 2, colSums(sq_loadings), FUN="/") * 100
signed_contrib <- contrib * sign(loadings)

contrib_df <- as.data.frame(signed_contrib)
contrib_df$codon <- codon_names

############################################################
# Add Combined PC1 + PC2 (display only)
############################################################

if(n_pcs >= 2){
  combined_name <- paste0(pc_names[1], "+", pc_names[2])
  contrib_df[[combined_name]] <-
    contrib_df[[pc_names[1]]] +
    contrib_df[[pc_names[2]]]
  pc_extended <- c(pc_names, combined_name)
} else {
  pc_extended <- pc_names
}

############################################################
# Select Top Codons
############################################################

top_list <- list()

# ---- Individual PCs (unchanged logic) ----
for(pc in pc_names){
  
  top_codons <- raw_loadings %>%
    as.data.frame() %>%
    mutate(codon = codon_names,
           abs_val = abs(.data[[pc]])) %>%
    arrange(desc(abs_val)) %>%
    slice(1:top_n) %>%
    pull(codon)
  
  top_list[[pc]] <- top_codons
}

# ---- PC1+PC2 selection FIXED ----
if(n_pcs >= 2){
  
  combined_metric <-
    abs(raw_loadings[, pc_names[1]]) +
    abs(raw_loadings[, pc_names[2]])
  
  combined_top <-
    codon_names[order(combined_metric, decreasing = TRUE)][1:top_n]
  
  top_list[[combined_name]] <- combined_top
}

# Final union of selected codons
selected_codons <- unique(unlist(top_list))

############################################################
# (Rest of your code unchanged)
############################################################
############################################################
# Genetic Code
############################################################

genetic_code <- c(
  TTT="Phe", TTC="Phe",
  TTA="Leu", TTG="Leu", CTT="Leu", CTC="Leu", CTA="Leu", CTG="Leu",
  ATT="Ile", ATC="Ile", ATA="Ile",
  ATG="Met",
  GTT="Val", GTC="Val", GTA="Val", GTG="Val",
  TCT="Ser", TCC="Ser", TCA="Ser", TCG="Ser", AGT="Ser", AGC="Ser",
  CCT="Pro", CCC="Pro", CCA="Pro", CCG="Pro",
  ACT="Thr", ACC="Thr", ACA="Thr", ACG="Thr",
  GCT="Ala", GCC="Ala", GCA="Ala", GCG="Ala",
  TAT="Tyr", TAC="Tyr",
  CAT="His", CAC="His",
  CAA="Gln", CAG="Gln",
  AAT="Asn", AAC="Asn",
  AAA="Lys", AAG="Lys",
  GAT="Asp", GAC="Asp",
  GAA="Glu", GAG="Glu",
  TGT="Cys", TGC="Cys",
  TGG="Trp",
  CGT="Arg", CGC="Arg", CGA="Arg", CGG="Arg", AGA="Arg", AGG="Arg",
  GGT="Gly", GGC="Gly", GGA="Gly", GGG="Gly"
)

############################################################
# Amino Acid Properties
############################################################

aa_properties <- read.csv(text="
AminoAcid,Hydropathy,Volume,Chemical,Charge,Polarity,Hbond
Ala,hydrophobic,verysmall,aliphatic,uncharged,nonpolar,none
Arg,hydrophilic,large,basic,positive,polar,donor
Asn,hydrophilic,small,amide,uncharged,polar,donor_acceptor
Asp,hydrophilic,small,acidic,negative,polar,acceptor
Cys,hydrophobic,small,sulfur,uncharged,nonpolar,none
Gln,hydrophilic,medium,amide,uncharged,polar,donor_acceptor
Glu,hydrophilic,medium,acidic,negative,polar,acceptor
Gly,neutral,verysmall,aliphatic,uncharged,nonpolar,none
His,neutral,medium,basic,positive,polar,donor_acceptor
Ile,hydrophobic,large,aliphatic,uncharged,nonpolar,none
Leu,hydrophobic,large,aliphatic,uncharged,nonpolar,none
Lys,hydrophilic,large,basic,positive,polar,donor
Met,hydrophobic,large,sulfur,uncharged,nonpolar,none
Phe,hydrophobic,verylarge,aromatic,uncharged,nonpolar,none
Pro,neutral,small,aliphatic,uncharged,nonpolar,none
Ser,neutral,verysmall,hydroxyl,uncharged,polar,donor_acceptor
Thr,neutral,small,hydroxyl,uncharged,polar,donor_acceptor
Trp,hydrophobic,verylarge,aromatic,uncharged,nonpolar,donor
Tyr,neutral,verylarge,aromatic,uncharged,polar,donor_acceptor
Val,hydrophobic,medium,aliphatic,uncharged,nonpolar,none
")

row.names(aa_properties) <- aa_properties$AminoAcid

############################################################
# Prepare Heatmap Data
############################################################

heatmap_data <- contrib_df %>%
  filter(codon %in% selected_codons)

# Add amino acid in brackets to row labels
aa_labels <- genetic_code[heatmap_data$codon]

row.names(heatmap_data) <- paste0(
  heatmap_data$codon,
  " (",
  aa_labels,
  ")"
)

heatmap_data <- heatmap_data[, pc_extended, drop=FALSE]

############################################################
# Annotation
############################################################

get_gc3 <- function(codon){
  third <- substr(codon,3,3)
  if(third %in% c("G","C")) return("GC3") else return("AT3")
}

#aa_vector  <- genetic_code[row.names(heatmap_data)]
codon_only <- sub(" \\(.*\\)", "", row.names(heatmap_data))
aa_vector  <- genetic_code[codon_only]
#gc3_vector <- sapply(row.names(heatmap_data), get_gc3)
gc3_vector <- sapply(codon_only, get_gc3)

annotation_row <- aa_properties[aa_vector, -1]
annotation_row$GC3 <- gc3_vector
rownames(annotation_row) <- rownames(heatmap_data)

############################################################
# Cleaner & More Digestible Annotation Colors
############################################################

ann_colors <- list(

  Hydropathy = c(
    hydrophobic = "#D55E00",
    hydrophilic = "#0072B2",
    neutral     = "#999999"
  ),

  Charge = c(
    positive  = "#E41A1C",
    negative  = "#377EB8",
    uncharged = "#999999"
  ),

  Polarity = c(
    polar    = "#1B9E77",
    nonpolar = "#E69F00"
  ),

  Volume = c(
    verysmall = "#F0F0F0",
    small     = "#D9D9D9",
    medium    = "#BDBDBD",
    large     = "#969696",
    verylarge = "#636363"
  ),

  Chemical = c(
    aliphatic = "#4DAF4A",
    aromatic  = "#984EA3",
    basic     = "#E41A1C",
    acidic    = "#377EB8",
    amide     = "#FF7F00",
    sulfur    = "#A65628",
    hydroxyl  = "#17BECF"
  ),

  Hbond = c(
    none           = "#CCCCCC",
    donor          = "#66C2A5",
    acceptor       = "#FC8D62",
    donor_acceptor = "#8DA0CB"
  ),

  GC3 = c(
    GC3 = "#2166AC",
    AT3 = "#B2182B"
  )
)

############################################################
# Column Labels (Cumulative in 3rd column)
############################################################

if(n_pcs >= 2){
  col_labels <- c(
    paste0(pc_names[1], "\n(", variance_percent[pc_names[1]], "%)"),
    paste0(pc_names[2], "\n(", variance_percent[pc_names[2]], "%)"),
    paste0(pc_names[1], "+", pc_names[2],
           "\n(", cumulative_selected, "%)")
  )
} else {
  col_labels <- paste0(
    pc_names,
    "\n(", variance_percent[pc_names], "%)"
  )
}

colnames(heatmap_data) <- col_labels

############################################################
# Plot
############################################################

max_abs <- max(abs(heatmap_data))
breaks_seq <- seq(-max_abs, max_abs, length.out = 101)

# More legend ticks (7 evenly spaced)
legend_ticks <- seq(-max_abs, max_abs, length.out = 7)

pdf(output_file, width=6, height=6)

pc_range_label <- paste0("PC1 - PC", n_pcs)

main_title <- sprintf(
  "PCA codon relative contribution (%s) | %s - SCFR >= %s bp",
  pc_range_label, sp, ln
)

pheatmap(
  heatmap_data,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  annotation_row    = annotation_row,
annotation_colors = ann_colors,
  annotation_names_row = TRUE,
  annotation_legend = TRUE,
  annotation_names_col = FALSE,
  display_numbers   = TRUE,
  number_format     = "%.1f",
  fontsize          = 8,
  fontsize_row      = 9,
  fontsize_col      = 9,
  fontsize_number   = 9,
  angle_col         = 90,
  legend            = TRUE,
  legend_breaks     = legend_ticks,
  legend_labels     = round(legend_ticks, 1),
  color = colorRampPalette(
    c("#6BAED6", "#BDD7E7", "white", "#FCBBA1", "#EF3B2C")
  )(100),
  breaks = breaks_seq,
  main   = main_title
)

dev.off()
