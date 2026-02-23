suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2){
  stop("Usage: Rscript pca_rscu_folder.R <input_folder> <output_folder>")
}

input_folder  <- args[1]
output_folder <- args[2]

rscu_file <- file.path(input_folder, "rscu_gc3_corrected_matrix.tsv")

if(!file.exists(rscu_file)){
  stop("File not found: ", rscu_file)
}

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

scores_out    <- file.path(output_folder, "pca_scores.tsv")
loadings_out  <- file.path(output_folder, "pca_loadings.tsv")
variance_out  <- file.path(output_folder, "explained_variance.tsv")

# ----------------------------
# Read data
# ----------------------------
rscu_df <- read_tsv(rscu_file, col_types = cols())

seq_ids  <- rscu_df[[1]]
rscu_mat <- rscu_df[,-1]

rownames(rscu_mat) <- seq_ids
rscu_mat <- as.matrix(rscu_mat)

# ----------------------------
# Scale + PCA
# ----------------------------
# ----------------------------
# Remove zero-variance columns
# ----------------------------
col_variances <- apply(rscu_mat, 2, var, na.rm = TRUE)

rscu_filtered <- rscu_mat[, col_variances > 0]

cat("Original codons:", ncol(rscu_mat), "\n")
cat("Codons retained for PCA:", ncol(rscu_filtered), "\n")

if (ncol(rscu_filtered) < 2) {
  stop("Not enough variable codons for PCA.")
}

# ----------------------------
# PCA (scale inside prcomp)
# ----------------------------
pca_res <- prcomp(rscu_filtered, center = TRUE, scale. = TRUE)
# ----------------------------
# Scores
# ----------------------------
scores_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("sequence")

write_tsv(scores_df, scores_out)

# ----------------------------
# Loadings
# ----------------------------
loadings_df <- as.data.frame(pca_res$rotation) %>%
  rownames_to_column("codon")

write_tsv(loadings_df, loadings_out)

# ----------------------------
# Explained Variance
# ----------------------------
eigenvalues <- pca_res$sdev^2
variance_ratio <- eigenvalues / sum(eigenvalues)
cumulative_variance <- cumsum(variance_ratio)

variance_df <- tibble(
  PC = paste0("PC", seq_along(eigenvalues)),
  eigenvalue = eigenvalues,
  variance_explained = variance_ratio,
  cumulative_variance = cumulative_variance
)

write_tsv(variance_df, variance_out)

cat("PCA completed.\n")
cat("Scores written to:", scores_out, "\n")
cat("Loadings written to:", loadings_out, "\n")
cat("Variance written to:", variance_out, "\n")
