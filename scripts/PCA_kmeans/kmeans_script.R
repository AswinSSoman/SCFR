#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(cluster)
  library(clusterSim)
})

# -----------------------------------
# ARGUMENTS
# -----------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript kmeans_script.R <input_folder> <output_folder>")
}

input_folder  <- args[1]
output_folder <- args[2]

pca_file <- file.path(input_folder, "pca_scores.tsv")

if (!file.exists(pca_file)) {
  stop("pca_scores.tsv not found in input folder.")
}

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# -----------------------------------
# LOAD DATA
# -----------------------------------

pca_scores <- read_tsv(pca_file, show_col_types = FALSE)

if (!"sequence" %in% colnames(pca_scores)) {
  stop("Column 'sequence' not found.")
}

pc_cols <- grep("^PC[0-9]+$", colnames(pca_scores), value = TRUE)

if (length(pc_cols) < 2) {
  stop("Less than 2 principal components found.")
}

cat("Detected", length(pc_cols), "principal components.\n")

X <- as.matrix(pca_scores[, pc_cols])

# -----------------------------------
# K RANGE
# -----------------------------------

k_values <- 2:10

sil_values  <- numeric(length(k_values))
dbi_values  <- numeric(length(k_values))
wcss_values <- numeric(length(k_values))

dmat <- dist(X)

# -----------------------------------
# OPTIMIZATION LOOP
# -----------------------------------

for (i in seq_along(k_values)) {

  k <- k_values[i]

  set.seed(42)
  km <- kmeans(X, centers = k, nstart = 50)

  wcss_values[i] <- km$tot.withinss

  sil <- silhouette(km$cluster, dmat)
  sil_values[i] <- mean(sil[, 3])

  dbi_values[i] <- clusterSim::index.DB(X, km$cluster)$DB
}

# -----------------------------------
# ELBOW DETECTION (CURVATURE METHOD)
# -----------------------------------

curvature <- rep(NA, length(wcss_values))

for (i in 2:(length(wcss_values) - 1)) {
  curvature[i] <- wcss_values[i - 1] -
                  2 * wcss_values[i] +
                  wcss_values[i + 1]
}

elbow_index <- which.max(curvature)
best_k <- k_values[elbow_index]

# -----------------------------------
# SAVE METRICS TSV
# -----------------------------------

metrics_df <- data.frame(
  k = k_values,
  Silhouette = sil_values,
  DBI = dbi_values,
  WCSS = wcss_values,
  Curvature = curvature
)

write_tsv(metrics_df,
          file.path(output_folder, "k_optimization_scores.tsv"))

# -----------------------------------
# FINAL KMEANS
# -----------------------------------

set.seed(42)
final_km <- kmeans(X, centers = best_k, nstart = 100)

cluster_assignments <- data.frame(
  sequence = pca_scores$sequence,
  Cluster = final_km$cluster
)

write_tsv(cluster_assignments,
          file.path(output_folder, "cluster_assignments.tsv"))

# -----------------------------------
# CLUSTER TXT FILES
# -----------------------------------

clusters <- sort(unique(final_km$cluster))

for (cl in clusters) {

  seq_ids <- cluster_assignments$sequence[
    cluster_assignments$Cluster == cl
  ]

  writeLines(seq_ids,
             file.path(output_folder,
                       paste0("cluster_", cl, "_sequences.txt")))
}

# -----------------------------------
# CLUSTER SUMMARY TSV
# -----------------------------------

summary_df <- cluster_assignments %>%
  group_by(Cluster) %>%
  summarise(
    n_sequences = n(),
    sequences = paste(sequence, collapse = ",")
  )

write_tsv(summary_df,
          file.path(output_folder, "cluster_summary.tsv"))

cat("Best k selected (Elbow / Curvature):", best_k, "\n")
cat("KMeans clustering complete.\n")
