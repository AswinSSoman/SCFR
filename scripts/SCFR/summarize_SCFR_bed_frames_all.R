#!/usr/bin/env Rscript

# Load required packages
suppressMessages(library(tools))

# Check for command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript summarize_bed_frames.R <input.bed>")
}

bed_file <- args[1]

# Check if file exists
if (!file.exists(bed_file)) {
  stop(paste("File not found:", bed_file))
}

# Read BED file
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bed_data) <- c("chrom", "start", "end", "frame")

# Clean and process
bed_data$start <- as.numeric(bed_data$start)
bed_data$end <- as.numeric(bed_data$end)
bed_data$length <- bed_data$end - bed_data$start
bed_data <- bed_data[complete.cases(bed_data[, c("chrom", "start", "end", "length", "frame")]), ]
bed_data <- bed_data[bed_data$length > 0, ]
frame_levels <- c("1", "2", "3", "-1", "-2", "-3")
bed_data$frame <- factor(as.character(bed_data$frame), levels = frame_levels)

# Function to get quantile label
get_quantile_label <- function(lengths, target_len) {
  q <- ecdf(lengths)
  prob <- q(target_len)
  return(paste0("Q", ceiling(prob * 100)))
}

# Function to compute summary stats
summarize_frame <- function(df, chrom_label, frame_label) {
  df_frame <- df[df$frame == frame_label & !is.na(df$length), ]
  lengths <- df_frame$length
  if (length(lengths) == 0) return(NULL)

  summary_stats <- data.frame(
    Chromosome = chrom_label,
    Frame = frame_label,
    N = length(lengths),
    Min = min(lengths),
    Q1 = quantile(lengths, 0.25),
    Median = median(lengths),
    Mean = mean(lengths),
    Q3 = quantile(lengths, 0.75),
    Max = max(lengths),
    SD = sd(lengths),
    P95 = quantile(lengths, 0.95),
    P99 = quantile(lengths, 0.99),
    Q_1Kb = get_quantile_label(lengths, 1000),
    Q_5Kb = get_quantile_label(lengths, 5000),
    Q_10Kb = get_quantile_label(lengths, 10000)
  )
  return(summary_stats)
}

# Summarize per chromosome and frame
chromosomes <- unique(bed_data$chrom)
summary_all <- list()

for (chrom in chromosomes) {
  chrom_data <- bed_data[bed_data$chrom == chrom, ]
  frame_summaries <- lapply(frame_levels, function(f) summarize_frame(chrom_data, chrom, f))
  frame_summaries <- do.call(rbind, frame_summaries)
  summary_all[[chrom]] <- frame_summaries
}

# Combine all per-chromosome summaries
summary_df <- do.call(rbind, summary_all)

# Write per-chromosome summary
prefix <- file_path_sans_ext(basename(bed_file))
out_file_chroms <- paste0(prefix, "_frame_summary_by_chromosome.csv")
write.csv(summary_df, file = out_file_chroms, row.names = FALSE)
cat("Per-chromosome summary written to:", out_file_chroms, "\n")

# Add global (all chromosomes) summary
combined_summary_list <- lapply(frame_levels, function(f) summarize_frame(bed_data, "ALL", f))
combined_summary_df <- do.call(rbind, combined_summary_list)

# Combine all into final file
summary_df_final <- rbind(summary_df, combined_summary_df)
out_file_final <- paste0(prefix, "_frame_summary_all.csv")
write.csv(summary_df_final, file = out_file_final, row.names = FALSE)
cat("Full summary (including ALL chromosomes) written to:", out_file_final, "\n")
