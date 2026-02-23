#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(optparse)
})

# ---------------------------
# Command-line arguments
# ---------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="CSV with strand asymmetry values"),
  make_option(c("-p", "--pdf"), type="character", default="asymmetry_plots.pdf", help="PDF file for plots"),
  make_option(c("-b", "--bed"), type="character", default="extreme_regions.bed", help="BED file for regions"),
  make_option(c("--min_region_size"), type="integer", default=100000, help="Min region size in bp"),
  make_option(c("--window_len"), type="integer", default=5, help="Sliding window size (#windows) for enrichment detection"),
  make_option(c("--min_hits"), type="integer", default=3, help="Min outliers in window to call region")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---------------------------
# Load CSV
# ---------------------------
df <- read_csv(opt$input, show_col_types = FALSE)

# ---------------------------
# Detect outliers per chromosome
# ---------------------------
df <- df %>%
  group_by(chrom) %>%
  mutate(
    Q1 = quantile(strand_count_asym, 0.25, na.rm = TRUE),
    Q3 = quantile(strand_count_asym, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    high_thresh = Q3 + 1.5 * IQR,
    low_thresh = Q1 - 1.5 * IQR,
    region_type = case_when(
      strand_count_asym >= high_thresh ~ "high",
      strand_count_asym <= low_thresh ~ "low",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# ---------------------------
# Sliding window detection of enriched extreme regions
# ---------------------------
detect_regions <- function(data, window_len = 5, min_hits = 3, region_type = "high") {
  indices <- which(data$region_type == region_type)
  if (length(indices) < min_hits) return(tibble())

  hits <- rep(FALSE, nrow(data))

  for (i in 1:(nrow(data) - window_len + 1)) {
    window_region_types <- data$region_type[i:(i + window_len - 1)]
    if (sum(window_region_types == region_type, na.rm = TRUE) >= min_hits) {
      hits[i:(i + window_len - 1)] <- TRUE
    }
  }

  # Extract contiguous hit regions
  data$hit_flag <- hits
  regions <- data %>%
    filter(hit_flag) %>%
    mutate(grp = cumsum(c(TRUE, diff(window_start) > 1.5 * median(window_end - window_start)))) %>%
    group_by(chrom, grp) %>%
    summarise(
      region_type = region_type[which.max(!is.na(region_type))],
      start = min(window_start),
      end = max(window_end),
      .groups = "drop"
    ) %>%
    filter(end - start >= opt$min_region_size) %>%
    select(chrom, start, end, region_type)

  return(regions)
}

all_regions <- bind_rows(
  lapply(split(df, df$chrom), function(chrom_df) {
    bind_rows(
      detect_regions(chrom_df, opt$window_len, opt$min_hits, "high"),
      detect_regions(chrom_df, opt$window_len, opt$min_hits, "low")
    )
  })
)

write_tsv(all_regions, opt$bed, col_names = FALSE)
cat("Extreme strand asymmetry regions written to:", opt$bed, "\n")

# ---------------------------
# Generate plots
# ---------------------------
pdf(opt$pdf, width = 10, height = 6)

for (chr in unique(df$chrom)) {
  chr_df <- df %>% filter(chrom == chr)
  chr_regions <- all_regions %>% filter(chrom == chr)

  p <- ggplot(chr_df, aes(x = window_start, y = strand_count_asym)) +
    geom_line(color = "grey30") +
    labs(
      title = paste("Strand Count Asymmetry with Sliding Outlier Regions -", chr),
      x = "Window Start Position (bp)",
      y = "Strand Count Asymmetry"
    ) +
    theme_minimal(base_size = 14)

  if (nrow(chr_regions) > 0) {
    p <- p +
      geom_rect(
        data = chr_regions,
        aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = region_type),
        alpha = 0.3,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c(high = "tomato", low = "steelblue"), name = "Region Type")
  }

  print(p)
}

dev.off()
cat("PDF plots saved to:", opt$pdf, "\n")
