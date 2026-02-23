#!/usr/bin/env Rscript

############################################
## Parse arguments
############################################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript shadow_length_enrichment.R <input.tsv> <output_prefix> [max_len] [fdr]")
}

INPUT_FILE <- args[1]
OUT_PREFIX <- args[2]
MAX_LEN <- ifelse(length(args) >= 3, as.numeric(args[3]), 500)
FDR_CUTOFF <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)
TOP_N <- 10

############################################
## Load libraries
############################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tidytext)
})

############################################
## Load data
############################################
df <- read_tsv(INPUT_FILE, col_types = cols())

df <- df %>%
  rename(
    upstream = upstream_len_in_scfr,
    downstream = downstream_len_in_scfr
  )

############################################
## Long format
############################################
shadow_long <- df %>%
  select(upstream, downstream) %>%
  pivot_longer(
    cols = everything(),
    names_to = "shadow_type",
    values_to = "shadow_length"
  ) %>%
  filter(!is.na(shadow_length),
         shadow_length >= 0)

############################################
## Frequency counts
############################################
shadow_counts <- shadow_long %>%
  count(shadow_type, shadow_length)

write_tsv(shadow_counts,
          paste0(OUT_PREFIX, "_shadow_length_counts.tsv"))

############################################
## Plot: full range
############################################
p_full <- ggplot(shadow_counts,
                 aes(x = shadow_length, y = n, color = shadow_type)) +
  geom_line() +
  geom_point(size = 0.7) +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Shadow length (bp)",
    y = "Count (log scale)",
    color = "Shadow type",
    title = "Shadow length frequency (full range)"
  )

x_max <- max(shadow_counts$shadow_length, na.rm = TRUE)

p_full <- p_full +
  scale_x_continuous(
    breaks = seq(0, x_max, by = 2000)   # adjust step if needed
  ) +
  theme(
    axis.text.x = element_text(size = 7)
  )

p_full <- p_full +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(paste0(OUT_PREFIX, "_shadow_length_full.pdf"),
       p_full, width = 8, height = 5, device = cairo_pdf)

############################################
## Plot: 1–MAX_LEN
############################################
p_zoom <- ggplot(shadow_counts %>% filter(shadow_length <= MAX_LEN),
                 aes(x = shadow_length, y = n, color = shadow_type)) +
  geom_line() +
  geom_point(size = 0.7) +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Shadow length (bp)",
    y = "Count (log scale)",
    color = "Shadow type",
    title = paste("Shadow length frequency (1-", MAX_LEN, "bp)", sep = "")
  )

p_zoom <- p_zoom +
  scale_x_continuous(
    breaks = seq(0, MAX_LEN, by = 25)   # dense but readable
  ) +
  theme(
    axis.text.x = element_text(size = 7)
  )

p_zoom <- p_zoom +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
  )

p_zoom <- p_zoom +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(paste0(OUT_PREFIX, "_shadow_length_1_", MAX_LEN, ".pdf"),
       p_zoom, width = 7, height = 5, device = cairo_pdf)

############################################
## Enrichment testing
############################################
enrichment_results <- shadow_counts %>%
  filter(shadow_length >= 1,
         shadow_length <= MAX_LEN) %>%
  group_by(shadow_type) %>%
  mutate(
    total = sum(n),
    K = n(),
    expected = total / K,
    enrichment_ratio = n / expected,
    pval = map_dbl(
      n,
      ~ binom.test(
        x = .x,
        n = unique(total),
        p = 1 / unique(K),
        alternative = "greater"
      )$p.value
    )
  ) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(pval, method = "BH")
  ) %>%
  arrange(shadow_type, padj)

write_tsv(enrichment_results,
          paste0(OUT_PREFIX, "_shadow_length_enrichment.tsv"))

############################################
## Significant enrichments
############################################
sig <- enrichment_results %>%
  filter(padj < FDR_CUTOFF)

write_tsv(sig,
          paste0(OUT_PREFIX, "_significant_enrichments.tsv"))

############################################
## Plot: top enriched lengths
############################################
p_top <- ggplot(
  sig %>%
    group_by(shadow_type) %>%
    slice_max(enrichment_ratio, n = TOP_N) %>%
#    slice_min(padj, n = TOP_N) %>%
    ungroup(),
  aes(
    x = reorder_within(shadow_length, -enrichment_ratio, shadow_type),
    y = enrichment_ratio,
    fill = shadow_type
  )
) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~shadow_type, scales = "free_x") +
  scale_x_reordered(labels = function(x) gsub("__.*", "", x)) +
  theme_classic() +
  labs(
    x = "Shadow length (bp)",
    y = "Enrichment ratio",
    title = paste("Top enriched shadow lengths (FDR <", FDR_CUTOFF, ")")
  )

p_top <- p_top +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

p_top <- p_top +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(paste0(OUT_PREFIX, "_top_enriched_shadow_lengths.pdf"),
       p_top, width = 13, height = 5, device = cairo_pdf)

############################################
## Volcano plot
############################################
p_volcano <- ggplot(enrichment_results %>% filter(shadow_length <= MAX_LEN),
                    aes(x = shadow_length,
                        y = -log10(padj),
                        color = shadow_type)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(FDR_CUTOFF),
             linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Shadow length (bp)",
    y = "-log10(FDR)",
    color = "Shadow type",
    title = paste("Shadow length enrichment volcano (0–", MAX_LEN, "bp)", sep = "")
  )

p_volcano <- p_volcano +
  coord_cartesian(
    expand = TRUE,
    clip = "off"
  ) +
  theme(
    plot.margin = margin(5.5, 5.5, 14, 5.5)
  )

p_volcano <- p_volcano +
  scale_x_continuous(
    breaks = seq(0, MAX_LEN, by = 25),   # more ticks
    minor_breaks = seq(0, MAX_LEN, by = 5)
  ) +
  scale_y_continuous(
    minor_breaks = waiver()
  ) +
  theme(
    axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)
  )

p_volcano <- p_volcano +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(paste0(OUT_PREFIX, "_shadow_length_volcano.pdf"),
       p_volcano, width = 7, height = 5, device = cairo_pdf)

############################################
## Plot: top enriched lengths: Window-wise
############################################
WINDOW_SIZE <- 10      # e.g. 5, 10, 20
MAX_SHADOW_LEN <- 200  # e.g. 50, 100, 500

# Define breaks and labels dynamically
breaks_vec <- c(seq(1, MAX_SHADOW_LEN + 1, by = WINDOW_SIZE), Inf)

labels_vec <- c(
  paste0(
    seq(1, MAX_SHADOW_LEN, by = WINDOW_SIZE),
    "–",
    pmin(seq(WINDOW_SIZE, MAX_SHADOW_LEN, by = WINDOW_SIZE), MAX_SHADOW_LEN)
  ),
  paste0(">", MAX_SHADOW_LEN)
)

# Bin shadow lengths
sig_binned <- sig %>%
  mutate(
    shadow_bin = cut(
      shadow_length,
      breaks = breaks_vec,
      labels = labels_vec,
      right = FALSE
    )
  )

# Summarise enrichment per bin
range_enrichment <- sig_binned %>%
  group_by(shadow_type, shadow_bin) %>%
  summarise(
    median_enrichment = median(enrichment_ratio, na.rm = TRUE),
    n_lengths = n(),
    .groups = "drop"
  )

# Plot
p_range <- ggplot(
  range_enrichment,
  aes(x = shadow_bin,
      y = median_enrichment,
      fill = shadow_type)
) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Shadow length range (bp)",
    y = "Median enrichment ratio",
    title = paste(
      "Enrichment of shadow-length ranges (",
      WINDOW_SIZE, "bp bins, FDR <", FDR_CUTOFF, ")",
      sep = ""
    )
  ) +
  theme(
    axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(paste0(OUT_PREFIX, "_top_enriched_shadow_length_bins.pdf"),
       p_range, width = 7, height = 5, device = cairo_pdf)

