#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript pca_variance_plot.R explained_variance.tsv output_plot.pdf species length")
}

input_file  <- args[1]
output_file <- args[2]
sp <- args[3]
ln <- as.numeric(args[4])

# -----------------------------
# Read Data
# -----------------------------
df <- read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(PC_num = as.numeric(gsub("PC", "", PC))) %>%
  arrange(PC_num)

n <- nrow(df)

# -----------------------------
# Broken-stick expectation
# -----------------------------
broken_stick <- sapply(1:n, function(k) {
  sum(1/(k:n)) / n
})

df$broken_stick <- broken_stick
df$retained <- df$variance_explained > df$broken_stick

n_retained_total <- sum(df$retained)
n_retained_top10 <- sum(df$retained[1:10])

# -----------------------------
# Plot top 10
# -----------------------------
plot_df <- df %>%
  slice(1:10) %>%
  mutate(
    PC_label = factor(paste0("PC", PC_num),
                      levels = paste0("PC", sort(PC_num))),
    var_label = percent(variance_explained, accuracy = 0.1),
    cum_label = percent(cumulative_variance, accuracy = 0.1)
  )

# -----------------------------
# Gentle conditional label offsets
# -----------------------------
threshold <- 0.01  # 1% closeness threshold

plot_df <- plot_df %>%
  mutate(
    var_close = c(FALSE, abs(diff(variance_explained)) < threshold),
    cum_close = c(FALSE, abs(diff(cumulative_variance)) < threshold),

    var_offset = ifelse(var_close,
                        ifelse(PC_num %% 2 == 0, 0.015, 0.028),
                        0.012),

    cum_offset = ifelse(cum_close,
                        ifelse(PC_num %% 2 == 0, 0.018, 0.032),
                        0.015)
  )

scale_factor <- max(plot_df$variance_explained) /
                max(plot_df$cumulative_variance)

y_max <- max(plot_df$variance_explained,
             plot_df$broken_stick,
             plot_df$cumulative_variance * scale_factor)

y_limit <- y_max * 1.25

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(plot_df, aes(x = PC_label)) +

  # Bars (same hue, darker if retained)
  geom_col(aes(y = variance_explained,
               fill = retained),
           width = 0.7,
           alpha = 0.95) +

  geom_line(aes(y = broken_stick,
                group = 1,
                color = "Broken-stick expectation"),
            linewidth = 1.2,
            linetype = "dashed") +

  geom_line(aes(y = cumulative_variance * scale_factor,
                group = 1,
                color = "Cumulative variance"),
            linewidth = 1.4) +

  geom_point(aes(y = cumulative_variance * scale_factor,
                 color = "Cumulative variance"),
             size = 3) +

  # Individual variance labels (very close to bars)
  geom_text(
    aes(y = variance_explained + var_offset,
        label = var_label),
    size = 3.8
  ) +

  # Cumulative variance labels (just above node circle)
  geom_text(
    aes(y = cumulative_variance * scale_factor + cum_offset,
        label = cum_label),
    size = 3.8,
    color = "#B22222"
  ) +

  scale_y_continuous(
    limits = c(0, y_limit),
    name = "Individual Variance Explained",
    labels = percent_format(accuracy = 1),
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Cumulative Variance",
                        labels = percent_format(accuracy = 1))
  ) +

  # Gradient-style coloring
  scale_fill_manual(
    values = c(
      "TRUE" = "#1B5E85",   # darker blue (retained)
      "FALSE" = "#5DA9E9"   # lighter blue
    ),
    labels = c("Below Broken-stick", "Exceeds Broken-stick"),
    name = ""
  ) +

  scale_color_manual(values = c(
    "Cumulative variance" = "#B22222",
    "Broken-stick expectation" = "#222222"
  )) +

  labs(
    title = paste0(sp, " - SCFR length >= ", ln),
    subtitle = paste0(
      n_retained_total, " PCs retained overall (Observed > Broken-stick)\n",
      n_retained_top10, " of top 10 exceed Broken-stick expectation"
    ),
    x = "Principal Component",
    color = ""
  ) +

  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(output_file, p, width = 9, height = 7, dpi = 300)

cat("Done. Retained PCs (Broken-stick):", n_retained_total, "\n")
