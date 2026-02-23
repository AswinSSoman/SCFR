#!/usr/bin/env Rscript

# =========================================
# GC Content Side-by-Side Boxplot per Species
# (Proper spacing: within species + between species)
# =========================================

library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# -------------------------
# Paths
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plot.R <input.txt> <output.pdf>")
}

input_file  <- args[1]
output_pdf  <- args[2]

COMBINED_TSV <- input_file

# -------------------------
# Read data
# -------------------------
df <- readr::read_tsv(COMBINED_TSV, show_col_types = FALSE)

# -------------------------
# Clean species names
# -------------------------
df$species <- df$species |>
  str_to_title() |>
  recode(
    "Borangutan" = "Bornean Orangutan",
    "Sorangutan" = "Sumatran Orangutan"
  )

# -------------------------
# Fixed order
# -------------------------
SPECIES_ORDER <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean Orangutan",
  "Sumatran Orangutan"
)

df$species <- factor(df$species, levels = SPECIES_ORDER, ordered = TRUE)
df$seq_type <- factor(
  df$seq_type,
  levels = c("upstream", "downstream", "exon", "intron", "exitron")
)

# -------------------------
# Species spacing (BIG GAP)
# -------------------------
species_pos <- seq(
  from = 1,
  by = 2,   # controls gap BETWEEN species
  length.out = length(SPECIES_ORDER)
)
names(species_pos) <- SPECIES_ORDER

# -------------------------
# Offsets (gap WITHIN species)
# -------------------------
offsets <- c(
  upstream   = -0.6,
  exon       = -0.3,
  exitron    =  0.0,
  intron     =  0.3,
  downstream =  0.6
)

df$xpos <- species_pos[as.character(df$species)] + offsets[df$seq_type]

# -------------------------
# Plot
# -------------------------
cap_width <- 0.125

p <- ggplot(df, aes(x = xpos, group = interaction(species, seq_type), fill = seq_type)) +
  
  geom_boxplot(
    stat = "identity",
    aes(
      ymin = .data$min,
      lower = q1,
      middle = .data$median,
      upper = q3,
      ymax = .data$max
    ),
    width = 0.25
  ) +
  
  # Whisker caps
  geom_segment(aes(x = xpos - cap_width, xend = xpos + cap_width,
                   y = .data$min, yend = .data$min),
               color = "red", linewidth = 0.8) +
  geom_segment(aes(x = xpos - cap_width, xend = xpos + cap_width,
                   y = .data$max, yend = .data$max),
               color = "red", linewidth = 0.8) +
  
  # Median & mean
  geom_segment(aes(x = xpos - cap_width, xend = xpos + cap_width,
                   y = .data$median, yend = .data$median),
               color = "darkgreen", linewidth = 0.4) +
  geom_segment(aes(x = xpos - cap_width, xend = xpos + cap_width,
                   y = .data$mean, yend = .data$mean),
               color = "blue", linetype = "dotted", linewidth = 0.75) +
  
  # Labels
  geom_text(aes(y = .data$median, label = round(.data$median, 1)),
            color = "darkgreen", vjust = 1.5, size = 2.4) +
  geom_text(aes(y = .data$mean, label = round(.data$mean, 1)),
            color = "blue", vjust = -0.7, size = 2.4) +
  geom_text(aes(y = .data$min, label = round(.data$min, 1)),
            color = "red", vjust = 1.6, size = 2.4) +
  geom_text(aes(y = .data$max, label = round(.data$max, 1)),
            color = "red", vjust = -0.6, size = 2.4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))+
  
  # N above max
  geom_text(aes(y = Inf, label = paste0(.data$N)),
            vjust = 4, hjust=1, fontface = "bold", size = 2.5, angle=45) +
  
  scale_x_continuous(
    breaks = species_pos,
    labels = sapply(SPECIES_ORDER, function(x) str_wrap(x, 10)),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_fill_manual(
    values = c(
      upstream = "pink",
      exon = "skyblue",
      intron = "khaki",
      downstream = "lightgray",
      exitron = "lightyellow"
    )
  ) +
  
  labs(
    title = args[3],
    x = "Species",
    y = args[4],
    fill = "Sequence type"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(hjust = 0.5, color="black"),
    legend.position = "top"
  )

# Dummy legend data
legend_df <- data.frame(
  statistic = c("Mean", "Median"),
  x = c(1, 2), y = c(0, 0)   # y values don't matter
)

# Add invisible geoms just to create legend
p1 <- p +
  geom_line(
    data = legend_df,
    aes(x = x, y = y, color = statistic, linetype = statistic), 
    inherit.aes = FALSE,
    size = 1,
    alpha = 1  # invisible on plot
  ) +
  scale_color_manual(
    name = "Statistic",
    values = c("Mean" = "blue", "Median" = "darkgreen")
  ) +
  scale_linetype_manual(
    name = "Statistic",
    values = c("Mean" = "dotted", "Median" = "solid")
  )

#print(p)

ggsave(output_pdf, p1, width=14, height=6, dpi=300)
