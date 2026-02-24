suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# -------------------------
# 1. Handle arguments
# -------------------------
if (interactive()) {
  # ---- RStudio usage ----
  INPUT_FILE  <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//exitron_pc_genes_txs//exitron_percentage_in_genes_txs.tsv"
  OUTPUT_FILE <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//exitron_pc_genes_txs//exitron_pc_in_genes_txs.pdf"
} else {
  # ---- Command-line usage ----
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 2) {
    stop("Usage: Rscript plot_code.R <input_file> <output_file>")
  }
  
  INPUT_FILE  <- args[1]
  OUTPUT_FILE <- args[2]
}

# -------------------------
# 2. Load data
# -------------------------
df <- read.delim(INPUT_FILE, stringsAsFactors = FALSE)

# -------------------------
# 3. Species naming + order
# -------------------------
df$Species <- df$Species |>
  str_to_lower() |>
  recode(
    "borangutan" = "Bornean Orangutan",
    "sorangutan" = "Sumatran Orangutan"
  ) |>
  str_to_title()

SPECIES_ORDER <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean Orangutan",
  "Sumatran Orangutan"
)

df$Species <- factor(df$Species, levels = SPECIES_ORDER)

plot_df <- df |>
  select(
    Species,
    pc_genes_with_exitron, pc_genes_wo_exitron,
    pc_txs_with_exitron, pc_txs_wo_exitron
  ) |>
  pivot_longer(
    -Species,
    names_to = c("Feature", "Exitron"),
    names_pattern = "pc_(genes|txs)_(with|wo)_exitron",
    values_to = "Percentage"
  ) |>
  mutate(
    Feature = recode(
      Feature,
      "genes" = "Genes",
      "txs"   = "Transcripts"
    ),
    Exitron = recode(
      Exitron,
      "with" = "With Exitrons",
      "wo"   = "Without Exitrons"
    ),
    Feature = factor(Feature, levels = c("Genes", "Transcripts"))
  )

# -------------------------
# Totals for annotation
# -------------------------
totals_df <- df |>
  select(Species, total_genes, total_tx) |>
  pivot_longer(
    -Species,
    names_to = "Feature",
    values_to = "Total"
  ) |>
  mutate(
    Feature = recode(
      Feature,
      "total_genes" = "Genes",
      "total_tx"    = "Transcripts"
    ),
    Feature = factor(Feature, levels = c("Genes", "Transcripts")),
    y = 115  # above with-exitron % label
  )

# -------------------------
# 5. Plot
# -------------------------
p <- ggplot(
  plot_df,
  aes(x = Feature, y = Percentage, fill = Exitron)
) +
  geom_col(width = 0.5) +
  
  # Total annotations
  geom_text(
    data = totals_df,
    aes(x = Feature, y = y, label = scales::comma(Total)),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
  ) +
  
  # Without Exitron: below bar
  geom_text(
    data = subset(plot_df, Exitron == "Without Exitrons"),
    aes(y = -3, label = sprintf("%.1f", Percentage)),
    size = 3.5
  ) +
  
  # With Exitron: above bar
  geom_text(
    data = subset(plot_df, Exitron == "With Exitrons"),
    aes(y = 103, label = sprintf("%.1f", Percentage)),
    size = 3.5
  ) +
  
  facet_wrap(
    ~ Species,
    nrow = 1,
    strip.position = "bottom",
    labeller = labeller(Species = label_wrap_gen(width = 0))
  ) +
  
  scale_y_continuous(
    limits = c(-8, 120),
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Percentage of Genes and Transcripts With and Without Exitrons",
    x = NULL,
    y = "Percentage (%)",
    fill = "Exitron Status"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(color="black", size=10, angle=45, hjust=1),
    plot.title = element_text(hjust=0, face="bold"),
    legend.title = element_text(face="bold")
  )

print(p)
ggsave(OUTPUT_FILE, plot=p, width=12, height=5, dpi=300)