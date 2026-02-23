#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# -----------------------------------------------------------
# 1. Read and Transform Data
# -----------------------------------------------------------
input_file  <- "all_species_exon_shadow_length_summary.tsv"
#output_pdf  <- "species_shadow_length_plot.pdf"

df <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE)

# Reshape data so each row is one box (SEU or SED)
df_seu <- df %>%
  select(species, starts_with("min-seu"), starts_with("q1-seu"), 
         starts_with("median-seu"), starts_with("q3-seu"), 
         starts_with("max-seu"), starts_with("mean-seu")) %>%
  rename(min = `min-seu`, q1 = `q1-seu`, median = `median-seu`, 
         q3 = `q3-seu`, max = `max-seu`, mean = `mean-seu`) %>%
  mutate(type = "seu")

df_sed <- df %>%
  select(species, starts_with("min-sed"), starts_with("q1-sed"), 
         starts_with("median-sed"), starts_with("q3-sed"), 
         starts_with("max-sed"), starts_with("mean-sed")) %>%
  rename(min = `min-sed`, q1 = `q1-sed`, median = `median-sed`, 
         q3 = `q3-sed`, max = `max-sed`, mean = `mean-sed`) %>%
  mutate(type = "sed")

df_plot <- bind_rows(df_seu, df_sed)

# Convert species to factor for consistent X-axis ordering
df_plot$species <- factor(df_plot$species)
df_plot$type <- factor(df_plot$type, levels = c("seu", "sed"))

# Calculate unique X-positions for the 'dodged' parallel look
df_plot <- df_plot %>%
  mutate(x_pos = as.numeric(species) + ifelse(type == "seu", -0.2, 0.2))

# Formatting labels (bp, kb, mb)
format_length <- function(x) {
  x <- as.numeric(x)
  if (is.na(x)) return("")
  if (x >= 1000000) return(paste0(round(x / 1000000, 1), " mb"))
  if (x >= 1000) return(paste0(round(x / 1000, 1), " kb"))
  return(paste0(x, " bp"))
}

df_plot$max_label  <- sapply(df_plot$max, format_length)
df_plot$min_label  <- sapply(df_plot$min, format_length)
df_plot$mean_label <- sapply(df_plot$mean, format_length)
df_plot$median_label <- sapply(df_plot$median, format_length)


# ------------------------------------
# 2. Build main plot
# ------------------------------------
p <- ggplot(df_plot, aes(x = x_pos)) +
  
  # FIX: Added 'group = x_pos' to ensure each row becomes its own box
  geom_boxplot(
    aes(ymin = pmax(min, 1), lower = q1, middle = median, 
        upper = q3, ymax = max, fill = type, group = x_pos),
    stat = "identity",
    color = "black",
    linewidth = 0.4,
    width = 0.35,
    outlier.shape = NA
  ) +
  
  # Connecting dashed lines (grouping by 'type' to connect species)
  geom_line(
    aes(y = mean, group = type),
    color = "darkgrey",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  # Red Caps (Min/Max)
  geom_segment(aes(x=x_pos, xend=x_pos, y=pmax(min, 1), yend=q1), linewidth=0.4) +
  geom_segment(aes(x=x_pos-0.1, xend=x_pos+0.1, y=pmax(min, 1), yend=pmax(min, 1)), color="red", linewidth=0.8) +
  geom_segment(aes(x=x_pos, xend=x_pos, y=q3, yend=max), linewidth=0.4) +
  geom_segment(aes(x=x_pos-0.1, xend=x_pos+0.1, y=max, yend=max), color="red", linewidth=0.8) +
  
  # Q1/Q3 thick edges
  geom_segment(aes(x=x_pos-0.17, xend=x_pos+0.17, y=q1, yend=q1), linewidth=1.1) +
  geom_segment(aes(x=x_pos-0.17, xend=x_pos+0.17, y=q3, yend=q3), linewidth=1.1) +
  
  # Median (dark green)
  geom_segment(aes(x = x_pos-0.15, xend = x_pos+0.15, y = median, yend = median), 
               color="darkgreen", linewidth=1.2) +
  
  # Mean (blue point + dotted line)
  geom_point(aes(y=mean), color="blue", size=1.8) +
  geom_segment(aes(x=x_pos-0.15, xend=x_pos+0.15, y=mean, yend=mean), 
               color="blue", linewidth=0.5, linetype="dotted") +
  
  # Annotations (Multiplicative offsets for log scale)
  geom_text(aes(y = max * 1.25, label = max_label), size=2.8, vjust=0, color="red") +
  geom_text(aes(y = pmax(min, 1) * 0.75, label = min_label), size=2.8, vjust=1, color="red") +
  geom_text(aes(x = x_pos + 0.2, y = mean * 1.5, label = mean_label), color="blue", size=2.8, vjust=0) +
  geom_text(aes(y = median * 0.75, label = mean_label), color="darkgreen", size=2.8, vjust=0) +
  
  scale_y_log10() +
  scale_x_continuous(breaks = 1:length(levels(df_plot$species)), 
                     labels = levels(df_plot$species)) +
  scale_fill_manual(values = c("seu" = "grey90", "sed" = "grey70")) +
  
  labs(x = "Species", y = "Exon Shadow Length (bp)", title = "Exon Shadow Length: seu vs sed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1), legend.position = "none")

p
# --------------------------------------------------------------
# 3. Legend and Output
# --------------------------------------------------------------
legend_plot <- ggplot() +
  geom_rect(aes(xmin=1, xmax=1.3, ymin=6.8, ymax=7.0), fill="grey90", color="black") +
  geom_text(aes(x=1.35, y=6.9, label="Upstream"), hjust=0) +
  geom_rect(aes(xmin=1, xmax=1.3, ymin=6.5, ymax=6.7), fill="grey70", color="black") +
  geom_text(aes(x=1.35, y=6.6, label="Downstream"), hjust=0) +
  geom_segment(aes(x=1, xend=1.3, y=6.2, yend=6.2), color="red", linewidth=1.2) +
  geom_text(aes(x=1.35, y=6.2, label="Min/Max"), hjust=0) +
  geom_segment(aes(x=1, xend=1.3, y=5.9, yend=5.9), color="darkgreen", linewidth=1.2) +
  geom_text(aes(x=1.35, y=5.9, label="Median"), hjust=0) +
  geom_segment(aes(x=1, xend=1.3, y=5.6, yend=5.6), color="blue", linetype="dotted", linewidth=1) +
  geom_text(aes(x=1.35, y=5.6, label="Mean"), hjust=0) +
  xlim(1, 2.5) + ylim(5, 7.5) + theme_void()

final_plot <- plot_grid(p, legend_plot, rel_widths = c(9, 1))
final_plot
#ggsave(output_pdf, final_plot, width = 11, height = 7)
