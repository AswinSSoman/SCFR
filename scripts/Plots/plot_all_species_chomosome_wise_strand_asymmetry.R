# 1. Load the required library
library(tidyverse)

# 2. Define the file name
file_name <- "SCFR_Supplmentary_Table.xlsx - Supplementary Table 4 (2).csv"

# 3. Read the data, skipping the initial header rows (skip = 2)
df <- read_csv(file_name, skip = 2)

# 4. Clean and Transform the data
df_long <- df %>%
  
  # Rename species for clarity (matching capitalization in data)
  mutate(Species = case_when(
    Species == "Borangutan" ~ "Bornean orangutan",
    Species == "Sorangutan" ~ "Sumatran orangutan",
    TRUE ~ Species
  )) %>%
  
  # FILTER: Remove ALL and NA rows from the chromosome column
  filter(!is.na(`Chromosome no.`), `Chromosome no.` != "ALL") %>%
  
  # X-AXIS ORDERING: Create custom numerical key for correct sorting (1, 2, ... X, Y)
  mutate(chr_order = case_when(
    # Sex chromosomes are guaranteed to be last
    `Chromosome no.` == "ChrX" ~ 998, 
    `Chromosome no.` == "ChrY" ~ 999,
    # Convert numeric chromosomes to numbers
    TRUE ~ as.numeric(gsub("Chr", "", `Chromosome no.`))
  )) %>%
  
  # Handle any remaining NA values (which result from non-numeric autosomes like 'ChrU1')
  # These are assigned a high value (900) to sort them after standard autosomes, but before X/Y.
  mutate(chr_order = if_else(is.na(chr_order), 900, chr_order)) %>%
  
  # Sort the data frame by the final robust numerical key
  arrange(Species, chr_order) %>%
  
  # Convert X-axis to factor based on sorted order
  mutate(`Chromosome no.` = factor(`Chromosome no.`, levels = unique(`Chromosome no.`))) %>%
  
  # REMOVE 'CHR' PREFIX FROM LABELS
  mutate(`Chromosome no.` = fct_relabel(`Chromosome no.`, ~ gsub("Chr", "", .x))) %>%
  
  # Remove temporary key and pivot to long format
  select(-chr_order, Species, `Chromosome no.`, strand_count_asym, strand_length_asym) %>%
  pivot_longer(
    cols = c(strand_count_asym, strand_length_asym),
    names_to = "Asymmetry_Type",
    values_to = "Asymmetry_Value"
  ) %>%
  
  # Set factor levels for dot layering (Small dot on top)
  mutate(Asymmetry_Type = factor(Asymmetry_Type, 
                                 levels = c("strand_count_asym", "strand_length_asym"))) %>%
  arrange(Asymmetry_Type) %>% 
  
  # Assign the custom, reduced sizes
  mutate(point_size = case_when(
    Asymmetry_Type == "strand_count_asym" ~ 1.5,  # Big dot
    TRUE ~ 1.0                                  # Small dot
  ))


# 5. Create the ggplot2 dot plot
plot <- ggplot(df_long, aes(x = `Chromosome no.`, y = Asymmetry_Value, 
                            color = Asymmetry_Type, 
                            shape = Asymmetry_Type,
                            size = point_size)) +
  
  geom_point() +
  
  # Custom color scale
  scale_color_manual(
    name = NULL, # Remove Title
    values = c("strand_count_asym" = "#004080", 
               "strand_length_asym" = "#FF7F0E"),
    labels = c("strand_count_asym" = "Strand count asymmetry", 
               "strand_length_asym" = "Strand length asymmetry")
  ) +
  
  # Custom shape scale
  scale_shape_manual(
    name = NULL, 
    values = c("strand_count_asym" = 19, 
               "strand_length_asym" = 1), 
    labels = c("strand_count_asym" = "Strand count asymmetry", 
               "strand_length_asym" = "Strand length asymmetry")
  ) +
  
  scale_size_identity(guide = "none") + 
  
  # Rename Y-AXIS
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1), 
                     name = "Strand asymmetry") +
  
  # Facet by Species
  facet_wrap(~ Species, 
             strip.position = "bottom", 
             scales = "free_x") + 
  
  # Customize appearance and Title
  labs(title = NULL, x = "Chromosome No.") +
  theme_minimal() +
  theme(
    # FINAL CHANGES: Small, spaced, and aligned X-axis labels
    axis.text.x = element_text(
      angle = 90, 
      hjust = 1, 
      vjust = 0.5, 
      size = 7, 
      margin = unit(c(0.2, 0, 0, 0), "cm")
    ), 
    axis.ticks.x = element_line(), 
    
    panel.border = element_rect(color = "lightgray", fill = NA, linewidth = 0.3),
    panel.spacing.x = unit(0.1, "lines"), 
    
    # LEGEND POSITION (Lower Right)
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"), 
    legend.box.background = element_rect(color = "black", size = 0.5), 
    legend.background = element_blank()
  ) +
  
  # Set a consistent size for legend points
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  )

# 6. Display the plot
print(plot)