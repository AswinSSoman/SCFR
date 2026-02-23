# Load necessary libraries
library(tidyverse)
library(scales) # For a better size scale if needed, though not strictly required here

# --- 1. Data Loading ---
# Read the TSV file. Do not use na.strings for 'inf' to allow R to convert it to numeric Inf.
data <- read.delim("all_species_summary_gene_deserts.tsv", sep = "\t", na.strings = c("-nan"))

# --- 2. Data Preparation and Cleaning ---
data_clean <- data %>%
  # Rename columns for clarity
  rename(
    Species_Abbr = Species,
    SCFR_Length = Query,
    Odds_Ratio = ratio,
    P_Value = two_tail_pvalue
  ) %>%
  
  # Rename species
  mutate(
    Species = case_when(
      Species_Abbr == "sorangutan" ~ "Sumatran orangutan",
      Species_Abbr == "borangutan" ~ "Bornean orangutan",
      TRUE ~ str_to_title(Species_Abbr) # Capitalize other species names (e.g., human -> Human)
    )
  ) %>%
  
  # Ensure Odds_Ratio is numeric. This step converts the string "inf" to numeric Inf.
  mutate(Odds_Ratio = as.numeric(Odds_Ratio))

# --- 3. Order Species (Y-axis) based on Odds Ratio at 7500 bp threshold ---
# Filter for 7500 bp threshold, then sort by Odds_Ratio descending.
species_order_df <- data_clean %>%
  filter(SCFR_Length == 7500) %>%
  # Arrange by Odds_Ratio, handling Inf by placing it at the top
  arrange(desc(Odds_Ratio))

# Extract the ordered species names
species_order <- species_order_df$Species

# Apply the new factor order to the Species column
data_clean$Species <- factor(data_clean$Species, levels = species_order)

# --- 4. Imputation and Transformation (Same as before, handling inf/nan/0) ---

# Impute Ratios and P-values for non-interpretable thresholds
data_clean <- data_clean %>%
  mutate(
    # Impute P-Value to 1 for non-interpretable thresholds (NegLog10_P_Value = 0)
    P_Value_Imputed = case_when(
      is.na(Odds_Ratio) | Odds_Ratio == 0 ~ 1,
      TRUE ~ P_Value
    ),
    # Impute OR to 1 for non-interpretable thresholds (Log2_Odds_Ratio = 0)
    Odds_Ratio_Imputed = case_when(
      is.na(Odds_Ratio) | Odds_Ratio == 0 ~ 1,
      TRUE ~ Odds_Ratio
    )
  ) %>%
  
  # Calculate transformed metrics
  mutate(
    Log2_Odds_Ratio = log2(Odds_Ratio_Imputed),
    NegLog10_P_Value = -log10(P_Value_Imputed)
  )

# --- Cap the infinite Log2_Odds_Ratio for plotting (Gorilla 5000 bp) ---
max_finite_log2_or <- max(data_clean$Log2_Odds_Ratio[is.finite(data_clean$Log2_Odds_Ratio)], na.rm = TRUE)
cap_value <- max_finite_log2_or + 1 
data_clean$Log2_Odds_Ratio[!is.finite(data_clean$Log2_Odds_Ratio)] <- cap_value 

# Convert SCFR_Length to a factor in the desired order
length_order <- c("500", "1000", "2500", "5000", "7500", "10000")
data_clean$SCFR_Length <- factor(data_clean$SCFR_Length, levels = length_order)

# --- 5. Plotting with increased font size ---
# Set the desired font size
axis_font_size <- 14
title_font_size <- 16
legend_font_size <- 12

p <- ggplot(data_clean, aes(x = SCFR_Length, y = Species)) +
  # Use geom_tile for the background color (OR)
  geom_tile(aes(fill = Log2_Odds_Ratio), color = "gray60", size = 0.5) +
  # Add geom_point for significance (P-value)
  geom_point(aes(size = NegLog10_P_Value), shape = 21, color = "black") +
  
  # Color scale for the Odds Ratio (Log2-transformed)
  scale_fill_gradient2(
    low = "blue",      
    mid = "white",     
    high = "red",      
    midpoint = 0,      
    name = expression(log[2](OR))
  ) +
  
  # Size scale for the P-value significance
  scale_size_continuous(
    name = expression(-log[10](italic(P))),
    range = c(1, 8), 
    breaks = c(0, 5, 20, 50, 100, 150), 
    limits = c(0, max(data_clean$NegLog10_P_Value, na.rm=TRUE)) 
  ) +
  
  # Labels and themes
  labs(
    title = "SCFR Enrichment/Depletion within Gene Deserts Across Primates",
    x = "SCFR Length Threshold (bp)",
    y = "Primate Species (Ordered by 7.5 kb OR)"
  ) +
  theme_minimal() +
  theme(
    # Increase X and Y axis text/label font size
    axis.text = element_text(size = axis_font_size),
    axis.title = element_text(size = axis_font_size, face = "bold"),
    
    # Increase title font size
    plot.title = element_text(hjust = 0.5, face = "bold", size = title_font_size),
    
    # Increase legend title/text font size
    legend.title = element_text(face = "bold", size = legend_font_size),
    legend.text = element_text(size = legend_font_size),
    
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate X labels
  )

# Print the plot object to display it
#print(p)

# Example to save the plot to a file (optional)
ggsave("fishers_test_all_species_gene_deserts.png", plot = p, width = 8, height = 6, dpi = 300)
