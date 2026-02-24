library(tidyverse)
library(stringr)

# 1. Load data
df <- read.csv(
  "C://Users//admin//Documents//internship//pathway_enrichment//enrichment_long_format_with_quantity.csv",
  stringsAsFactors = FALSE
)

# 2. Fix species names
df$species <- df$species %>%
  str_to_lower() %>%
  recode(
    "borangutan" = "Bornean Orangutan",
    "sorangutan" = "Sumatran Orangutan"
  )

# Capitalize remaining species properly
df$species <- ifelse(
  df$species %in% c("Bornean Orangutan", "Sumatran Orangutan"),
  df$species,
  str_to_title(df$species)
)

# 3. Enforce species order
species_order <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean Orangutan",
  "Sumatran Orangutan"
)

df$species <- factor(df$species, levels = species_order)

# 4. Create fg-bg pair
df$fg_bg_pair <- paste(df$foreground, df$background, sep = " | ")

# 5. Loop over fg-bg combinations
fg_bg_pairs <- unique(df$fg_bg_pair)

for (pair in fg_bg_pairs) {
  
  df_sub <- df %>% filter(fg_bg_pair == pair)
  
  # Keep pathways enriched in at least one species
  enriched_pathways <- df_sub %>%
    group_by(property) %>%
    summarise(any_enrichment = any(enrichment_quantity > 0), .groups = "drop") %>%
    filter(any_enrichment) %>%
    pull(property)
  
  if (length(enriched_pathways) == 0) next
  
  df_plot <- df_sub %>% filter(property %in% enriched_pathways)
  
  # Wrap pathway names (20 chars, multiline)
  df_plot$property_wrapped <- str_wrap(df_plot$property, width = 10)
  
  # -------------------------
  # Square heatmap (pink gradient)
  # -------------------------
  p <- ggplot(df_plot, aes(x = species, y = property_wrapped)) +
    geom_tile(aes(fill = enrichment_quantity),
              color = "white") +
    geom_text(
      aes(label = sprintf("%.2f", enrichment_quantity)),
      color = "black",
      size = 3
    )+ 
    scale_fill_gradientn(
      colors = c("#F8F4EC", "#F3B7CE", "#ED79AF", "#E83C91"),
      name = "Enrichment\nQuantity"
    ) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) + 
    labs(
      title = paste("Foreground–Background:", pair),
      x = "Species",
      y = "Pathway"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(hjust = 0.8, vjust = 0.8),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(p)
  
  # 6. Save plot (tall height to avoid y-axis overlap)
  ggsave(
    filename = paste0("C://Users//admin//Documents//internship//pathway_enrichment//results//heatmap_", gsub(" \\| ", "_", pair), ".pdf"),
    plot = p,
    width = 10,
    height = 8 + 1.25*length(species_order),
   dpi = 300
   )
}