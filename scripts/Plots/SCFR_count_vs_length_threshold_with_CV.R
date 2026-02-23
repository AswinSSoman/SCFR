# Load the necessary libraries
library(ggplot2)
library(readr)
library(scales) 
library(dplyr) 

# --- 1. Load the Data ---
data_file <- "/home/abhishek/Downloads/SCFR_Supplmentary_Table.xlsx - Supplementary Table 5.csv"
df <- read_csv(data_file, skip = 2)

# --- 2. Data Preprocessing ---
df$Window <- as.factor(df$Window)

# Apply the full species name mapping and order
species_map <- c(
  "human" = "Human", "chimpanzee" = "Chimpanzee", "bonobo" = "Bonobo", "gorilla" = "Gorilla",
  "sorangutan" = "Sumatran orangutan", "borangutan" = "Bornean orangutan", "gibbon" = "Gibbon"
)
df$Species <- dplyr::recode(df$Species, !!!species_map)
species_order_full <- c("Human", "Chimpanzee", "Bonobo", "Gorilla", 
                        "Sumatran orangutan", "Bornean orangutan", "Gibbon")
df$Species <- factor(df$Species, levels = species_order_full)

# --- 3. Calculate Coefficient of Variation (CV) per Window ---
cv_df <- df %>%
  group_by(Window) %>%
  summarise(
    # Calculate CV (StDev / Mean) * 100
    CV = sd(Total_No_SCFR) / mean(Total_No_SCFR) * 100,
    # Determine the Y-position for the label (max + 5% of range)
    y_max = max(Total_No_SCFR) + (max(Total_No_SCFR) - min(Total_No_SCFR)) * 0.05
  ) %>%
  ungroup() %>%
  # Create the label string for display
  mutate(CV_Label = paste0("CV: ", format(CV, digits = 2), "%"))

# Set fixed X-position for the CV label (e.g., last species, 'Gibbon')
cv_df$Species <- "Gibbon" 

# --- 4. Define Professional Color Palette (8 colors) ---
professional_colors_8 <- c(
  "0" = "#000000", "100" = "#1F456E", "500" = "#8B4513", "1000" = "#4682B4", 
  "2500" = "#B22222", "5000" = "#696969", "7500" = "#008B8B", "10000" = "#D2B48C" 
)

# --- 5. Create the Faceted Plot ---
p <- ggplot(
  df, 
  aes(
    x = Species, 
    y = Total_No_SCFR
  )
) +
  # Add points and lines
  geom_point(aes(color = Window), size = 2) + 
  geom_line(aes(color = Window, group = Window), linewidth = 1) + 
  
  # 🚨 ADD CV LABEL to each facet
  geom_text(
    data = cv_df, 
    aes(y = y_max, label = CV_Label), 
    x = "Gibbon", # Position the label near the end of the X-axis
    hjust = 1,    # Align to the right
    vjust = 0,    # Push label slightly up
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  # Apply professional color palette
  scale_color_manual(values = professional_colors_8) +
  
  # Continuous (Linear) Scale with Dynamic Expansion
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = 0.1) 
  ) +
  
  # FACETING: Create 8 separate plots in one panel, allowing free Y-axis scales
  facet_wrap(~ Window, scales = "free_y", ncol = 4, 
             labeller = labeller(Window = function(x) paste(x, "bp"))) +
  
  # Add Labels and Title
  labs(
    title = "Total Number of SCFRs Across Species by Window Stretch with CV",
    x = "Species",
    y = "SCFR Count (Linear Scale)",
    color = "Window Stretch (bp)"
  ) +
  
  # Professional Theme with Centering Adjustments
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5), 
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# --- 6. Save the Plot ---
ggsave("SCFR_Faceted_Plot_FINAL_WITH_CV.png", plot = p, width = 14, height = 8, units = "in", dpi = 300)

print("Plot saved as SCFR_Faceted_Plot_FINAL_WITH_CV.png")
print(p)