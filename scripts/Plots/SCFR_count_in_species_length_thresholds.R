# Load the necessary libraries
library(ggplot2)
library(readr)
library(scales) 
library(dplyr) 

# --- 1. Load the Data ---
# 🚨 Use the correct path to your Downloads folder
data_file <- "/home/abhishek/Downloads/SCFR_Supplmentary_Table.xlsx - Supplementary Table 5.csv"
df <- read_csv(data_file, skip = 2)

# --- 2. Data Preprocessing (Updated Species Names) ---

# Define the mapping from short names to full professional names
species_map <- c(
  "gibbon" = "Gibbon",
  "gorilla" = "Gorilla",
  "human" = "Human",
  "bonobo" = "Bonobo",
  "chimpanzee" = "Chimpanzee",
  "borangutan" = "Bornean orangutan",
  "sorangutan" = "Sumatran orangutan"
  
)

# Apply the mapping to the Species column
df$Species <- dplyr::recode(df$Species, !!!species_map)

# Define the order of species on the X-axis using the new full names
species_order_full <- c("Gibbon", "Gorilla", "Human", "Bonobo", "Chimpanzee", "Bornean orangutan",
                        "Sumatran orangutan")
df$Species <- factor(df$Species, levels = species_order_full)

# Convert 'Window' to factor for coloring/grouping
df$Window <- as.factor(df$Window)


# --- 3. Create Label Data Frame ---
# Filter the data to get only the last point of each line (Species == 'Gibbon') 
label_df <- df %>% 
  filter(Species == 'Sumatran orangutan')

# --- 4. Define Professional Color Palette ---
# Define a professional palette (8 colors) using shades of black, blue, and brown/red
professional_colors <- c(
  "0" = "#000000",        # Black (0 bp)
  "100" = "#1F456E",      # Dark Navy Blue
  "500" = "#8B4513",      # Saddle Brown
  "1000" = "#4682B4",     # Steel Blue
  "2500" = "#B22222",     # Firebrick Red/Brown
  "5000" = "#696969",     # Dim Gray
  "7500" = "#008B8B",     # Dark Cyan
  "10000" = "#D2B48C"     # Tan/Light Brown
)

# --- 5. Create the Plot ---
p <- ggplot(
  df, 
  aes(
    x = Species, 
    y = Total_No_SCFR, 
    group = Window, 
    color = Window
  )
) +
  # Add points and lines
  geom_point(size = 3) + 
  geom_line(linewidth = 1) + 
  
  # Apply professional color palette
  scale_color_manual(values = professional_colors) +
  
  # Corrected Logarithmic Scale (starts at 1)
  scale_y_log10(
    labels = scales::comma,
    limits = c(1, NA), # Start Y-axis at 1 to handle log transformation correctly
    breaks = scales::log_breaks(n = 10, base = 10) # Set prominent log breaks
  ) +
  
  # Add labels at the end of each line
  geom_text(
    data = label_df,
    aes(label = paste(Window, "bp")),
    size = 4,
    hjust = 0,                      
    nudge_x = 0.1,                 
    show.legend = FALSE            
  ) +
  
  # Expand X-axis limits to prevent labels from being cut off
  scale_x_discrete(expand = expansion(add = c(0.5, 1.5))) +
  
  # Add Labels and Title (Y-axis uses bquote for log subscript)
  labs(
    title = "Total Number of SCFRs Across Species for Different Length Threshold",
    x = "Species",
    # 🚨 CHANGE 2: Using bquote for log subscript in Y-axis title
    y = bquote(bold("SCFR Count " (log[10] * " scale"))),
    color = "Window Stretch (bp)"
  ) +
  
  # Professional Theme (White background, prominent axes, no grids)
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5), 
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

# --- 6. Save the Plot ---
# Save the plot to a high-resolution PNG file
ggsave("SCFR_Line_Plot_Final.png", plot = p, width = 12, height = 7, units = "in", dpi = 300)

print("Plot saved as SCFR_Line_Plot_Final.png")
print(p)