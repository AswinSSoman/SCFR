#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(cowplot) # Needed for plot_grid

args <- commandArgs(trailingOnly = TRUE)

# ---------------------------------------------------------
# Input file
# ---------------------------------------------------------
if (length(args) < 2) {
    stop("Usage: Rscript plot.R <input.txt> <output.pdf>")
}

input_file  <- args[1]
output_pdf  <- args[2]
species <- args[3] # Used for the main plot title

# Function to format length values into bp, kb, or mb
format_length <- function(x) {
    # Convert to numeric in case it's not already
    x <- as.numeric(x)

    # Check for MB (1,000,000 bp or more)
    if (x >= 1000000) {
        return(paste0(round(x / 1000000, 1), " mb"))
    }
    # Check for KB (1,000 bp or more, but less than 1,000,000 bp)
    else if (x >= 1000) {
        return(paste0(round(x / 1000, 1), " kb"))
    }
    # Less than 1,000 bp
    else {
        return(paste0(x, " bp"))
    }
}

# -----------------------------------------------------------
# Read input file and preprocess
# -----------------------------------------------------------
df <- read.table(input_file, header = TRUE, sep = "\t")

# Ensure 'Length_threshold' is treated as a factor for discrete plotting, ordering numerically.
df$Length_threshold <- factor(df$Length_threshold, levels = df$Length_threshold[order(as.numeric(as.character(df$Length_threshold)))])

# Apply the formatting function
df$max_kb <- sapply(df$max, format_length)
df$min_kb <- sapply(df$min, format_length)
df$mean_kb <- sapply(df$mean, format_length)

num_categories <- length(levels(df$Length_threshold))

# --- DEFINING THE BROKEN AXIS PARAMETERS ---
# Break the axis into a high plot (just the max values) and a low plot (all other data)
break_limit_low <- 150000 # Upper limit for the main, lower plot
break_limit_high_min <- 150000 # Lower limit for the upper plot
# -------------------------------------------

# Function to create the base plot structure
create_base_plot <- function(data) {
    p <- ggplot(data, aes(
        x = Length_threshold,
        ymin = min, lower = q1, middle = median, upper = q3, ymax = max
    )) +

    # Grey box
    geom_boxplot(
        stat="identity",
        fill="grey85",
        color="black",
        linewidth=0.4,
        width=0.5,
        outlier.shape = NA
    ) +

    # Dashed Light Curve through the means
    geom_line(
        aes(y = mean, group = 1),
        color = "darkgrey",
        linetype = "dashed",
        linewidth = 0.5
    ) +

    # Thicker Q1/Q3 edges
    geom_segment(aes(x=as.numeric(Length_threshold)-0.25,
                     xend=as.numeric(Length_threshold)+0.25,
                     y=q1, yend=q1),
                 linewidth=1.2) +
    geom_segment(aes(x=as.numeric(Length_threshold)-0.25,
                     xend=as.numeric(Length_threshold)+0.25,
                     y=q3, yend=q3),
                 linewidth=1.2) +

    # Thin whiskers + caps (Will be selectively filtered by coord_cartesian)
    geom_segment(aes(x=as.numeric(Length_threshold), xend=as.numeric(Length_threshold),
                     y=min, yend=q1), color="black", linewidth=0.4) +
    geom_segment(aes(x=as.numeric(Length_threshold)-0.12, xend=as.numeric(Length_threshold)+0.12,
                     y=min), color="red", linewidth=0.8) +
    # Q3 to Max (Will show up in p_high or be capped at break_limit_low in p_low)
    geom_segment(aes(x=as.numeric(Length_threshold), xend=as.numeric(Length_threshold),
                     y=q3, yend=max), color="black", linewidth=0.4) +
    geom_segment(aes(x=as.numeric(Length_threshold)-0.12, xend=as.numeric(Length_threshold)+0.12,
                     y=max), color="red", linewidth=0.8) +


    # Median (dark green)
    geom_segment(aes(
        x = as.numeric(Length_threshold)-0.22,
        xend = as.numeric(Length_threshold)+0.22,
        y = median, yend = median
    ), color="darkgreen", linewidth=1.2) +

    # Mean (blue dotted) + point
    geom_point(aes(y=mean, color="Mean"), size=2) +
    geom_segment(aes(
        x=as.numeric(Length_threshold)-0.22, xend=as.numeric(Length_threshold)+0.22,
        y=mean, yend=mean, color="Mean"
    ), linewidth=0.5, linetype="dotted") +

    scale_color_manual(values=c("Mean"="blue")) +

    # Labels --------------------------------

    # Mean label
    geom_text(aes(
        x = as.numeric(Length_threshold) + 0.35,
        y = mean,
        label = paste0(mean_kb)
    ), color="blue", size=3.5) +

    # Min label
    geom_text(aes(
        y = min -250,
        label = paste0(min_kb)
    ), size=3.5) +

    # Max label (Only visible in the high plot)
    geom_text(aes(
        y = max +550,
        label = paste0(max_kb)
    ), size=3.5) +

    labs(
        x = "Length Threshold (bp)",
        y = "SCFR - Gene Desert Overlap (bp)",
        title = species,
        color = ""
    ) +

    theme_bw() +
    theme(
        legend.position = "none"
    )
    return(p)
}

# ------------------------------------
# 1. High Plot (for Max values)
# ------------------------------------
p_high <- create_base_plot(df) +
    # Restrict y-axis to the top part of the break
    coord_cartesian(xlim = c(0.1, num_categories + 0.8),
                    ylim = c(break_limit_high_min, max(df$max) + 50000),
                    expand = FALSE) +
    # Remove all axis elements at the bottom of the plot
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt") # Adjust margins to touch
    )

# Filter out all labels *except* the max labels that appear in this range
p_high_layers <- lapply(p_high$layers, function(layer) {
    if (inherits(layer$geom, "GeomText") && !identical(layer$mapping$y, rlang::sym("max"))) {
        layer$data <- layer$data[FALSE, ] # Remove all non-max labels
    }
    return(layer)
})
p_high$layers <- p_high_layers
p_high <- p_high + labs(title = element_blank(), y = element_blank()) # Remove redundant titles

# ------------------------------------
# 2. Low Plot (for Min, Q1, Median, Q3, Mean)
# ------------------------------------
p_low <- create_base_plot(df) +
    # Restrict y-axis to the main, lower part of the data
    coord_cartesian(xlim = c(0.1, num_categories + 0.8),
                    ylim = c(0, break_limit_low),
                    expand = FALSE) +
    # Keep x-axis elements visible
    theme(
        axis.text.x = element_text(angle=35, hjust=1, size=11),
        axis.text.y = element_text(hjust=1, size=11),
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt") # Adjust margins to touch
    )

# Filter out the max label and any other labels/segments that cross the break point
p_low_layers <- lapply(p_low$layers, function(layer) {
    if (inherits(layer$geom, "GeomText") && identical(layer$mapping$y, rlang::sym("max"))) {
        layer$data <- layer$data[FALSE, ] # Remove max label
    }
    # For segments/points, we rely on coord_cartesian to clip them, but we want to ensure
    # the upper whiskers and max cap are not drawn if they exceed the limit,
    # and the mean/median labels are not drawn if the value exceeds the limit.
    # The 'y' aes in GeomText is what we care about here.
    # The mean label y is filtered by coord_cartesian.
    return(layer)
})
p_low$layers <- p_low_layers


# ----------------------------------------------------------------
# Build external legend (No changes needed here as it's static)
# ----------------------------------------------------------------

legend_plot <- ggplot() +
    geom_segment(aes(x=1, xend=1.4, y=6.2, yend=6.2), color="red", linewidth=1.2) +
    geom_text(aes(x=1.45, y=6.2, label="Min / Max"), hjust=0, size=3.5)

legend_plot <- legend_plot +
    geom_segment(aes(x=1, xend=1.4, y=6.4, yend=6.4), linewidth=1.2) +
    geom_text(aes(x=1.45, y=6.4, label="Q1 / Q3"), hjust=0, size=3.5)

legend_plot <- legend_plot +
    geom_segment(aes(x=1, xend=1.4, y=6.6, yend=6.6), color="darkgreen", linewidth=1.2) +
    geom_text(aes(x=1.45, y=6.6, label="Median"), hjust=0, size=3.5)

legend_plot <- legend_plot +
    geom_segment(aes(x=1, xend=1.4, y=6.8, yend=6.8), color="blue", linetype="dotted", linewidth=0.7) +
    geom_text(aes(x=1.45, y=6.8, label="Mean"), hjust=0, size=3.5)

legend_plot <- legend_plot +
    xlim(1, 3) + ylim(1, 7) +
    theme_void()

# --------------------------------------------------------------
# Combine main plot (p_high + p_low) + external legend
# --------------------------------------------------------------

# Combine the high and low plots vertically
broken_y_axis_plot <- plot_grid(
    p_high, p_low,
    align = "v",
    rel_heights = c(1, 4), # Assign more space to the lower plot
    nrow = 2
)

# Combine the broken axis plot and the legend horizontally
final_plot <- plot_grid(
    broken_y_axis_plot, legend_plot,
    rel_widths = c(6, 1),
    nrow = 1
)

# Set the main title using a title placeholder
title_grob <- ggdraw() + draw_label(
    species,
    fontface = 'bold',
    x = 0,
    hjust = 0
) + theme(
    plot.margin = margin(0, 0, 0, 7)
)
final_plot <- plot_grid(title_grob, final_plot, ncol = 1, rel_heights = c(0.05, 1))

# --------------------------------------------------------------
# SAVE as high-resolution PDF
# --------------------------------------------------------------

ggsave(
    filename = output_pdf,
    plot = final_plot,
    width = 9,
    height = 8,
    units = "in",
    dpi = 800,
    device = cairo_pdf
)

cat("Saved PDF:", output_pdf, "\n")
