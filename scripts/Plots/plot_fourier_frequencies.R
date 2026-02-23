suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript plot_freq.R <input.tsv> <output.pdf> <coding_status>")
}

input_file  <- args[1]
output_pdf  <- args[2]
cd <- args[3]

# -----------------------------
# Smart Reader
# -----------------------------
smart_read <- function(file) {

  df <- tryCatch(
    read.table(file, sep = "\t", header = TRUE,
               fill = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(df) || ncol(df) == 1) {
    df <- read.table(file, header = TRUE,
                     fill = TRUE, stringsAsFactors = FALSE)
  }

  return(df)
}

# -----------------------------
# Filename Parser
# -----------------------------
parse_filename <- function(file_path) {

  fname <- basename(file_path)

  fname <- sub("_with_peaks\\.tsv$", "", fname)

  # Handles BOTH formats
  if (grepl("_", fname)) {

    parts <- strsplit(fname, "_")[[1]]
    species <- parts[1]
    region  <- parts[2]

  } else {

    species <- sub("(upstream|downstream|exitron)$", "", fname)
    region  <- sub(species, "", fname)
  }

  list(species = species, region = region)
}

# -----------------------------
# Local Maxima
# -----------------------------
find_local_maxima <- function(y) {
  which(diff(sign(diff(y))) == -2) + 1
}

# -----------------------------
# Peak Detection
# -----------------------------
compute_peaks <- function(values) {

  dens <- density(values, na.rm = TRUE)

  # Find all local maxima
  local_max_idx <- find_local_maxima(dens$y)

  if (length(local_max_idx) == 0) {
    return(list(
      peak_x = NA, peak_y = NA,
      second_peak_x = NA, second_peak_y = NA,
      low_point_x = NA
    ))
  }

  # Extract peak positions and heights
  peak_data <- data.frame(
    x = dens$x[local_max_idx],
    y = dens$y[local_max_idx]
  )

  # Sort by height (descending)
  peak_data <- peak_data[order(-peak_data$y), ]

  # Peak 1
  peak_x <- peak_data$x[1]
  peak_y <- peak_data$y[1]

  # Peak 2 (first sufficiently separated peak)
  min_sep <- 0.03   # adjust if needed

  second_peak_x <- NA
  second_peak_y <- NA

  if (nrow(peak_data) > 1) {

    for (i in 2:nrow(peak_data)) {

      if (abs(peak_data$x[i] - peak_x) > min_sep) {

        second_peak_x <- peak_data$x[i]
        second_peak_y <- peak_data$y[i]
        break
      }
    }
  }

  # Dip between Peak1 and Peak2 (optional)
  low_point_x <- NA

  if (!is.na(second_peak_x)) {

    between_idx <- which(
      dens$x > min(peak_x, second_peak_x) &
      dens$x < max(peak_x, second_peak_x)
    )

    if (length(between_idx) > 0) {
      low_idx <- between_idx[which.min(dens$y[between_idx])]
      low_point_x <- dens$x[low_idx]
    }
  }

  list(
    peak_x = peak_x,
    peak_y = peak_y,
    second_peak_x = second_peak_x,
    second_peak_y = second_peak_y,
    low_point_x = low_point_x
  )
}


# -----------------------------
# Plot Builder
# -----------------------------
build_plot <- function(values, species, region) {

  peaks <- compute_peaks(values)
  dens  <- density(values, na.rm = TRUE)

  p <- ggplot(data.frame(x = values), aes(x = x)) +
    geom_density(fill = "skyblue", alpha = 0.5) +

    geom_vline(xintercept = peaks$peak_x,
               color = "red", linetype = "dashed", linewidth = 1)

  if (!is.null(peaks$second_peak_x))
    p <- p +
      geom_vline(xintercept = peaks$second_peak_x,
                 color = "blue", linetype = "dashed", linewidth = 1)

  if (!is.null(peaks$low_point_x))
    p <- p +
      geom_vline(xintercept = peaks$low_point_x,
                 color = "green", linetype = "dashed")

  # ---------------------------------------------------
  # CREATE LABEL DATAFRAME
  # ---------------------------------------------------

  label_df <- data.frame(
    x = numeric(0),
    y = numeric(0),
    label = character(0),
    color = character(0),
    stringsAsFactors = FALSE
  )

  # Peak 1
  label_df <- rbind(label_df,
                    data.frame(
                      x = peaks$peak_x,
                      y = peaks$peak_y,
                      label = paste0("Peak1: ", round(peaks$peak_x, 3)),
                      color = "red"
                    ))

  # Peak 2 + Delta
  if (!is.null(peaks$second_peak_x)) {

    label_df <- rbind(label_df,
                      data.frame(
                        x = peaks$second_peak_x,
                        y = peaks$second_peak_y,
                        label = paste0("Peak2: ",
                                       round(peaks$second_peak_x, 3)),
                        color = "blue"
                      ))

    diff_val <- abs(peaks$peak_y - peaks$second_peak_y)

    label_df <- rbind(label_df,
                      data.frame(
                        x = mean(c(peaks$peak_x,
                                   peaks$second_peak_x)),
                        y = max(dens$y),
                        label = paste0("Peak diff = ",
                                       round(diff_val, 3)),
                        color = "black"
                      ))
  }

  # Dip
  if (!is.null(peaks$low_point_x)) {

    dip_y <- dens$y[which.min(abs(dens$x -
                                  peaks$low_point_x))]

    label_df <- rbind(label_df,
                      data.frame(
                        x = peaks$low_point_x,
                        y = dip_y,
                        label = paste0("Dip: ",
                                       round(peaks$low_point_x, 3)),
                        color = "darkgreen"
                      ))
  }

  # ---------------------------------------------------
  # ADD AUTOMATICALLY REPELLED LABELS
  # ---------------------------------------------------

  p <- p +
    geom_text_repel(
      data = label_df,
      aes(x = x, y = y, label = label, color = color),
      show.legend = FALSE,
      size = 4,
      fontface = "bold",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    scale_color_identity()

  p +
    ggtitle(paste(species, "-", region, "-", cd, "-", "Distribution")) +
    xlab("Frequency") +
    ylab("Density") +
    theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
)
}

# -----------------------------
# Main Execution
# -----------------------------
df <- smart_read(input_file)

if (!"Frequencies" %in% colnames(df)) {
  stop("Column 'Frequencies' not found.")
}

values <- as.numeric(unlist(
  strsplit(as.character(df$Frequencies), ";")
))

info <- parse_filename(input_file)

plot <- build_plot(values, info$species, info$region)

ggsave(output_pdf, plot = plot, width = 8, height = 6)

cat("Saved:", output_pdf, "\n")
