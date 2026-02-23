suppressPackageStartupMessages({
  library(stats)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript extract_fourier_peaks.R <input.tsv> <output.tsv>")
}

input_file  <- args[1]
output_file <- args[2]

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
# Filename Parser (species only)
# -----------------------------
parse_species <- function(file_path) {

  fname <- basename(file_path)
  fname <- sub("_with_peaks\\.tsv$", "", fname)

  if (grepl("_", fname)) {
    parts <- strsplit(fname, "_")[[1]]
    species <- parts[1]
  } else {
    species <- fname
  }

  return(species)
}

# -----------------------------
# Local Maxima Detector
# -----------------------------
find_local_maxima <- function(y) {
  which(diff(sign(diff(y))) == -2) + 1
}

# -----------------------------
# Peak Detection
# -----------------------------
compute_peaks <- function(values) {

  dens <- density(values, na.rm = TRUE)
  local_max_idx <- find_local_maxima(dens$y)

  if (length(local_max_idx) == 0) {
    return(list(
      peak_x = NA,
      second_peak_x = NA,
      peak_diff = NA,
      low_point_x = NA,
      modality = "No clear peak"
    ))
  }

  peak_data <- data.frame(
    x = dens$x[local_max_idx],
    y = dens$y[local_max_idx]
  )

  peak_data <- peak_data[order(-peak_data$y), ]

  peak_x <- peak_data$x[1]
  peak_y <- peak_data$y[1]

  min_sep <- 0.03

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

  peak_diff <- ifelse(!is.na(second_peak_y),
                      abs(peak_y - second_peak_y),
                      NA)

  modality <- ifelse(is.na(second_peak_x),
                     "Unimodal",
                     "Bimodal")

  list(
    peak_x = peak_x,
    second_peak_x = second_peak_x,
    peak_diff = peak_diff,
    low_point_x = low_point_x,
    modality = modality
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

species <- parse_species(input_file)

peaks <- compute_peaks(values)

summary_df <- data.frame(
  species = species,
  highest_peak = round(peaks$peak_x, 5),
  second_highest_peak = ifelse(is.na(peaks$second_peak_x), NA,
                               round(peaks$second_peak_x, 5)),
  peak_difference = ifelse(is.na(peaks$peak_diff), NA,
                           round(peaks$peak_diff, 5)),
  dip = ifelse(is.na(peaks$low_point_x), NA,
               round(peaks$low_point_x, 5)),
  modal_of_distribution = peaks$modality,
  stringsAsFactors = FALSE
)

write.table(summary_df,
            file = output_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("Summary saved to:", output_file, "\n")
