library(Biostrings)
library(stringr)

# Define the required order of influenza genome segments
segment_order <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")

# Set the date threshold to filter out older samples (applied only if dates exist)
date_threshold <- as.Date("2024-12-05")

# Read the FASTA file containing influenza sequences
input_file <- "four_sequences.fas"
sequences <- readDNAStringSet(input_file)

# Initialize lists to store sample sequences and their corresponding dates
samples <- list()
dates <- list()
has_dates <- FALSE  # Track if any headers contain dates

# Function to extract sample ID, segment, and date from the header
parse_header <- function(header) {
  # Check if header contains '|' (indicating detailed format)
  if (grepl("\\|", header)) {
    parts <- unlist(strsplit(header, "\\|"))
    sample_id <- parts[2]  # Extract sample identifier
    segment <- parts[length(parts)]  # Extract segment name
    
    # Try to extract a date (assuming the second-last field is the date)
    sample_date <- suppressWarnings(as.Date(parts[length(parts) - 1], format = "%Y-%m-%d"))
    if (!is.na(sample_date)) has_dates <<- TRUE  # Track if we have dates
    
  } else { 
    # Simple format: U2501319_0001_PB2
    parts <- unlist(strsplit(header, "_"))
    sample_id <- paste(parts[1:(length(parts) - 1)], collapse = "_")  # All but last part
    segment <- parts[length(parts)]  # Last part is the segment
    sample_date <- NA  # No date available
  }
  
  return(list(sample_id = sample_id, segment = segment, sample_date = sample_date))
}

# Process each sequence in the FASTA file
for (i in seq_along(names(sequences))) {
  header <- names(sequences)[i]  # Extract header
  parsed <- parse_header(header)  # Parse the header
  
  sample_id <- parsed$sample_id
  segment <- parsed$segment
  sample_date <- parsed$sample_date
  
  # Only process valid segments
  if (segment %in% segment_order) {
    # Initialize storage for a new sample
    if (!sample_id %in% names(samples)) {
      samples[[sample_id]] <- setNames(rep(NA_character_, length(segment_order)), segment_order)
      dates[[sample_id]] <- sample_date  # Store date (if available)
    }
    # Store the sequence under the corresponding segment
    samples[[sample_id]][segment] <- as.character(sequences[i])
  }
}

# Filter by date threshold **only if dates are available**
if (has_dates) {
  samples <- samples[!is.na(dates) & dates >= date_threshold]
  dates <- dates[!is.na(dates) & dates >= date_threshold]
}

# Write concatenated sequences to the output FASTA file
output_file <- "foursequences_concatenated.fasta"
output_lines <- c()

for (sample_id in names(samples)) {
  # Extract and sort sequences according to segment_order
  sorted_sequences <- samples[[sample_id]][segment_order]
  concatenated_sequence <- paste(sorted_sequences, collapse = "")
  
  # Construct FASTA header
  if (has_dates && !is.na(dates[[sample_id]])) {
    fasta_header <- paste0(">", sample_id, "|", dates[[sample_id]])
  } else {
    fasta_header <- paste0(">", sample_id)  # No date available
  }
  
  # Append to output lines
  output_lines <- c(output_lines, fasta_header, concatenated_sequence)
}

# Write the final sorted and concatenated sequences to file
writeLines(output_lines, output_file)

# Print a message to indicate successful writing
cat("Sorted and concatenated sequences written to", output_file, "\n")
