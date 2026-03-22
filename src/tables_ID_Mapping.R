library(dplyr)
library(stringr)
set.seed(666)

# This is run from `sync_publication.sh`
# It gathers all files for publication then makes a list of study IDs.
# Filenames for publication are in `./sync_list_table_raw`.
# These are then replaced with a sequential count; 1,2,3...n.
# The new tables are send to the publication dir and a file of the translation is saved.

# Function to extract file paths from the sync list
get_file_paths <- function(file_list_path, base_dir = "../data") {
  # Read lines from the file list
  lines <- readLines(file_list_path)

  # Filter out comment lines and empty lines
  valid_lines <- lines[!grepl("^#", lines) & lines != ""]

  # Prepend the base directory to create full paths
  file_paths <- paste0(base_dir, "/", valid_lines)

  return(file_paths)
}

# Path to your sync list
sync_list_path <- "./sync_list_table_raw"  # Update this path to the actual location

# Call the function to get the file paths
files <- get_file_paths(sync_list_path)
#
# # Initialize a vector to store unique IDs from all files
# all_unique_ids <- character()
#
# # Loop through each file to collect all unique sample_IDs
# for (file in files) {
#   data <- read.csv(file)
#   # data <- read.table(file)
#   # all_unique_ids <- union(all_unique_ids, unique(data[[2]]))  # Assuming sample_IDs are in column 2
#   all_unique_ids <- union(all_unique_ids, unique(data$sample.id))
# }
#
# # Generate a global mapping of unique sample_IDs to sequential pseudo_IDs
# id_mapping <- setNames(seq_along(all_unique_ids), all_unique_ids)
#
# # Convert the mapping to a data frame for saving
# mapping_df <- data.frame(Original_ID = names(id_mapping), Index_ID = id_mapping, row.names = NULL)
#
# # Save the mapping to a CSV file
# write.csv(mapping_df, "../data/ID_Mapping.csv", row.names = FALSE, quote = FALSE)
#
# # Apply the ID mapping to each file
# for (file in files) {
#   data <- read.csv(file)
#
#   # Fix any sample id column
#   colnames(data)[colnames(data) == 'Sample.ID'] <- 'sample.id'
#   print(names(data))
#
#   data$sample.id <- id_mapping[data$sample.id]
#   print(is.na(data$sample.id))
#   # Write the modified dataframe back to a CSV file
#   write.csv(data, paste0(sub(".csv", "", file), "_indexedstudyIDs.csv"), row.names = FALSE, quote = FALSE)
# }
#
#

# When reading and unifying IDs, trim whitespace:
all_unique_ids <- character()

for (file in files) {
  # data <- read.csv(file)

  data <- read.table(file, sep =  ",", header = TRUE)
  # Standardise column name
  colnames(data)[colnames(data) == 'Sample.ID'] <- 'sample.id'
  # Trim spaces and convert to lower case (if appropriate)
  data$sample.id <- str_trim(data$sample.id)
  # Optionally, if sample IDs should be case-insensitive, convert them to lower case:
  data$sample.id <- tolower(data$sample.id)
  
  all_unique_ids <- union(all_unique_ids, unique(data$sample.id))
}

# Generate the mapping with cleaned IDs
id_mapping <- setNames(seq_along(all_unique_ids), all_unique_ids)

# Save mapping for reference
mapping_df <- data.frame(Original_ID = names(id_mapping), Index_ID = id_mapping, row.names = NULL)
write.csv(mapping_df, "../data/ID_Mapping.csv", row.names = FALSE, quote = FALSE)

# Apply the mapping to each file
for (file in files) {
  # data <- read.csv(file)
  data <- read.table(file, sep =  ",", header = TRUE)
  # Standardise column name
  colnames(data)[colnames(data) == 'Sample.ID'] <- 'sample.id'
  colnames(data)[colnames(data) == 'X'] <- 'Variant'
  
  # Trim and lower-case as before
  data$sample.id <- tolower(str_trim(data$sample.id))
  
  # Map IDs; any IDs not matching will become NA
  data$sample.id <- id_mapping[data$sample.id]
  
  names(data) <- gsub("\\.", " ", names(data))
  
  # Optionally, print IDs that resulted in NA for debugging:
  if (any(is.na(data$sample.id))) {
    warning("NA values found in sample.id mapping for file: ", file)
    print(unique(data$sample.id[is.na(data$sample.id)]))
  }
  
  # write.csv(data, paste0(sub(".csv", "", file), "_indexedstudyIDs.csv"),
  #           row.names = FALSE, quote = FALSE)

  colnames(data)[colnames(data) == 'sample id'] <- 'Sample ID'  
  write.table(data, file = paste0(sub(".csv", "", file), "_indexedstudyIDs.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
}
