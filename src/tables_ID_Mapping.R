library(dplyr)

# List of file names
files <- c("Table_1_ACMGuru_singlecase_df_report_cohort_data.csv", "Table_2_ACMGuru_post_ppi_genetic_df_report_main_text_22_586.csv")

# Initialize a vector to store unique IDs from all files
all_unique_ids <- character()

# Loop through each file to collect all unique sample_IDs
for (file in files) {
  data <- read.csv(file)
  all_unique_ids <- union(all_unique_ids, unique(data[[2]]))  # Assuming sample_IDs are in column 2
}

# Generate a global mapping of unique sample_IDs to sequential pseudo_IDs
id_mapping <- setNames(seq_along(all_unique_ids), all_unique_ids)

# Convert the mapping to a data frame for saving
mapping_df <- data.frame(Original_ID = names(id_mapping), Index_ID = id_mapping, row.names = NULL)

# Save the mapping to a CSV file
write.csv(mapping_df, "ID_Mapping.csv", row.names = FALSE)

# Apply the ID mapping to each file
for (file in files) {
  data <- read.csv(file)
  data[[2]] <- id_mapping[data[[2]]]
  # Write the modified dataframe back to a CSV file
  write.csv(data, paste0(sub(".csv", "", file), "_indexedstudyIDs.csv"), row.names = FALSE)
}
