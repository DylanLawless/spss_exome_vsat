
geneset_MCL_ID <- c(22, 586, 836 )

df_report_main_text <- readRDS(paste0("../../data/ACMGuru_post_ppi/df_report_main_text_", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

df_singlecase_main_text <- readRDS(paste0("../../data/ACMGuru_singlecase/df_report_main_text.Rds"))

library(dplyr)

# Get unique sample.ids from each data frame
unique_samples_ppi <- unique(df_report_main_text$sample.id)
unique_samples_singlecase <- unique(df_singlecase_main_text$sample.id)

# Find the intersection (overlap) of the two lists
overlap_samples <- intersect(unique_samples_ppi, unique_samples_singlecase)

# Text setup for file output
text_sample_report <- paste("Unique sample.ids from df_report_main_text:", "Count:", length(unique_samples_ppi))
text_sample_singlecase <- paste("Unique sample.ids from df_singlecase_main_text:", "Count:", length(unique_samples_singlecase))
text_sample_overlap <- paste("Overlap of sample.ids between the two data frames:", "Count of Overlapping sample.ids:", length(overlap_samples))

# File setup
sample_summary_file <- "../../data/cohort_summary_curated/cohort_summary_overlap_detailed.txt"

# Check if the file exists and remove if it does to start fresh
if (file.exists(sample_summary_file)) {
  file.remove(sample_summary_file)
}

# Open a file connection for writing in append mode
file_conn <- file(sample_summary_file, open = "a")

# Write texts to the file
writeLines(text_sample_report, file_conn)
writeLines(text_sample_singlecase, file_conn)
writeLines(text_sample_overlap, file_conn)

capture.output(print(overlap_samples), file = file_conn)


# Close the file connection
close(file_conn)




# SYMBOL ----

library(dplyr)

# Assuming df_report_main_text and df_singlecase_main_text are your data frames, already loaded

# Get unique SYMBOLs from each data frame
unique_symbols_ppi <- unique(df_report_main_text$SYMBOL)
unique_symbols_singlecase <- unique(df_singlecase_main_text$SYMBOL)

# Find the intersection (overlap) of the two lists of SYMBOLs
overlap_symbols <- intersect(unique_symbols_ppi, unique_symbols_singlecase)

# Output the results including counts for SYMBOLs
print("Unique SYMBOLs from df_report_main_text:")
print(unique_symbols_ppi)
print(paste("Count:", length(unique_symbols_ppi)))

print("Unique SYMBOLs from df_singlecase_main_text:")
print(unique_symbols_singlecase)
print(paste("Count:", length(unique_symbols_singlecase)))

print("Overlap of SYMBOLs between the two data frames:")
print(overlap_symbols)
print(paste("Count of Overlapping SYMBOLs:", length(overlap_symbols)))
