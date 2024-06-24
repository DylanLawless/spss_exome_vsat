library(dplyr)
library(knitr)
library(kableExtra)
library(jsonlite)

# Define the key columns for genetic variation concepts
key_columns <- c("CHROM", "POS", "REF", "ALT")
key_metadata <- paste(key_columns, "Metadata", sep = "_")

# Enhance dataframe with visible metadata columns for easy reference
df_enhanced <- df_report_set %>%
  mutate(
    CHROM_Metadata = "Type: Chromosome; Cardinality: 1:1; Value Set: SNOMED CT: 91272006, LOINC:48000-4",
    POS_Metadata = "Type: Genomic Position; Cardinality: 1:1; Value Set: GENO:0000902",
    REF_Metadata = "Type: Reference Allele; Cardinality: 1:1; Value Set: string",
    ALT_Metadata = "Type: Alternate Allele; Cardinality: 1:1; Value Set: string"
  )

# Generate footnotes dynamically based on key metadata columns
metadata_descriptions <- c(
  "CHROM_Metadata" = "Metadata for CHROM column: Type: Chromosome, Cardinality: 1:1, Value Set: SNOMED CT: 91272006, LOINC:48000-4",
  "POS_Metadata" = "Metadata for POS column: Type: Genomic Position, Cardinality: 1:1, Value Set: GENO:0000902",
  "REF_Metadata" = "Metadata for REF column: Type: Reference Allele, Cardinality: 1:1, Value Set: string",
  "ALT_Metadata" = "Metadata for ALT column: Type: Alternate Allele, Cardinality: 1:1, Value Set: string"
)

# Create the HTML table with footnotes for metadata
html_table <- df_enhanced %>%
  kable("html", escape = FALSE) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"), 
    full_width = FALSE,
    font_size = 14, 
    position = "left"
  ) %>%
  column_spec(1, bold = TRUE) %>% 
  row_spec(0, bold = TRUE, background = "#D3D3D3", color = "black") %>% 
  # add_header_above(c(" " = 1, "Genetic Variation Info" = 4)) %>% 
  scroll_box(width = "100%", height = "500px") 

# Apply red color to key columns and their metadata
for (col in key_columns) {
  col_index <- which(names(df_enhanced) == col)
  html_table <- html_table %>%
    column_spec(col_index, color = "red", bold = TRUE)
}

# Adding footnotes for key metadata columns
html_table <- html_table %>%
  footnote(
    general = metadata_descriptions[key_metadata],
    general_title = "Key Metadata Descriptions",
    symbol = "*"
  )

# Save the HTML table to a file
output_file <- "enhanced_variant_table.html"
writeLines(as.character(html_table), output_file)

# Convert the enhanced dataframe to JSON format
json_data <- toJSON(df_enhanced, pretty = TRUE)
json_output_file <- "enhanced_variant_data.json"
writeLines(json_data, json_output_file)

# Optionally, open the HTML file in the default system browser
if (Sys.info()["sysname"] == "Windows") {
  shell(paste("start", output_file))
} else {
  system(paste("open", output_file))
}

# Print completion messages
cat("The HTML table has been generated and saved successfully.\n")
cat("The JSON data has been generated and saved successfully.\n")









































# version 0 ----
# get example variant
output_directory <- "../../data/"
df_report_set <- readRDS(file=paste0(output_directory, "example_variant.Rds"))
df_report_set |> names()
df <- df_report_set

# Define metadata (as described previously)
metadata <- list(
  genetic_variation = list(
    type = "concept",
    cardinality = NA,
    value_set = "SO:0001060 |sequence variant|; GENO:0000476 |variant|"
  ),
  genomic_position = list(
    type = "composedOf",
    cardinality = "1:1",
    value_set = "GENO:0000902 |genomic feature location|"
  ),
  chromosome = list(
    type = "composedOf",
    cardinality = "1:1",
    value_set = "SNOMED CT: 91272006 |Chromosome (cell structure)|; LOINC:48000-4 |Chromosome|"
  ),
  reference_allele = list(
    type = "composedOf",
    cardinality = "1:1",
    value_set = "string"
  ),
  alternate_allele = list(
    type = "composedOf",
    cardinality = "1:1",
    value_set = "string"
  )
)

# Function to format data with metadata
format_with_metadata <- function(data_row, metadata) {
  fields <- names(data_row)
  output <- vector("list", length(fields))

  for (i in seq_along(fields)) {
    field_name <- fields[i]
    output[[i]] <- sprintf('"%s | %s | %s | %s"',
                           data_row[[field_name]],
                           metadata[[field_name]]$type,
                           metadata[[field_name]]$cardinality,
                           metadata[[field_name]]$value_set)
  }

  return(paste(output, collapse = ","))
}

names(metadata)
names(df)

# match concept names
df$chromosome <- df$CHROM
df$genomic_position <- df$POS
df$reference_allele <- df$REF
df$alternate_allele <- df$ALT


# simple structure table ----

# Using your existing script to format each variant with its metadata
df_formatted <- apply(df, 1, format_with_metadata, metadata)

# Convert formatted data back into a dataframe
df_table <- data.frame(row_id = 1:nrow(df), formatted_data = df_formatted, stringsAsFactors = FALSE)

# Write to CSV
write.csv(df_table, file = paste0(output_directory, "formatted_variant_data.csv"), row.names = FALSE)


# Variant Metadata Report ----
df_report_set <- readRDS(file=paste0(output_directory, "example_variant.Rds"))
df_report_set |> names()

# Add a metadata description to each field in the dataframe
df_report_set <- df_report_set %>%
  mutate(
    CHROM = paste(CHROM, " | Chromosome | 1:1 | SNOMED CT: 91272006 |Chromosome (cell structure)|; LOINC:48000-4 |Chromosome|"),
    POS = paste(POS, " | Genomic Position | 1:1 | GENO:0000902 |genomic feature location|"),
    REF = paste(REF, " | Reference Allele | 1:1 | string"),
    ALT = paste(ALT, " | Alternate Allele | 1:1 | string")
  )

kable(df_report_set)

# readbale ---
library(DT)

datatable(
  df_report_set,
  options = list(
    columnDefs = list(
      list(targets = c(2, 3, 4, 5),  # Assuming these are the positions of CHROM, POS, REF, ALT
           render = JS(
             "function(data, type, row, meta){
                return type === 'display' ? 
                '<span title=\"' + data + '\">' + data.split(' | ')[0] + '</span>' : 
                data;
              }"
           ))
    )
  )
)


# cols ---
# Add metadata as separate columns for clarity
df_report_set <- df_report_set %>%
  mutate(
    Chromosome_Metadata = "Chromosome | 1:1 | SNOMED CT: 91272006 |Chromosome (cell structure)|; LOINC:48000-4 |Chromosome|",
    Position_Metadata = "Genomic Position | 1:1 | GENO:0000902 |genomic feature location|",
    Ref_Metadata = "Reference Allele | 1:1 | string",
    Alt_Metadata = "Alternate Allele | 1:1 | string"
  )

# Display with kable
table <- kable(df_report_set, format = "html")  # Use format="html" to enable HTML rendering if needed


library(knitr)
library(kableExtra)

# Assuming df_report_set is already created and available
html_table <- kable(df_report_set, format = "html", escape = FALSE) #%>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Save the HTML table to a file
output_file <- "./table.html"
writeLines(html_table, output_file)

# Optionally, open the HTML file in the default system browser
if (Sys.info()["sysname"] == "Windows") {
  shell(paste("start", output_file))
} else {
  system(paste("open", output_file))
}



