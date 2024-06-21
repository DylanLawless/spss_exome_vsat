# Set the working directory (adjust the path as necessary)
setwd("~/web/ACMGuru/data")

# Load the data from the CSV file
df <- read.csv("example_variant.csv", stringsAsFactors = FALSE)

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

df$chromosome <- df$seqnames
df$genomic_position <- df$start
df$refcodon
df$alt

# Apply the formatting function to each row
df$output <- apply(df, 1, format_with_metadata, metadata)

# Write to CSV with enriched data
write.csv(df$output, "path_to_enriched_output_file.csv", row.names = FALSE, quote = FALSE)
