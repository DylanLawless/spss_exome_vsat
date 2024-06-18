source("ACMGuru_gene_illustrate_vcurrent.R")
capture.output(sessionInfo(), file = "../session_info.txt")

# Protein Structure ----
# https://nvelden.github.io/NGLVieweR/articles/NGLVieweR.html
# install.packages("NGLVieweR")
library(NGLVieweR)

# uniprot to PDB ----

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")

# Load the library
library(biomaRt)

# List available marts
listMarts()

# Choose the Ensembl Mart
ensembl = useMart("ensembl")

# List available datasets
listDatasets(ensembl)

# Choose the appropriate dataset (for example, 'hsapiens_gene_ensembl' for human)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

# merge grouped_df_max and df_report_position on the 'seqid' column to make a new dataframe (assuming 'seqid' column is common to both dataframes).
grouped_df_max_report_position <- merge(grouped_df_max, df_report_position)

# Define your UniProt IDs
Accessions <- grouped_df_max_report_position$seqid

# Get the PDB IDs
pdb_ids <- getBM(
  filters = "uniprotswissprot", 
  attributes = c("uniprotswissprot", "pdb"),
  values = Accessions,
  mart = ensembl
)

print(pdb_ids)

# get consistent col name for protein ID (uniprot)
# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
colnames(grouped_df_max_report_position)[colnames(grouped_df_max_report_position) == 'seqid'] <- 'uniprotswissprot'

# merge
grouped_df_max_report_position_pdb <- merge(
  grouped_df_max_report_position, 
  pdb_ids,
  by = "uniprotswissprot")

# Find a "best" PDB for each uniprotswissprot ID using longest sequence ----

# Define the URL and output file
# file from SIFTS
url <- "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
output_file <- "pdb_chain_uniprot.tsv"
output_dir <- "../../data/ACMGuru_gene_illustrate_protein_structure/"

# Download the file
# system2("curl", args = c("-v -# -o", paste0(output_dir, output_file, ".gz"), url))

# Read the SIFTS data
sifts_data <- read.csv(paste0(output_dir, output_file), sep = "\t",  header = TRUE, comment.char = "#")

names(sifts_data)

# Filter the SIFTS data to include only the rows that match the IDs in your data
matched_data <- sifts_data %>% 
  filter(SP_PRIMARY %in% grouped_df_max_report_position_pdb$uniprotswissprot)

# First make sure the SP_BEG and SP_END columns are numeric
matched_data <- matched_data %>% 
  mutate(across(c(SP_BEG, SP_END), as.numeric))

# Calculate the chain lengths
matched_data <- matched_data %>% 
  mutate(CHAIN_LENGTH = SP_END - SP_BEG + 1)

# This 3 part code block filters 'matched_data' to select representative protein chain records.
# The selection criteria are: 
# 1) Maximum 'CHAIN_LENGTH' 
# 2) Maximum count of unique 'PDB' entries per 'SP_PRIMARY' 
# 3) Maximum count of unique 'CHAIN' entries per 'PDB'.
# The result is 'longest_chain', retaining the rows with the top values in these three categories per 'SP_PRIMARY'.

# Counting the number of unique PDBs and CHAINs for each SP_PRIMARY and PDB, respectively
matched_data <- matched_data %>%
  add_count(SP_PRIMARY, name = "PDB_COUNT") %>%
  group_by(PDB) %>%
  mutate(CHAIN_COUNT = n_distinct(CHAIN)) %>%
  ungroup()

# Identify the top group for each SP_PRIMARY
top_group <- matched_data %>%
  group_by(SP_PRIMARY) %>%
  arrange(desc(CHAIN_LENGTH), desc(PDB_COUNT), desc(CHAIN_COUNT), .by_group = TRUE) %>%
  filter(row_number() == 1) %>%
  dplyr::select(SP_PRIMARY, PDB)

# Filtering to keep the rows of the top group for each SP_PRIMARY
longest_chain <- matched_data %>%
  semi_join(top_group, by = c("SP_PRIMARY", "PDB"))


colnames(longest_chain)[colnames(longest_chain) == 'SP_PRIMARY'] <- 'uniprotswissprot'

colnames(longest_chain)[colnames(longest_chain) == 'PDB'] <- 'pdb'

# Convert PDB IDs to upper case to match
longest_chain$pdb <- toupper(longest_chain$pdb)

longest_chain$longest_chain <- "Yes"

names(longest_chain)
names(grouped_df_max_report_position_pdb)

# Full outer join of the dataframes
grouped_df_max_report_position_pdb_longest <- merge(grouped_df_max_report_position_pdb, longest_chain, by = c("uniprotswissprot", "pdb"), all = TRUE)

# # Filter the rows where longest_chain == "Yes"
grouped_df_max_report_position_pdb_longest <-
  grouped_df_max_report_position_pdb_longest |>
  filter(longest_chain == "Yes")
# 
# # unique(longest_chain$uniprotswissprot)
# # unique(grouped_df_max_report_position_pdb$uniprotswissprot)
# # unique(grouped_df_max_report_position_pdb_longest$uniprotswissprot)
# 
# # nrow(longest_chain)
# # nrow(grouped_df_max_report_position_pdb)
# # nrow(grouped_df_max_report_position_pdb_longest)
# 
# # Create a function to plot each protein structure
# plot_protein <- function(pdb_id, Protein_position) {
#   NGLVieweR(pdb_id) |>
#     stageParameters(backgroundColor = "white") |>
#     addRepresentation("cartoon",
#                       param = list(colorScheme = "residueindex")) |>
#     addRepresentation("ball+stick", param = list(
#       colorScheme = "element",
#       colorValue = "red",
#       sele = as.character(Protein_position)  # 'sele' needs to be a character
#     )) |>
#     addRepresentation("label",
#                       param = list(
#                         sele = as.character(Protein_position),
#                         labelType = "format",
#                         labelFormat = "[%(resname)s]%(resno)s",
#                         labelGrouping = "residue",
#                         color = "white",
#                         fontFamiliy = "sans-serif",
#                         xOffset = 1,
#                         yOffset = 0,
#                         zOffset = 0,
#                         fixedSize = FALSE,
#                         radiusType = 2,
#                         radiusSize = 6,
#                         showBackground = TRUE,
#                         backgroundColor="black"
#                         # backgroundOpacity=0.5
#                       )
#     )
# }
# 
# # apply function to data and save html ----
# library(htmlwidgets)
# 
# # note: to get every PDB structure for genes instead of only the longest one, use grouped_df_max_report_position_pdb
# 
# short_test <- grouped_df_max_report_position_pdb_longest |> head(50)
# 
# apply(grouped_df_max_report_position_pdb_longest,
#       1, function(row) {
#   if (row['pdb'] != "") {
#     plot <- plot_protein(row['pdb'], row['Protein_position'])
#     saveWidget(plot, file=paste0(output_dir,row['SYMBOL'], "_",row['uniprotswissprot'], "_", row['pdb'], ".html"))
#   }
# })


# Group data by PDB and create a list of unique Protein_position values for each group
grouped_df <- grouped_df_max_report_position_pdb_longest %>%
  group_by(pdb, SYMBOL, uniprotswissprot) %>%
  summarise(Protein_position = list(unique(Protein_position)))


# Define your plot_protein function to iterate over all positions
plot_protein <- function(pdb_id, Protein_positions) {
  plot <- NGLVieweR(pdb_id) |> stageParameters(backgroundColor = "white") |> 
    addRepresentation("cartoon", param = list(colorScheme = "residueindex"))
  
  for(position in Protein_positions){
    plot <- plot |> 
      addRepresentation("ball+stick", 
                        param = list(colorScheme = "element",
                                     colorValue = "red",
                                     sele = as.character(position)
                        )
      ) |> 
      addRepresentation("label", 
                        param = list(sele = paste0("protein and ", as.character(position)),
                                     labelType = "format",
                                     labelFormat = "[%(resname)s]%(resno)s",
                                     labelGrouping = "residue",
                                     color = "white",
                                     fontFamiliy = "sans-serif",
                                     xOffset = 1,
                                     yOffset = 0,
                                     zOffset = 0,
                                     fixedSize = FALSE,
                                     radiusType = 2,
                                     radiusSize = 6,
                                     showBackground = TRUE,
                                     backgroundColor="black"
                        )
      )
  }
  
  return(plot)
}


library(htmlwidgets)
library(webshot2) # requires chrome on system

# get PNG with Webshot? ("Yes" or "No") - this is significantly slower
use_webshot <- "Yes"
# use_webshot <- "No"

# Apply the function to data and save HTML
apply(grouped_df,
      1, function(row) {
        if (row[['pdb']] != "") {
          cat("PDB: ", row[['pdb']], "\n")  # print pdb id
          cat("Positions: ", unlist(row[['Protein_position']]), "\n")  # print positions
          plot <- plot_protein(row[['pdb']], unlist(row[['Protein_position']]))
          file_html = paste0(output_dir,
                             row[['SYMBOL']], "_",
                             row[['uniprotswissprot']], "_", 
                             row[['pdb']], ".html")
          htmlwidgets::saveWidget(plot, 
                                  file=file_html,
                                  selfcontained = FALSE)
          
          if (use_webshot == "Yes") {
            # Convert the HTML file to a PNG image
            file_png = paste0(output_dir,
                              row[['SYMBOL']], "_",
                              row[['uniprotswissprot']], "_", 
                              row[['pdb']], ".png")
            webshot2::webshot(url = file_html, 
                              file = file_png, 
                              vwidth = 1024, 
                              vheight = 768,
                              expand = -50)  # Negative value to trim borders
          }
        }
      })

# Webshot ----
# install.packages("webshot2")
# remotes::install_github("rstudio/webshot2")

# library(webshot2)

# Combined html and PNG snapshot --- note that it is super slow and memory intensive.

# It generates the filenames for the HTML and PNG files.
# It saves the plot as a standalone HTML file with saveWidget().
# It uses webshot2 to take a screenshot of the saved HTML and saves it as a PNG file.
# It removes the temporary HTML file with unlink().

# apply(grouped_df_max_report_position_pdb_longest,
#       1, function(row) {
#         if (row['pdb'] != "") {
#           plot <- plot_protein(row['pdb'], row['Protein_position'])
#           file_html <- paste0(output_dir,row['SYMBOL'], "_",row['uniprotswissprot'], "_", row['pdb'], ".html")
#           file_png <- paste0(output_dir,row['SYMBOL'], "_",row['uniprotswissprot'], "_", row['pdb'], ".png")
#           saveWidget(plot, file = file_html)
#           webshot2::webshot(url = file_html, file = file_png, vwidth = 1024, vheight = 768)
#           # unlink(file_html)
#         }
#       }
# )
