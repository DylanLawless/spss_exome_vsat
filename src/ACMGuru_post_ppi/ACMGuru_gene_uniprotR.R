 # source("ACMGuru_gene_illustrate_vcurrent.R")
# this is slow, call from Rds

# uniprotR ----
# install.packages('UniprotR')
library(UniprotR)
# install.packages("UniprotR")

# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# 
# BiocManager::install("UniprotR", force = TRUE)

# The function `PlotEnrichedPathways` is missing from CRAN version so use github version
# if (!require(devtools)) install.packages("devtools")
# library(devtools)
# install_github("Proteomicslab57357/UniprotR")

geneset_MCL_ID <- c(22, 586)
grouped_df_max <- readRDS(paste0("../../data/ACMGuru_post_ppi/acmguru_gene_illustrate_grouped_df_max", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

df_report <- readRDS(paste0("../../data/ACMGuru_post_ppi/acmguru_gene_illustrate_df_report", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

df_report |> dplyr::select(SYMBOL) |> unique()
df_report |> dplyr::select(seqid) |> unique()
grouped_df_max|> dplyr::select(seqid) |> unique()

rm(list=setdiff(ls(), c("grouped_df_max", "df_report", "geneset_MCL_ID")))


# Specify file path
# filename = "your_file_name.jpeg"
output_dir = "../../data/ACMGuru_post_ppi/"
# output_file_path = paste(output_dir, filename, sep = "")

#Read Accessions from csv file , Note : Accessions must be in the first column. 
# Accessions <-GetAccessionList("Accessions.csv") 
# head(Accessions)

Accessions <- grouped_df_max$seqid

#Get Taxonomy Information 
# TaxaObj <- GetNamesTaxa(Accessions) 
#Get Gene ontolgy Information 
# GeneOntologyObj <- GetProteinGOInfo(Accessions) 
# GetProteinFunction
# ProteinFunction <- GetProteinFunction(Accessions)
# saveRDS(TaxaObj, file=paste(output_dir, "ontology_taxa/TaxaObj.Rds", sep = ""))
# saveRDS(GeneOntologyObj, file=paste(output_dir, "ontology_taxa/GeneOntologyObj.Rds", sep = ""))
# saveRDS(ProteinFunction, file=paste(output_dir, "ontology_taxa/ProteinFunction.Rds", sep = ""))

# load local copy
TaxaObj <- readRDS(file=paste(output_dir, "ontology_taxa/TaxaObj.Rds", sep = ""))
GeneOntologyObj <- readRDS(file=paste(output_dir, "ontology_taxa/GeneOntologyObj.Rds", sep = ""))
ProteinFunction <- readRDS(file=paste(output_dir, "ontology_taxa/ProteinFunction.Rds", sep = ""))

# loop ----
# # List of plotting functions
# plot_functions <- list(
#   PlotChromosomeInfo = function() PlotChromosomeInfo(TaxaObj),
#   PlotGOBiological = function() PlotGOBiological(GeneOntologyObj, Top = 10),
#   PlotGOMolecular = function() Plot.GOMolecular(GeneOntologyObj, Top = 20),
#   PlotGOSubCellular = function() Plot.GOSubCellular(GeneOntologyObj),
#   PlotGoInfo = function() PlotGoInfo(GeneOntologyObj),
#   PlotGOAll = function() PlotGOAll(GOObj = GeneOntologyObj, Top = 10, directorypath = getwd(), width = 8, height = 5),
#   PathwayEnr = function() Pathway.Enr(Accessions),
#   PlotEnrichedGO = function() PlotEnrichedGO(Accs = Accessions, Path = getwd(), theme = "lancet", width = 9, height = 5),
#   PlotEnrichedPathways = function() PlotEnrichedPathways(Accs = Accessions, Path = getwd(), theme = "jama", w = 9, h = 5)
# )
# 
# # List of output file names
# file_names <- c("chromosome.jpeg", "biological.jpeg", "molecular.jpeg", "subcellular.jpeg", "info.jpeg", "all.jpeg", "pathway.jpeg", "enrichedGO.jpeg", "enrichedPathways.jpeg")
# 
# # Output directory
# output_dir <- "../../data/ACMGuru_post_ppi/uniprotr/"
# 
# # Loop through the plotting functions and save each plot
# for(i in seq_along(plot_functions)) {
#   output_file_path <- paste(output_dir, file_names[i], sep = "")
#   
#   # Open a JPEG graphics device
#   jpeg(output_file_path)
#   
#   # Generate the plot
#   plot_functions[[i]]()
#   
#   # Close the graphics device
#   dev.off()
# }

#Visualize Chromosomes localization
p_chr <- PlotChromosomeInfo(TaxaObj)

#Plot Biological process information top 10 go terms  
# p_gob <- PlotGOBiological(GeneOntologyObj, Top = 10) 
# p_gom <- Plot.GOMolecular(GeneOntologyObj, Top = 20)
# p_gsc <- Plot.GOSubCellular(GeneOntologyObj) 
# p_goa <- PlotGOAll(GOObj = GeneOntologyObj, Top = 10, directorypath = getwd(), width = 8, height = 5)

#Combine Gene ontology plots into one plot 
p_goi <- PlotGoInfo(GeneOntologyObj)

# Enrichment analysis using KEGG, Reactome of protein list
p_kr <- Pathway.Enr(Accessions)

#For ready graphs for publications 
# Enrichment analysis using KEGG, Reactome of protein list
PlotEnrichedGO(Accs = Accessions, Path = output_dir, theme = "lancet", width = 9, height = 5)
p_erp <- PlotEnrichedPathways(Accs = Accessions, Path = output_dir, theme = "jama", w = 9, h = 5)
# UniprotR::PlotEnrichedPathways()
# UniprotR::PlotEnrichedGO()

# patchwork ----
library(patchwork)
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (p_kr / p_erp) + plot_annotation(tag_levels = 'A')

# Output directory
output_dir <- "../../data/ACMGuru_post_ppi/"

# This one plot has everything
ggsave(paste(output_dir, "uniprotr_combined_go_plots_", output_ID, ".pdf", sep = "") , plot = p_goi, width = 16, height = 10 )

ggsave(paste(output_dir, "uniprotr_combined_keeg_reactome_plots_", output_ID, ".pdf", sep = "") , plot = patch1, width = 10, height = 8 )


# These are individual plots in different combination
# ggsave(paste(output_dir, "p1.pdf", sep = "") , plot = patch1, width = 16, height = 6 )
# ggsave(paste(output_dir, "p2.pdf", sep = "") , plot = patch2, width = 16, height = 6 )
# ggsave(paste(output_dir, "p4.pdf", sep = "") , plot = p_goa, width = 16, height = 10 )


# get discussion ----
names(TaxaObj)
names(GeneOntologyObj)
# GeneOntologyObj[2:10,4]

GeneOntologyObj$seqid  <- rownames(GeneOntologyObj)
grouped_df_max_GO <- merge(grouped_df_max, GeneOntologyObj)

TaxaObj$seqid  <- rownames(TaxaObj)
grouped_df_max_GO_taxa <- merge(grouped_df_max_GO, TaxaObj)

ProteinFunction$seqid  <- rownames(ProteinFunction)
grouped_df_max_GO_taxa_funct <- merge(grouped_df_max_GO_taxa, ProteinFunction)

df_report_grouped_grouped_df_max_GO_taxa_funct <- merge(df_report, grouped_df_max_GO_taxa_funct)

# names(ProteinFunction)

df_report_discussion <-
  df_report_grouped_grouped_df_max_GO_taxa_funct |>
  dplyr::select(SYMBOL,
         Protein.names,
         "Gene.Ontology..molecular.function.",
         "Function..CC."
         ) |>
  unique()

readr::write_tsv(df_report_discussion, file=(paste0("../../data/ACMGuru_post_ppi/df_report_discussion_", paste(geneset_MCL_ID, collapse="_"), ".tsv"))
          )

write.csv(df_report_discussion, file=(paste0("../../data/ACMGuru_post_ppi/df_report_discussion_", paste(geneset_MCL_ID, collapse="_"), ".csv")),
          row.names = FALSE)

