# source("ACMGuru_gene_illustrate_vcurrent.R")
# this is slow, call from Rds

library(dplyr) 
library(ggplot2)
 
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

# # Install the remotes package if you don't have it
# if (!requireNamespace("remotes", quietly = TRUE)) {
  # install.packages("remotes")
# }
# 
# # Install a specific version of UniprotR
# remotes::install_version("UniprotR", version = "2.4.0")

output_dir <- "../../data/ACMGuru_post_ppi/"
images_dir <- "../../images/ACMGuru_post_ppi_uniprotr/"

geneset_MCL_ID <- c(22, 586, 836)

print("Note here we import grouped_df_max_pathway_id")
grouped_df_max <- readRDS(paste0("../../data/ACMGuru_post_ppi/acmguru_gene_illustrate_grouped_df_max", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

df_report <- readRDS(paste0("../../data/ACMGuru_post_ppi/acmguru_gene_illustrate_df_report", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

df_report |> dplyr::select(SYMBOL) |> unique()
df_report |> dplyr::select(seqid) |> unique()
grouped_df_max|> dplyr::select(seqid) |> unique()

rm(list=setdiff(ls(), c("output_dir", "images_dir", "grouped_df_max", "df_report", "geneset_MCL_ID")))

# add pathway_id to grouped_df_max
grouped_df_max$seqid |> unique()
grouped_df_max$SYMBOL |> unique()

Accessions_full <- df_report |> dplyr::select(SYMBOL, seqid, pathway_id) |> unique()
# Accessions <- df_report |> dplyr::select(seqid) |> unique() |> as.vector()
Accessions <- grouped_df_max$seqid

# tmp <- df_report |> select(SYMBOL, pathway_id, seqid) |> unique()
# tmp2 <- merge(grouped_df_max, tmp, all.x = T)
# grouped_df_max <- tmp2

# Specify file path
# output_dir = "../../data/ACMGuru_post_ppi/"

#Read Accessions from csv file , Note : Accessions must be in the first column. 
# Accessions <-GetAccessionList("Accessions.csv") 
# head(Accessions)

# Accessions <- grouped_df_max |> select(seqid, pathway_id)

# Download ----
#Get Taxonomy Information 
print("TaxaObj <- GetNamesTaxa(Accessions)")
# TaxaObj <- GetNamesTaxa(Accessions)

# Get Gene ontolgy Information
print("GeneOntologyObj <- GetProteinGOInfo(Accessions)")
# GeneOntologyObj <- GetProteinGOInfo(Accessions)

# GetProteinFunction
print("ProteinFunction <- GetProteinFunction(Accessions)")
# ProteinFunction <- GetProteinFunction(Accessions)

# save ----
# saveRDS(TaxaObj, file=paste(output_dir, "ontology_taxa/TaxaObj.Rds", sep = ""))
# saveRDS(GeneOntologyObj, file=paste(output_dir, "ontology_taxa/GeneOntologyObj.Rds", sep = ""))
# saveRDS(ProteinFunction, file=paste(output_dir, "ontology_taxa/ProteinFunction.Rds", sep = ""))

# read local copy ----
TaxaObj <- readRDS(file=paste(output_dir, "ontology_taxa/TaxaObj.Rds", sep = ""))
GeneOntologyObj <- readRDS(file=paste(output_dir, "ontology_taxa/GeneOntologyObj.Rds", sep = ""))
ProteinFunction <- readRDS(file=paste(output_dir, "ontology_taxa/ProteinFunction.Rds", sep = ""))

TaxaObj$seqid <- rownames(TaxaObj)
GeneOntologyObj$seqid <- rownames(GeneOntologyObj)
ProteinFunction$seqid <- rownames(ProteinFunction)

# add pathway id ----

tmp <- df_report |> dplyr::select(SYMBOL, pathway_id, seqid) |> unique()
TaxaObj <- merge(TaxaObj, tmp, all.x = T)
GeneOntologyObj <- merge(GeneOntologyObj, tmp, all.x = T)
ProteinFunction <- merge(ProteinFunction, tmp, all.x = T)

rm(tmp)

TaxaObj <- TaxaObj |> filter(!is.na(pathway_id))
GeneOntologyObj <- GeneOntologyObj |> filter(!is.na(pathway_id))
ProteinFunction <- ProteinFunction |> filter(!is.na(pathway_id))

# Pathway.Enr(Accessions)
# str(Accessions)


# Plotting ----
# Plot Biological process information top 10 go terms
# p_gob <- PlotGOBiological(GeneOntologyObj, Top = 10) 
# p_gom <- Plot.GOMolecular(GeneOntologyObj, Top = 20)
# p_gsc <- Plot.GOSubCellular(GeneOntologyObj) 
# p_goa <- PlotGOAll(GOObj = GeneOntologyObj, Top = 10, directorypath = getwd(), width = 8, height = 5)

geneset_MCL_ID_str <- paste(geneset_MCL_ID, collapse = "_")

#Visualize Chromosomes localization
p_chr <- PlotChromosomeInfo( TaxaObj )
p_chr
ggsave(paste(images_dir, "uniprotr_p_chr_merged_", geneset_MCL_ID_str, ".pdf", sep = "") , plot = p_chr, width = 12, height = 8 )
  
#Combine Gene ontology plots into one plot 
p_goi <- PlotGoInfo(GeneOntologyObj)
p_goi
ggsave(paste(images_dir, "uniprotr_p_goi_merged_", geneset_MCL_ID_str, ".pdf", sep = "") , plot = p_goi, width = 12, height = 10 )
  
# Enrichment analysis using KEGG, Reactome of protein list
Accessions_match <- Accessions_full |> dplyr::select(seqid) |> unique() |> as.list()
p_kr <- Pathway.Enr(Accessions_match)
p_kr
ggsave(paste(images_dir, "uniprotr_p_kr_merged_", geneset_MCL_ID_str, ".pdf", sep = "") , plot = p_kr, width = 10, height = 6 )

# loop on pathways ----
for (pathway in geneset_MCL_ID) {

print(paste("Running:", pathway))

#Visualize Chromosomes localization
print(paste("chr..."))
p_chr <- PlotChromosomeInfo( (TaxaObj |> filter(pathway_id == pathway)))
# p_chr
ggsave(paste(images_dir, "uniprotr_p_chr_pathway_ID_", pathway, ".pdf", sep = "") , plot = p_chr, width = 15, height = 8 )

#Combine Gene ontology plots into one plot 
print(paste("goi..."))
p_goi <- PlotGoInfo( (GeneOntologyObj |> filter(pathway_id == pathway)))
# p_goi
ggsave(paste(images_dir, "uniprotr_p_goi_pathway_ID_", pathway, ".pdf", sep = "") , plot = p_goi, width = 18, height = 8 )

# Enrichment analysis using KEGG, Reactome of protein list
print(paste("filt..."))
Accessions_match <- Accessions_full |> filter(pathway_id == pathway) |> dplyr::select(seqid) |> unique() |> as.list()
print(paste("kr..."))
p_kr <- Pathway.Enr(Accessions_match)
# p_kr
ggsave(paste(images_dir, "uniprotr_p_kr_pathway_ID_", pathway, ".pdf", sep = "") , plot = p_kr, width = 8, height = 6 )

}

# get discussion ----
names(TaxaObj)
names(GeneOntologyObj)
# GeneOntologyObj[2:10,4]

# seqid used to merege which were added earlier
grouped_df_max_GO <- merge(grouped_df_max, GeneOntologyObj)
grouped_df_max_GO_taxa <- merge(grouped_df_max_GO, TaxaObj)
grouped_df_max_GO_taxa_funct <- merge(grouped_df_max_GO_taxa, ProteinFunction)

df_report_grouped_grouped_df_max_GO_taxa_funct <- merge(df_report, grouped_df_max_GO_taxa_funct)

df_report_discussion <-
  df_report_grouped_grouped_df_max_GO_taxa_funct |>
  dplyr::select(pathway_id, SYMBOL,
         Protein.names,
         "Gene.Ontology..molecular.function.",
         "Function..CC."
         ) |>
  unique()

df_report_discussion <-
  df_report_discussion |> arrange(pathway_id, SYMBOL)
  
readr::write_tsv(df_report_discussion, file=(paste0("../../data/ACMGuru_post_ppi/df_report_discussion_", paste(geneset_MCL_ID, collapse="_"), ".tsv"))
          )

write.csv(df_report_discussion, file=(paste0("../../data/ACMGuru_post_ppi/df_report_discussion_", paste(geneset_MCL_ID, collapse="_"), ".csv")),
          row.names = FALSE)

