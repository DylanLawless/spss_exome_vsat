# ACMGuru ----
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scico) # devtools::install_github("thomasp85/scico") # scico_palette_show()
library(grid)
library(forcats) # new facet labels
library(ggrepel)
library(ggpubr) # For ggarrange
library(cowplot) # For get_legend
library(patchwork) # plots in panels
library(knitr)

file_suffix <- "ACMGuru_singlecase_"
output_directory <- "ACMGuru_singlecase/"

# setting population filter threshold ----
# The same style is applied in stand_alone_vcf_to_table but the threshold is there set to 1 for generality and could instead be passed from this runner script. You can see from the following tests what effect this step has. This is a sanity test to see if our scores are as expected or if gnomad freq is basically the same proxy, which it is not at the final candidate stage. 
gnomad_freq_thresh_local <- ""
# no threshold = 66 candidate variants
# gnomad_rare_thresh <- 1400 # 46 OK
# gnomad_rare_thresh <- 140 # 19 candidate variants ?
# gnomad_rare_thresh <- 1.4 # 43 candidate variants OK
# gnomad_total <- 140000 # approx number in gnomad
# gnomad_rare_thresh <- 140 # approx carriers
# gnomad_freq_thresh_local <- gnomad_rare_thresh/gnomad_total
# gnomad_freq_thresh_local

# acmg ----
# For reference
df_acmg <- fread("../../ref/acmg_criteria_table.txt", sep = "\t", header = TRUE, fill=TRUE)

df_acmg_caveat <- fread("../../ref/acmg_criteria_table_caveats.txt", sep = "\t", header = TRUE)

# iuis ----

iuis <- read.table(
  file = "../../ref/10875_2022_1289_MOESM2_ESM_DLcleaned.tsv",
  sep = "\t",
  fill = TRUE,  # To handle rows with fewer columns
  header = TRUE # Change this based on whether the first line is a header
)

colnames(iuis)[colnames(iuis) == 'Gene.symbol'] <- 'SYMBOL'

# varsome ----
# LE = less than equal to, GE = greater than equal to
varsome <- read.table(file = "../../ref/varsome_calibrated_insilico_thresholds.tsv", sep="\t", header = TRUE)

# qv ----
# Define the chromosome identifiers
chromosomes <- c(1:22, "X")
# chromosomes <- c(21:22, "X") # TEMP TEST
# chromosomes <- c(21:22, "X")

# Generate file names using paste0 and the chromosome identifiers
file_list <- paste0(
  "../../data/", output_directory, "bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1_chr_", 
  chromosomes, 
  ".vcf.gz"
)

df_pathway_list <- list()
# for (f in 21) {
for (f in 1:length(file_list)) {
	cat("Now analysing", f, "\n")
	source("../stand_alone_vcf_to_table/stand_alone_vcf_to_table.R")

	# qv clean ----
	df$cohort_pheno <- df$sample

	# "setpt" = controls "0" / not "setpt" = cases "1"
	df$cohort_pheno[grep("^setpt", df$sample)] <- "0"
	df$cohort_pheno[!grepl("^setpt", df$sample)] <- "1"

	# frequency for cases and controls
	df_genotype_frequency <- df %>%
		dplyr::select(sample, rownames, genotype) %>% 
		unique() %>% # this is import to count genomic positions once rather than transcripts
		mutate(cohort_pheno = ifelse(grepl("^setpt", sample), "0", "1")) %>%
		group_by(rownames, cohort_pheno) %>%
		dplyr::summarize(genotype_total_frequency = sum(genotype)/n(), .groups = "drop") %>%
		pivot_wider(names_from = cohort_pheno, values_from = genotype_total_frequency, names_prefix = "frequency_in_")  %>%
		mutate(is_frequency_in_0_less = ifelse(frequency_in_0 < frequency_in_1, "Yes", "No"))

	df <- df |> filter(genotype > 0) # Keep only variants
	df <- merge(df, df_genotype_frequency, all.x=TRUE)
	rm(df_genotype_frequency)
	
	df <- df |> filter(IMPACT %in% c("HIGH", "MODERATE"))
	
	df <- df |> dplyr::select(-"ClinVar.x",
									  - "ClinVar_CLNSIG.x",
									  - "ClinVar_CLNREVSTAT.x",
									  - "ClinVar_CLNDN.x") # annotation duplicates
	
	df <- df |> distinct()
	df <- df |> filter(cohort_pheno == 1)
	df <- df |> filter(AC < 10)
	
	df_pathway_list[[f]] <- df
}

cat("\nFinished vcf_to_table")

# dfx <- df_pathway_list[[23]] # check chr X

df_pathway <- do.call(rbind, df_pathway_list)
df <- df_pathway
df <- df |> filter(!is.na(SYMBOL)) # clean out unassigned
hold <- df
df <- hold
# filter gnomad

# print(paste("Filtering with gnomad_freq_thresh_local", gnomad_freq_thresh_local))
# df <- df |> filter(gnomAD_AF < gnomad_freq_thresh_local)


# print("!!! PS1 test !!\n\n")
# print(hold$CLIN_SIG |> unique())
# print(hold$ClinVar_CLNDN.y |> unique())
# print(hold |> dplyr::select(CLIN_SIG, ClinVar_CLNDN.y) |> unique())
# hold_x <- hold |> dplyr::select(CLIN_SIG, ClinVar_CLNDN.y) |> unique()
# hold_x <- hold |> filter(genotype == 2)


# TEST BREAK 
# hold <- hold |> filter(CHROM == "chrX")
# hold$CHROM |> unique()

# saveRDS(df, "./df.Rds")
# df <- readRDS("./df.Rds")

# tidy rm trash
rm(list=setdiff(ls(), c("gnomad_freq_thresh_local", "df",  "df_acmg", "df_acmg_caveat", "hold", "iuis", "varsome", "output_directory", "file_suffix")))
gc()


# iuis merge ----
df <- merge(df, iuis, by="SYMBOL", all.x=TRUE) |> dplyr::select(SYMBOL, Inheritance, everything())

# summary ----
# library(Hmisc)
df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
df$AC <- as.numeric(df$AC)
df$AF.x <- as.numeric(df$AF.x)
temp <- df |> ungroup() |> dplyr::select(genotype, Inheritance, IMPACT, Consequence, AF.x, AC, gnomAD_AF, HGVSc) |> unique()

temp |> 
  group_by(genotype) |>
  summarise(n())

temp |> 
  group_by(Inheritance) |>
  summarise(n())

temp |> 
  group_by(IMPACT) |>
  summarise(n())

temp |> 
  group_by(Consequence) |>
  summarise(n())

temp |> 
  group_by(AF.x) |>
  summarise(n())

temp |> 
  group_by(AC) |>
  summarise(n())

temp |>
  ungroup() |>
  dplyr::select(HGVSc) |>
  unique() |>
  summarise(n())

# Test for the Borghesi_paper variants present -----
# Filter for any of the known variants:
# WAS p.Glu131Lys, CYBB p.Gly364Arg, CFH p.Pro503Ala

Borghesi_test <- df |> 
  filter(
    (SYMBOL == "WAS" & HGVSp == "ENSP00000365891.4:p.Glu131Lys") |
      (SYMBOL == "CYBB" & HGVSp == "ENSP00000367851.4:p.Gly364Arg") |
      (SYMBOL == "CFH")
  ) |>
  dplyr::select(sample, 
                SYMBOL, 
                rownames,
                chr,
                HGVSp,
                HGVSc,
                Consequence,
                "IMPACT", 
                "genotype", 
                "Inheritance", 
                "CLIN_SIG", 
                "gnomAD_AF", 
  ) |> 
  arrange(SYMBOL,
          sample)

Borghesi_test |> dplyr::select(SYMBOL, sample, genotype, gnomAD_AF) |> unique()

# There were three "known pathogenic" variants reported; WAS p.Glu131Lys, CYBB p.Gly364Arg, CFH p.Pro503Ala. CFH was filtered out during normal QC stages. The CYBB variant is present for 5 samples and 1 has genotype 2 (hemizygous male on X chromosome), however it is flagged as benign on CliVar and removed by MAF (0.003943) using newer gnomAD database. The WAS variant is present for 4 samples and 1 has genotype 2 (hemizygous male on X chromosome) however it is not identified as pathogenic by ACMG criteria (and nCliVar likely_benign) and and removed by MAF (0.002591) using newer gnomAD database. All other variants are basically do not have sufficient evidence of pathogenicity with updated criteria.

rm(Borghesi_test)

# df_desc <- describe(temp)
# df_desc
# latex(df_desc, file = "./df_desc.tex")
# rm(temp)

# fix missing chromosome info ----
# chr comes from dbnsfp (I believe), but is a sensible header name so we will reuse it
# Remove "chr" prefix from seqname and create a new 'chr' column
df <- df %>%
  mutate(chr = sub("chr", "", seqnames)) 

# comp_het_flag ----
# flag for comp het. WARNING NOT PHASE CHECKED
df <- df %>%
	group_by(sample, SYMBOL) %>%
	mutate(comp_het_flag = ifelse(n() > 1, 1, NA)) 

# same flag for genotype == 2 (homozygous)
df <- df %>%
	mutate(comp_het_flag = ifelse(is.na(comp_het_flag) & genotype == 2, 1, comp_het_flag)) %>%
	ungroup() %>%
	dplyr::select(comp_het_flag, everything())

# Update the flag for chromosome X observations to include a check for genotype == 2
# NB:: this must be screened using sex info
df <- df %>%
  mutate(comp_het_flag = ifelse(chr == "X" & genotype == 2, 1, comp_het_flag)) %>%
  dplyr::select(comp_het_flag, everything())


# print(n = 30, df |> filter(SYMBOL == "WAS") |> dplyr::select(IMPACT, SYMBOL, genotype))

# source acmg filters ----
print("We are adding PS3 now")

source("../ACMG_filters/distribution_variables.R")
source("../ACMG_filters/acmg_filters.R")

# plot scores ----
# 
p.acmg_score <- df |> 
  arrange(ACMG_total_score) |>
  # ggplot(aes(x = as.character(ACMG_total_score), fill= as.numeric(ACMG_total_score) )) +
  ggplot(aes(x = as.factor(ACMG_total_score), fill= as.numeric(ACMG_total_score) )) +
	geom_histogram(stat='count', bins = length(acmg_scores), color="black") +
	theme_minimal() +
	xlab("ACMG score") +
	ylab("No. variants") +
	geom_text(stat='count', aes(label=..count.., y=..count..+300), color = "black") + 
	guides(fill=FALSE) +
	scale_fill_scico(palette = 'bamako', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.acmg_score 

ggsave(paste("../../images/", output_directory, file_suffix, "acmg_score.pdf", sep = "") ,plot = p.acmg_score )


# panel ----
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (
	(p.criteria_gene_total) / ( p.variants_per_criteria | p.criteria_per_sample ) / ( p.pathogenicity_distributions | p.acmg_score)
) + plot_annotation(tag_levels = 'A')
ggsave(paste("../../images/", output_directory, file_suffix, "patch1.pdf", sep = "") ,plot = patch1  + plot_annotation(tag_levels = 'A'), width = 8, height = 10 )

patch2 <- (
	(p.criteria_gene_total) / ( p.variants_per_criteria | p.criteria_per_sample ) / ( p.pathogenicity_distributions | p.acmg_score)
)  | (p.pathogenicity_distributions_engines_threshold) + plot_annotation(tag_levels = 'A')
# patch2
ggsave(paste("../../images/", output_directory, file_suffix, "patch2.pdf", sep = "") ,plot = patch2 + plot_annotation(tag_levels = 'A'), width = 16, height = 10 )
 
# plot order
# p.criteria_count_each_gene
# p.criteria_gene_total
# p.variants_per_criteria
# p.criteria_per_sample
# p.pathogenicity_distributions
# p.pathogenicity_distributions_engines_threshold
# p.acmg_score

# * * Report * *----
# df_report <- df |> filter(ACMG_count > 0)
df_report <- df |> filter(ACMG_total_score > 0)
# see: iuis_iei_table.R for reactable

# clean up the result data for merging
df_report <- df_report |> dplyr::select(sample, everything())

# clean IDs 
df_report <- separate(df_report, sample, into = c("V1", "V2", "V3", "V4", "V5"))

df_report <- df_report |> mutate(V1 = ifelse(V1 == "raw", NA, V1))

df_report <- df_report |>
  unite(V1, V2, col = "sample.id", sep = "", na.rm = TRUE) 

df_report <- df_report |> dplyr::select(-V3, -V4, -V5)

# Collapse Inheritance column if more than one ----
# Aggregating the dataframe to collapse Inheritance while keeping all columns
df_report <- df_report %>%
  group_by(SYMBOL, rownames, sample.id) %>%
  mutate(Inheritance = paste(unique(Inheritance), collapse = " / ")) %>%
  ungroup() %>%
  distinct(SYMBOL, rownames, sample.id, .keep_all = TRUE)

# collect columns where evidence was used
list_of_used_columns <- c()
list_of_used_columns <- c(list_of_used_columns,
                          "IMPACT", 
                          "genotype", 
                          "Inheritance", 
                          "CLIN_SIG", 
                          "gnomAD_AF", 
                          "comp_het_flag",
                          # pathgenicty predics
                          "Strong_pathogenic_GE",
                          "Moderate_pathogenic_GE",
                          "Supporting_pathogenic_GE"
                          # ACMG_PP3 columns set 1:
                          #"CADD_PHRED", "REVEL_rankscore", "MetaLR_pred",
                          #"MutationAssessor_pred", "SIFT_label", "PolyPhen_label"
                          )

# df_report |> names()

print("Updating `ACMG_total_score > 9` to `ACMG_total_score > 9`")
df_report_main_text <- df_report |> 
  filter(ACMG_total_score > 5 ) |>
  filter(ACMG_total_score > 9 ) |>
  dplyr::select(sample.id, 
                # ACMG_score, 
                # ACMG_count, 
                # ACMG_highest, 
                SYMBOL, 
                rownames,
                # chr,
                HGVSp,
                HGVSc,
                Consequence,
                # IMPACT,                   
                genotype,
                Inheritance,
                gnomAD_AF,
                # comp_het_flag  
                # list_of_used_columns
                ACMG_total_score
                ) |> 
  arrange(SYMBOL, sample.id,
          desc(ACMG_total_score),
          sample.id)

colnames(df_report_main_text)[colnames(df_report_main_text) == 'ACMG_total_score'] <- 'ACMG score'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'rownames'] <- 'Variant GRCh38'
# colnames(df_report_main_text)[colnames(df_report_main_text) == 'Strong_pathogenic_GE'] <- 'Strong_patho'
# colnames(df_report_main_text)[colnames(df_report_main_text) == 'Moderate_pathogenic_GE'] <- 'Moder_patho'
# colnames(df_report_main_text)[colnames(df_report_main_text) == 'Supporting_pathogenic_GE'] <- 'Suppor_patho'

saveRDS(df_report, file=paste0("../../data/", output_directory, "df_report.Rds")
)

saveRDS(df_report_main_text, file=paste0("../../data/", output_directory, "df_report_main_text.Rds"))

geneset_MCL_ID <- "" #ignore pathway level info
write.csv(df_report_main_text,  paste0("../../data/", output_directory, "ACMGuru_singlecase_genetic_df_report_main_text.csv"))

write.csv(df_report,  paste0("../../data/", output_directory, "ACMGuru_singlecase_genetic_df_report.csv"))


# Compare to borghesi VUS supplemental table 4 ----
# Define the vector of VUS genes
borghesi_vus_genes <- c("BACH2",
"BCL11B",
"C4A",
"CFH",
"CHD7",
"CXCR4",
"CYBB",
"Gene",
"IL17F",
"IRAK1",
"IRF2BP2",
"KDM6A",
"KMT2D",
"NFAT5",
"NFKB2",
"PLCG2",
"POLA1",
"PTEN",
"SAMD9",
"SEMA3E",
"STAT3",
"TBX1",
"TCF3",
"TNFRSF13B",
"WAS")

# Filter the dataframe for the specified genes
df_borghesi_vus_genes <- df_report %>% 
	dplyr::filter(SYMBOL %in% borghesi_vus_genes) %>% 
	dplyr::filter(ACMG_total_score > 5)

df_borghesi_vus_genes_report_main <- df_borghesi_vus_genes |> 
	filter(ACMG_score > 2 ) |>
	dplyr::select(sample.id, 
								ACMG_total_score,
								# ACMG_score, 
								ACMG_count, 
								ACMG_highest, 
								SYMBOL, 
								rownames,
								chr,
								HGVSp,
								HGVSc,
								Consequence,
								list_of_used_columns
	) |> 
	arrange(SYMBOL,
					desc(ACMG_total_score),
					 sample.id
					) %>%
	distinct()


write.csv(df_borghesi_vus_genes_report_main,  paste0("../../data/", output_directory, "ACMGuru_singlecase_genetic_df_borghesi_vus_genes_report_main.csv"))

# Merge clinical data ----
# Similar steps are done in cohort_summary_post_singlecase.R
# 
# 
# setwd("../ACMGuru_singlecase/")
# # source("ACMGuru_singlecase_vcurrent.R")
# setwd("../cohort_summary_curated")
# 
# # load clinical info
# samples <- read.csv("../../data/cohort_summary_curated/SAMPLE_LIST", header = F)
# samples$sample <- samples$V1
# samples <- samples |> dplyr::select(-V1)
# 
# # Pheno ----
# # Create new column "cohort_pheno"
# samples$cohort_pheno <- samples$sample
# 
# # Replace any value that starts with "setpt" with "0" in the "cohort_pheno" column
# samples$cohort_pheno[grep("^setpt", samples$sample)] <- "0"
# 
# # Replace any value that does not start with "setpt" with "1" in the "cohort_pheno" column
# samples$cohort_pheno[!grepl("^setpt", samples$sample)] <- "1"
# 
# # clean IDs
# samples <- separate(samples, sample, into = c("V1", "V2", "V3", "V4", "V5"))
# 
# samples <- samples |>
# 	mutate(V1 = ifelse(V1 == "raw", NA, V1))
# 
# samples <- samples |>
# 	unite(V1, V2, col = "sample.id", sep = "", na.rm = TRUE)
# 
# samples <- samples |> filter(cohort_pheno == 1)
# 
# # Clinical data
# df <- read.csv("../../data/cohort_summary_curated/sepsis_v2.csv")
# names(df)
# df <- df |> dplyr::select(
# 	-exome_dataset_1,
# 	-exome_dataset_1.1,  
# 	-exome_dataset_2_path,
# 	-exome_dataset_1_path,
# 	-sqlpkey,
# 	-personal.id)
# 
# df$sample.id <- gsub("-", "", df$sample.id)
# 
# df <- merge(samples, df, by = "sample.id", all.x = TRUE)
# 
# missing_samples <- subset(df, is.na(study.site))
# missing_sample_ids <- missing_samples$sample.id
# 
# 
# 
# df_report_main_text_clinical <-
# 	merge(df_report_main_text, df_dedup, by = "sample.id", all.x = TRUE)
# 
# 



# Merge clinical data ----

setwd("../ACMGuru_singlecase/")
# source("ACMGuru_singlecase_vcurrent.R")
setwd("../cohort_summary_curated")

# This step is not required to make a generic phenotype from sample list, only inthe statistical test of PPI
# # load clinical info
# samples <- read.csv("../../data/cohort_summary_curated/SAMPLE_LIST", header = F)
# 
# samples$sample <- samples$V1
# samples <- samples |> dplyr::select(-V1)
# 
# # Pheno ----
# # Create new column "cohort_pheno"
# samples$cohort_pheno <- samples$sample
# 
# # Replace any value that starts with "setpt" with "0" in the "cohort_pheno" column
# samples$cohort_pheno[grep("^setpt", samples$sample)] <- "0"
# 
# # Replace any value that does not start with "setpt" with "1" in the "cohort_pheno" column
# samples$cohort_pheno[!grepl("^setpt", samples$sample)] <- "1"
# 
# # clean IDs
# samples <- separate(samples, sample, into = c("V1", "V2", "V3", "V4", "V5"))
# 
# samples <- samples |>
# 	mutate(V1 = ifelse(V1 == "raw", NA, V1))
# 
# samples <- samples |>
# 	unite(V1, V2, col = "sample.id", sep = "", na.rm = TRUE)
# 
# samples <- samples |> filter(cohort_pheno == 1)

# Clinical data ----
df <- read.csv("../../data/cohort_summary_curated/sepsis_v2.csv")

df <- df |> dplyr::select(
	-exome_dataset_1,
	-exome_dataset_1.1,  
	-exome_dataset_2_path,
	-exome_dataset_1_path,
	-sqlpkey,
	-personal.id)

df$sample.id <- gsub("-", "", df$sample.id)
df_v1 <- df

# Clinical data - most updated?
df <- read.csv("../../data/cohort_summary_curated/20200623_v3_clical_data_for_gwas/spss_gwas_episode.csv")
df$sample.id <- gsub("-", "", df$sample.id)
df_v2 <- df 

# Merge all clinical data ----
# This step meges two datasets and clean up columns which are duplicated as a result.
# Merge the two clinical datasets to get best coverage.

# Merge the data frames while keeping all rows from both
df_v3 <- merge(df_v1, df_v2, by = "sample.id", all = TRUE)

# List of columns that end with .x or .y
columns_x <- grep("\\.x$", names(df_v3), value = TRUE)
columns_y <- grep("\\.y$", names(df_v3), value = TRUE)

# Remove .x or .y suffix and find common base names to identify matches
base_names_x <- sub("\\.x$", "", columns_x)
base_names_y <- sub("\\.y$", "", columns_y)

# For each matched column name, choose which version (.x or .y) to keep
for (base_name in intersect(base_names_x, base_names_y)) {
	# Example logic: keep the .y version if not all are NA, otherwise keep .x
	if (all(is.na(df_v3[[paste0(base_name, ".y")]]))) {
		df_v3[[base_name]] <- df_v3[[paste0(base_name, ".x")]]
	} else {
		df_v3[[base_name]] <- df_v3[[paste0(base_name, ".y")]]
	}
	# Drop the now unnecessary .x and .y columns
	df_v3[[paste0(base_name, ".x")]] <- NULL
	df_v3[[paste0(base_name, ".y")]] <- NULL
}

# Remove the original .x and .y columns if they were not processed (not paired)
remaining_x <- setdiff(columns_x, paste0(base_names_x[base_names_x %in% base_names_y], ".x"))
remaining_y <- setdiff(columns_y, paste0(base_names_y[base_names_y %in% base_names_x], ".y"))
df_v3[remaining_x] <- NULL
df_v3[remaining_y] <- NULL

# Merge clinical with main genetics tables ----
# df <- merge(samples, df, by = "sample.id", all.x = TRUE)

# missing_samples <- subset(df, is.na(study.site))
# missing_sample_ids <- missing_samples$sample.id


df_report_main_text_clinical <-
	merge(df_report_main_text, df, by = "sample.id", all.x = TRUE)

rownames(df_report_main_text_clinical) <- NULL

df_report_main_text_clinical <- df_report_main_text_clinical |>
  ungroup() |>
  arrange(SYMBOL,
          sample.id)

saveRDS(df_report_main_text_clinical, file=paste0("../../data/", output_directory, "df_report_main_text_clinical.Rds"))

geneset_MCL_ID <- "" #ignore pathway level info
write.csv(df_report_main_text_clinical,  paste0("../../data/", output_directory, "ACMGuru_singlecase_genetic_df_report_main_text_clinical.csv"))


df_report_main_text_clinical_short <- df_report_main_text_clinical |>
  dplyr::select(sample.id:'ACMG score',sex,age.at.bc,pathogen.grp, focus.grp, picu)

write.csv(df_report_main_text_clinical_short,  paste0("../../data/", output_directory, "ACMGuru_singlecase_genetic_df_report_main_text_clinical_short.csv"))



