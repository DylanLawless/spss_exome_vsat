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

file_suffix <- "singlecase_"
output_directory <- "singlecase"

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
chromosomes <- c(1:22, "X", "Y")
chromosomes <- c(1:22, "X")

# Generate file names using paste0 and the chromosome identifiers
file_list <- paste0(
  "../../data/", output_directory, "/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1_chr_", 
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
		summarize(genotype_total_frequency = sum(genotype)/n(), .groups = "drop") %>%
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

df_pathway <- do.call(rbind, df_pathway_list)
df <- df_pathway
df <- df |> filter(!is.na(SYMBOL)) # clean out unassigned
hold <- df

# saveRDS(df, "./df.Rds")
# df <- readRDS("./df.Rds")

# tidy rm trash
rm(list=setdiff(ls(), c("df",  "df_acmg", "df_acmg_caveat", "hold", "iuis", "varsome", "output_directory", "file_suffix")))
gc()
df <- hold

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

  # source acmg filters ----
source("../ACMG_filters/acmg_filters.R")

# plot scores ----
# 
p.acmg_score <- df |> 
	ggplot(aes(x = as.character(ACMG_total_score), fill= as.numeric(ACMG_total_score) )) +
	geom_histogram(stat='count', bins = length(acmg_scores), color="black") +
	theme_minimal() +
	xlab("ACMG score") +
	ylab("No. variants") +
	geom_text(stat='count', aes(label=..count.., y=..count..+300), color = "black") + 
	guides(fill=FALSE) +
	scale_fill_scico(palette = 'bamako', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.acmg_score 
ggsave(paste("../../images/ACMGuru_singlecase_", file_suffix, "acmg_score.pdf", sep = "") ,plot = p.acmg_score )


# panel ----
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (
	(p.criteria_gene_total) / ( p.variants_per_criteria | p.criteria_per_sample ) / ( p.pathogenicity_distributions | p.acmg_score)
) + plot_annotation(tag_levels = 'A')
ggsave(paste("../../images/ACMGuru_singlecase_", file_suffix, "patch1.pdf", sep = "") ,plot = patch1  + plot_annotation(tag_levels = 'A'), width = 8, height = 10 )

patch2 <- (
	(p.criteria_gene_total) / ( p.variants_per_criteria | p.criteria_per_sample ) / ( p.pathogenicity_distributions | p.acmg_score)
)  | (p.pathogenicity_distributions_engines_threshold) + plot_annotation(tag_levels = 'A')
# patch2
ggsave(paste("../../images/ACMGuru_singlecase_", file_suffix, "patch2.pdf", sep = "") ,plot = patch2 + plot_annotation(tag_levels = 'A'), width = 16, height = 10 )
 
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

df_report |> names()

df_report_main_text <- df_report |> 
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
          sample.id)

colnames(df_report_main_text)[colnames(df_report_main_text) == 'Strong_pathogenic_GE'] <- 'Strong_patho'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'Moderate_pathogenic_GE'] <- 'Moder_patho'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'Supporting_pathogenic_GE'] <- 'Suppor_patho'

saveRDS(df_report, file=paste0("../../data/", output_directory, "/df_report.Rds")
)

saveRDS(df_report_main_text, file=paste0("../../data/", output_directory, "/df_report_main_text.Rds"))

geneset_MCL_ID <- "" #ignore pathway level info
write.csv(df_report_main_text,  paste0("../../data/", output_directory, "/ACMGuru_singlecase_genetic_df_report_main_text.csv"))

write.csv(df_report,  paste0("../../data/", output_directory, "/ACMGuru_singlecase_genetic_df_report.csv"))

