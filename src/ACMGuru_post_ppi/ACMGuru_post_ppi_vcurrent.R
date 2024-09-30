# ACMGuru ----

# https://varsome.com/about/resources/germline-implementation/
# https://mart.ensembl.org/info/genome/variation/prediction/protein_function.html
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scico) # devtools::install_github("thomasp85/scico")
# scico_palette_show()
library(knitr)
library(ggpubr) # For ggarrange
library(cowplot) # For get_legend
library(gridExtra)
library(grid)
library(forcats) # new facet labels
library(ggrepel)
library(patchwork)

# make a loop for all genesets:
# geneset_MCL_ID <- "22"

geneset_MCL_ID <- c(22, 586)
geneset_MCL_ID[[1]]
geneset_MCL_ID[[2]]

# geneset_MCL_ID <- "586"
file_suffix <- paste("ACMGuru_post_ppi_MCL_ID_", paste(geneset_MCL_ID, collapse="_"), "_", sep = "")
output_directory <- "ACMGuru_post_ppi/"

f1 <- paste("../../data/post_ppi/bcftools_gatk_norm_maf01.recode_vep_conda_impact_MCL_", geneset_MCL_ID[[1]], ".vcf.gz", sep = "")

f2 <- paste("../../data/post_ppi/bcftools_gatk_norm_maf01.recode_vep_conda_impact_MCL_", geneset_MCL_ID[[2]], ".vcf.gz", sep = "")

file_list <- c(f1, f2)

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

# for (f in 6) {
# file_list <- c(
# 	paste0("../../data/ACMGuru_post_ppi/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1_chr_", 1:22, ".vcf.gz")
# 	# ,"../data/annotation/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad_chr_X.vcf.gz", 
# 	#  "../data/annotation/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad_chr_Y.vcf.gz"
# )

df_pathway_list <- list()
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

df_pathway <- do.call(rbind, df_pathway_list)
df <- df_pathway
df <- df |> filter(!is.na(SYMBOL)) # clean out unassigned
hold <- df

# saveRDS(df, "./df.Rds")
# df <- readRDS("./df.Rds")

rm(list=setdiff(ls(), c("df",  "df_acmg", "df_acmg_caveat", "geneset_MCL_ID", "file_suffix", "hold", "iuis", "varsome", "file_suffix", "output_directory")))
gc()
df <- hold

# iuis merge ----
df <- merge(df, iuis, by="SYMBOL", all.x=TRUE) |> dplyr::select(SYMBOL, Inheritance, everything())

# summary ----
# library(Hmisc)
df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
df$AC <- as.numeric(df$AC)
df$AF.x <- as.numeric(df$AF.x)
df_summaries <- df |> ungroup() |> dplyr::select(genotype, Inheritance, IMPACT, Consequence, AF.x, AC, gnomAD_AF, HGVSc) |> unique()

# Number of total variants:
# 857 heterozygous, 10 homozygous.

df |> 
  group_by(genotype) |>
  summarise(variants = n()) 

# Number of unique variants:
# 440 heterozygous, 10 homozygous.
df_summaries |> 
  group_by(genotype) |>
  summarise(variants = n()) 
# |> kable("latex", booktabs = TRUE)


# Inheritance pattern of known IEI:
# AR 7, NA 443
df_summaries |> 
  group_by(Inheritance) |>
  summarise(n())
# |> kable("latex", booktabs = TRUE)

df_summaries |> 
  ungroup() |>
  group_by(Consequence, IMPACT, genotype) |>
  summarise(unique_variants = n()) |>
  arrange(IMPACT, unique_variants, genotype) 
# |>  kable("latex", booktabs = TRUE)

# Calculate the tally per group for "IMPACT":
#   IMPACT   genotype unique_variants

df_summaries %>%
  ungroup() %>%
  group_by(IMPACT, genotype) %>%
  summarise(unique_variants = n(), .groups = 'drop') %>%
  arrange(IMPACT, desc(unique_variants), genotype)

# 42 high impact variants frameshift, stop gained or lost, start lost, splice acceptor or donor), 398 Moderate variants (missense, inframe deletion), and 10 homozygous moderate.

df_summaries %>%
  ungroup() %>%
  group_by(IMPACT, Consequence) %>%
  dplyr::select(Consequence, IMPACT) |>
  unique()

# df_summaries |> 
#   group_by(Consequence, IMPACT, genotype) |>
#   summarise( "unique variants" = n_distinct(HGVSc)) |>
#   arrange(IMPACT, "unique variants", genotype) |>
#   kable("latex", booktabs = TRUE)

# In the enriched PPI we identified 443 unique variants, some of which are present in more than one patient for a total of 867 variant observations.
# The number 
#  1   249

df_summaries |> 
  group_by(AC) |>
  summarise(count = n())
# |> kable("latex", booktabs = TRUE)

df_summaries |> 
  group_by(AC, genotype) |>
  summarise(count = n()) %>%
  arrange(desc(genotype), AC, count)


df_summaries |> 
  group_by(AC) |>
  summarise(count = n()) |>
  filter(AC > 1) |>
  summarise(total  = sum(count))

df_summaries |>
  ungroup() |>
  dplyr::select(HGVSc) |>
  unique() |>
  summarise(n())
# |>  kable("latex", booktabs = TRUE)

# plot AC per var
df_summaries$rownames <- rownames(df_summaries)

# rank by AC
df_summaries_grouped <- df_summaries |>
  dplyr::select(HGVSc, AC, genotype) |>
  arrange(AC) |>
  mutate(Rank = row_number())


ac_count_per_var <- df_summaries_grouped |>
  ggplot(aes(x = Rank, y = AC, fill=as.factor(genotype) )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey", "#ee5d6c"), name = "Carrier\ngenotype", 
                    guide = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  xlab("Unique variant\n(arranged by allele count)") +
  ylab("Alelle count")

ac_count_per_var

ggsave(paste("../../images/", output_directory, file_suffix, "ac_count_per_var.pdf", sep = "") ,plot = ac_count_per_var )

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
p.acmg_score <- df |> 
	ggplot(aes(x = as.character(ACMG_total_score), fill= as.numeric(ACMG_total_score) )) +
	geom_histogram(stat='count', bins = length(acmg_scores), color="black") +
	theme_minimal() +
	xlab("ACMG score") +
	ylab("No. variants") +
	geom_text(stat='count', aes(label=..count.., y=..count..+50), color = "black") + 
	guides(fill=FALSE) +
	scale_fill_scico(palette = 'bamako', direction = -1) # batlowK, acton, lajolla, lapaz, turku
p.acmg_score 
ggsave(paste("../../images/", output_directory, file_suffix, "acmg_score.pdf", sep = "") ,plot = p.acmg_score )


# panel ----
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (
	(p.criteria_gene_total) / ( p.variants_per_criteria | p.criteria_per_sample ) / ( p.pathogenicity_distributions | p.acmg_score)
)  | (p.pathogenicity_distributions_engines_threshold) + plot_annotation(tag_levels = 'A')
patch1
ggsave(paste("../../images/", output_directory, file_suffix, "patch1.pdf", sep = "") ,plot = patch1 + plot_annotation(tag_levels = 'A'), width = 16, height = 10 )
 
# plot order
# p.criteria_count_each_gene
# p.criteria_gene_total
# p.variants_per_criteria
# p.criteria_per_sample
# p.pathogenicity_distributions
# p.pathogenicity_distributions_engines_threshold
# p.acmg_score


# For pathways, summarise set ----
# p.var_per_gene <- 
df_summary_vpg <- df |> 
  dplyr::select(SYMBOL, rownames) |> unique() |>
  group_by(SYMBOL) |> 
  summarise(var_per_gene=n())


 df_summary_nc <- df |> 
  dplyr::select(SYMBOL, rownames, sample) |> unique() |>
  group_by(SYMBOL, rownames) |> 
  summarise(n_carriers=n())

df_summary_unq <- df |> 
	dplyr::select(ACMG_score, ACMG_highest, SYMBOL, rownames, HGVSc, HGVSp) |> unique() |>
	unique()

temp <- merge(df_summary_vpg, df_summary_nc, by="SYMBOL")
df_summary_unq_vpg_nc <- merge(temp, df_summary_unq)
rm(temp)

var_per_gene <- df_summary_unq_vpg_nc |>
  dplyr::select(SYMBOL, var_per_gene) |>
  unique() |>
  ggplot(aes(x=SYMBOL, y=var_per_gene)) +
  geom_point(aes(fill=var_per_gene), color="black", shape = 21) +
  xlab("Gene symbol") +
  ylab("No. variants") +
  theme_minimal()  +
  scale_fill_scico(palette = 'lapaz', direction = 1,
                   name = "Variants\nper gene",) + # batlowK, acton, lajolla, lapaz, turku
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

var_per_gene

ggsave(paste("../../images/", output_directory, file_suffix, "var_per_gene.pdf", sep = "") ,plot = var_per_gene + plot_annotation(tag_levels = 'A'), width = 9, height = 3 )

# joint figure ----
var_per_gene_ac_count_per_var <- (ac_count_per_var / var_per_gene)

ggsave(paste("../../images/", output_directory, file_suffix, "var_per_gene_ac_count_per_var.pdf", sep = "") ,plot = var_per_gene_ac_count_per_var + plot_annotation(tag_levels = 'A'), width = 8, height = 5 )

# Report ----
df_report <- df
# df_report <- df |> filter(ACMG_count > 0)
# df_report <- df |> filter(ACMG_score > 0)
# see: iuis_iei_table.R for reactable

df_report |> dplyr::select(ACMG_total_score, ACMG_count, ACMG_highest, SYMBOL, rownames, Protein_position, CDS_position) |> arrange(desc(ACMG_count))

# t <- df_report |> filter(gnomAD_AF < 1e-4) 

# save a copy to archipelago
df$chr
df$rownames
df$MCL_ID <-  paste(geneset_MCL_ID, collapse="_")

df_archi <- df |> dplyr::select(rownames, MCL_ID)

df_archi <- df_archi %>%
	separate(rownames, into = c("CHR", "BP_variant"), sep = ":") |>
	separate(BP_variant, into = c("BP", "variant"), sep = "_") |>
	mutate(CHR = str_replace(CHR, "chr", ""))

saveRDS(df_archi, paste0("../../data/archipelago/archipelago", paste(geneset_MCL_ID, collapse="_"), ".R"))

# clean up the VSAT result data for merging
df_report_sample_vsat <- df_report |> dplyr::select(sample, everything())

# clean IDs 
df_report_sample_vsat <- separate(df_report_sample_vsat, sample, into = c("V1", "V2", "V3", "V4", "V5"))

df_report_sample_vsat <- df_report_sample_vsat |> mutate(V1 = ifelse(V1 == "raw", NA, V1))

df_report_sample_vsat <- df_report_sample_vsat |>
  unite(V1, V2, col = "sample.id", sep = "", na.rm = TRUE)

df_report_sample_vsat <- df_report_sample_vsat |> filter(cohort_pheno == 1)
df_report_sample_vsat <- df_report_sample_vsat |> dplyr::select(sample.id)
df_report_sample_vsat <- df_report_sample_vsat |> unique()
df_report_sample_vsat$group <- "VSAT_contributer"

# saveRDS(df_report_sample_vsat, file = "../../data/ACMGuru_post_ppi/df_report_sample_vsat.Rds")

saveRDS(df_report_sample_vsat, paste0("../../data/ACMGuru_post_ppi/df_report_sample_vsat_", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

# * * Report * *----
# df_report <- df |> filter(ACMG_count > 0)
# df_report <- df |> filter(ACMG_total_score > 2)
df_report <- df |> filter(ACMG_total_score >= 0)
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
  # filter(ACMG_score > 2 ) |>
  dplyr::select(sample.id, 
                ACMG_total_score,
                ACMG_count, 
                ACMG_highest, 
                rownames, 
                CHROM, REF, ALT, 
                POS, start, end, width, 
                Gene, SYMBOL, HGNC_ID, 
                HGVSp, HGVSc, Consequence,  
                IMPACT, genotype,
                Feature_type, Feature, BIOTYPE, VARIANT_CLASS, CANONICAL,
                list_of_used_columns
  ) |> 
  arrange(SYMBOL,
          desc(ACMG_total_score),
          sample.id)

colnames(df_report_main_text)[colnames(df_report_main_text) == 'Strong_pathogenic_GE'] <- 'Strong_patho'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'Moderate_pathogenic_GE'] <- 'Moder_patho'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'Supporting_pathogenic_GE'] <- 'Suppor_patho'

# saveRDS(df_report, file="../../data/singlecase/df_report.Rds")

saveRDS(df_report_main_text,  paste0("../../data/ACMGuru_post_ppi/df_report_main_text_", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

# collapse Inheritance column 
# df_report_main_text_inhrt <- 
#   df_report_main_text %>%
#   # group_by(across(-c(Inheritance, sample.id))) %>%
#   group_by(across(-c(Inheritance))) %>%
#   summarise(Inheritance = paste(Inheritance, collapse = ", "), .groups = 'drop')

write.csv(df_report_main_text,  paste0("../../data/ACMGuru_post_ppi/ACMGuru_post_ppi_genetic_df_report_main_text_", paste(geneset_MCL_ID, collapse="_"), ".csv"))

write.csv(df_report,  paste0("../../data/ACMGuru_post_ppi/ACMGuru_post_ppi_genetic_df_report_", paste(geneset_MCL_ID, collapse="_"), ".csv"))

df_report |> 
  filter(ACMG_total_score > 1) |> 
  dplyr::select(SYMBOL) |> unique()

# Go now to cohort_summary_curated_r/cohort_summary_post_ppi.R

