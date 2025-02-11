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
gnomad_freq_thresh_local <- ""

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
	source("../stand_alone_vcf_to_table_no_filters/stand_alone_vcf_to_table_no_filters.R")

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
	# df <- df |> filter(AC < 10)
	
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

print("Guru with no filters returns 6848 variants")


# END HERE ----


gc()


# iuis merge ----
df <- merge(df, iuis, by="SYMBOL", all.x=TRUE) |> dplyr::select(SYMBOL, Inheritance, everything())

# summary ----
# library(Hmisc)
df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
df$AC <- as.numeric(df$AC)
df$AF.x <- as.numeric(df$AF.x)
temp <- df |> ungroup() |> dplyr::select(genotype, Inheritance, IMPACT, Consequence, AF.x, AC, gnomAD_AF, HGVSc) |> unique()


# fix missing chromosome info ----
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


# Validation tests ----
## Borghesi paper validation test ----
# Test for the Borghesi_paper variants present -----
# Filter for any of the known variants:
# WAS p.Glu131Lys, CYBB p.Gly364Arg, CFH p.Pro503Ala

Borghesi_test <- df |> 
  filter(
    (SYMBOL == "WAS" & HGVSp == "ENSP00000365891.4:p.Glu131Lys" & genotype == "2") |
      (SYMBOL == "CYBB" & HGVSp == "ENSP00000367851.4:p.Gly364Arg" & genotype == "2") |
      (SYMBOL == "CFH" & HGVSp == "ENSP00000356399.4:p.Pro503Ala")
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
                ACMG_total_score
  ) |> 
  arrange(SYMBOL,
          sample)

Borghesi_test |> dplyr::select(rownames, SYMBOL, sample, genotype, HGVSp, gnomAD_AF, Inheritance, ACMG_total_score) |> unique()

# There were three "known pathogenic" variants reported; WAS p.Glu131Lys, CYBB p.Gly364Arg, CFH p.Pro503Ala. CFH was filtered out during normal QC stages. The CYBB variant is present for 5 samples and 1 has genotype 2 (hemizygous male on X chromosome), however it is flagged as benign on CliVar and removed by MAF (0.003943) using newer gnomAD database. The WAS variant is present for 4 samples and 1 has genotype 2 (hemizygous male on X chromosome) however it is not identified as pathogenic by ACMG criteria (and ClinVar likely_benign) and and removed by MAF (0.002591) using newer gnomAD database. All other variants are basically do not have sufficient evidence of pathogenicity with updated criteria.

# rm(Borghesi_test)

# Next we wanted to see if other supplemental variants of unknown significance were still present in the input dataset.

# I ran liftover of a large file that has genomic coordinates. That file does not include the amino acid positions as reported in the paper Sup table 4. but it has genes taht seem to match and many more. file: hglft_genome_ae749_8fb250.bed

df_val <- df 
df_borgh_lift <- read.table(file = "hglft_genome_ae749_8fb250.bed")
df_borgh_unlifted <- read.table(file = "all0.01_Borghesi_table_1_171222_larger_variant_slim.tsv", header = TRUE)
df_borgh <- read.table(file="public_borghesi_Supplementary_Table_4.tsv", header = TRUE, sep = "\t")
df_borgh_unlifted$GENE |> unique() |> length()
# This contains 1024 variants in 303 genes (i.e. all variants below AF 0.1 that overlap with the PID gene list at that time).
df_borgh_lift$chr_position <- paste(df_borgh_lift$V1, df_borgh_lift$V2, sep = ":")
df_val$chr_position <- gsub("_.*", "", df_val$rownames)
df_borgh_lift <- df_borgh_lift$chr_position |> unique()
df_validate <- df_val |> filter(chr_position %in% df_borgh_lift) |> unique()
df_validate$rownames |> unique() |> length()
df_validate$rownames |> length()
# This contains 411 variants matches in our cohort

df_validate <- as.data.frame(df_validate)
df_validate <- df_validate %>% 
  tibble::rownames_to_column("rownames")

## Count hits ----
# Generate the bar plot
matches_unfil_total <- df_validate %>%
dplyr::select(rownames, ACMG_total_score) |> 
  nrow()

matches_unfil_uniq <- df_validate %>%
  dplyr::select(rownames, ACMG_total_score) |>
  unique()|>
  nrow()

matches_total <- df_validate %>%
  dplyr::select(rownames, ACMG_total_score) |> 
  filter(ACMG_total_score > 0) |>
  unique()|>
  nrow()

matches_uniq <- df_validate %>%
  dplyr::select(rownames, ACMG_total_score) |> 
  filter(ACMG_total_score > 5) |>
  unique() |>
  nrow()

matches_unfil_total
matches_unfil_uniq
matches_total
matches_uniq

subtitle_text <-  paste0("Total matches = ", matches_unfil_total, 
                  ", unique variants = ", matches_unfil_uniq, 
                  ",\n evidence score >0 = ", matches_total,
                  ", unique variants evidence score >0 = ", matches_uniq
)

subtitle_text

library(ggplot2); theme_set(theme_bw())
p_df_validate <- df_validate %>%
  dplyr::select(rownames, ACMG_total_score) %>%
  group_by(ACMG_total_score) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = ACMG_total_score, y = count)) +
  geom_bar(stat = "identity", fill= "darkblue", color = "black") +
  labs(x = "Total score", y = "Count of variants after new QC/QV\noverlapping with previous study", title = "Count of variants per interpretation score",
# subtitle = subtitle_text
) +
  geom_vline(xintercept = 5.5, linetype = "dotted")

ggsave(paste("../../images/", output_directory, file_suffix, "validation_check_score.pdf", sep = "") ,plot = p_df_validate )

# END -----


df_val <- df_val |> dplyr::select(SYMBOL, HGVSp, rownames)
df_borgh$source <- "borgh"
df_val$source <- "law"


# df_borgh_gene <- df_borgh$SYMBOL |> unique()
# df_validate_genes <- df_val |> filter(SYMBOL %in% df_borgh_gene) |> unique()

head(df_borgh)
head(df_val)



