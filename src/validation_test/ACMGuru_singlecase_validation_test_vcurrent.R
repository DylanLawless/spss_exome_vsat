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

# Other variants in known patient ----

patient_list <- c("BE056", "LU013", "LAU164")

# Create a regular expression from the patient_list
pattern <- paste(patient_list, collapse="|")

# Filter df where 'sample' partially matches any pattern in patient_list
df_known_patients <- df %>%
  filter(grepl(pattern, sample))


df_known_patients$sample  |> unique()




df_report_main_text <- df_known_patients |> 
  dplyr::select(sample, 
                SYMBOL, 
                rownames,
                HGVSp,
                HGVSc,
                Consequence,
                genotype,
                Inheritance,
                gnomAD_AF,
                ACMG_total_score
  ) |> 
  arrange(SYMBOL, sample,
          desc(ACMG_total_score),
          sample)

# Add a new column with the maximum ACMG_total_score for each SYMBOL
df_report_main_text <- df_report_main_text %>%
  group_by(SYMBOL) %>%
  mutate(max_score = max(ACMG_total_score)) %>%
  ungroup()

# Arrange the data frame based on sample, then descending max_score, SYMBOL, and then by desc(ACMG_total_score)
df_report_main_text <- df_report_main_text %>%
  unique() %>%
  arrange(sample, desc(max_score), SYMBOL, desc(ACMG_total_score))

df_report_main_text |> nrow()

colnames(df_report_main_text)[colnames(df_report_main_text) == 'ACMG_total_score'] <- 'ACMG score'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'rownames'] <- 'Variant GRCh38'

# write.csv(df_report_main_text,  paste0("../../data/", output_directory, "ACMGuru_singlecase_genetic_df_report_main_text.csv"))

# are other XL genotype == 2 ? ----

df_x <- df |> filter(Inheritance == "XL") |> filter(genotype == 2) |> 
  dplyr::select(sample, 
                SYMBOL, 
                rownames,
                HGVSp,
                HGVSc,
                Consequence,
                genotype,
                Inheritance,
                gnomAD_AF,
                ACMG_total_score
  ) |> 
  unique()

df_x <- df_x %>%
  group_by(SYMBOL) %>%
  mutate(max_score = max(ACMG_total_score)) %>%
  ungroup()

df_x <- df_x %>%
  unique() %>%
  arrange(sample, desc(max_score), SYMBOL, desc(ACMG_total_score))

df_x$sample |> unique()
df_x |> filter(ACMG_total_score > 5) |> dplyr::select(sample) |> unique()

# It is possible to encounter events were a patient has a high evidence interpretation score in a variant and a lower evidence score in what might be considered a better candidate. 
# For example, in the known PID disease genes we may wish to prioritise a homozygous X-linked variant because of the potentital importance based on "intuition".
# We had 42 cases where a X-linked disease gene was detected with genotype 2 (compount het, homozygous, hemizyous).
# Of these, seven had scores >=6 and all of these are reported in our main text table 1. 
# Therefore, while in our study these were reported as top candidates for the patient, it might always be so. 
# Our analysis passes the first challenge of automated evidence quantification and prioritises interpretation. 
# The remaining hurdle is the countless nuanced criteria to prioritise disease-speicific priors in order to make the final determination. 
# Today, this still often depends on the researchers intuition and both research and commercial software generally provide automated of this final result for "research use only".
# We continue to improve such downstream progress for future work (e.g. https://github.com/DylanLawless/heracles is currently under development). 

# Specific screen for variants that pass QV ----

df_long <- read.table(file = "../../data/ACMGuru_singlecase/ACMGuru_singlecase_genetic_df_report.csv", sep = ",", header = T)

df_borgh <- read.table(file="public_borghesi_Supplementary_Table_4.tsv", header = TRUE, sep = "\t")

df_borgh$SYMBOL <- df_borgh$Gene
df_borgh$HGVSp <- df_borgh$HGVS.p
df_borgh$source <- "Borgh"
df_long$source <- "Lawl"

df_long <- df_long |> dplyr::select(SYMBOL, HGVSp, source)
df_borgh <- df_borgh |> dplyr::select(SYMBOL, HGVSp, source)
df_long$HGVSp <- sub(".*:", "", df_long$HGVSp)
df_m <- merge(df_borgh, df_long, by = c("SYMBOL", "HGVSp")) |> unique()

# Of the variants that pass QV filters, 17 reported in our main text table with new evidence ACMG score.
df_m$gene_variant <- paste(df_m$SYMBOL, df_m$HGVSp, sep=" ")
cat(paste(df_m$gene_variant, collapse=", "), "\n")
