# AMCGuru ----
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
rm(list=setdiff(ls(), c("df",  "df_acmg", "df_acmg_caveat", "file_suffix", "hold", "iuis", "varsome", "output_directory", "file_suffix")))
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

# # PVS1 ----
# # PVS1 are null variants where IMPACT=="HIGH" and inheritance match, in gene where LoF cause disease.
# df$ACMG_PVS1 <- NA
# df <- df %>% dplyr::select(ACMG_PVS1, everything())
# df$ACMG_PVS1 <- ifelse(df$IMPACT == "HIGH" & df$genotype == 2, "PVS1", NA) # homozygous
# df$ACMG_PVS1 <- ifelse(df$IMPACT == "HIGH" & df$Inheritance == "AD", "PVS1", df$ACMG_PVS1) # dominant
# # df |> filter(ACMG_PVS1 == "PVS1")
# 
# # include comp_het if both HIGH impact. WARNING NOT PHASE CHECKED
# df <- df |>
# 	group_by(sample, SYMBOL) |>
# 	mutate(ACMG_PVS1 = ifelse(ACMG_PVS1 == "PVS1", "PVS1", 
# 									  ifelse(sum(IMPACT == "HIGH" & comp_het_flag == 1) >= 2 & IMPACT == "HIGH", "PVS1", ACMG_PVS1))) %>%
# 	ungroup() 
# # df |> filter(ACMG_PVS1 == "PVS1")
# 
# # PS1 ----
# # PS1 Same amino acid change as a previously established pathogenic variant regardless of nucleotide change. Note to keep splicing variant as PSV1 (these are covered by IPACT HIGH).
# df$ACMG_PS1 <- NA
# df <- df %>% dplyr::select(ACMG_PS1, everything())
# df$ACMG_PS1 <- ifelse(df$CLIN_SIG == "pathogenic", "PS1", NA)
# # df |> filter(ACMG_PS1 == "PS1")
# 
# # PS2 skip ----
# # PS2 De novo (both maternity and paternity confirmed) in a patient with the disease and no family history
# # Skip due to no parental genetics.
# 
# # PS3 skip ----
# # PS3 Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.
# # df$ACMG <- ifelse(uniprot == "pathogenic" |
# # 							pubmed == "pathogenic",
# # 						"PS3")
# 
# # PS4 skip ----
# # The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls
# # Skip do to statistical analysis separately
# 
# # PS5 ----
# # The user has additional (value) strong pathogenic evidence
# df$ACMG_PS5 <- NA
# df <- df |> dplyr::select(ACMG_PS5, everything())
# 
# # comp_het with at least 1 HIGH impact. WARNING NOT PHASE CHECKED
# df <- df %>%
# 	group_by(sample, SYMBOL) %>%
# 	mutate(ACMG_PS5 = ifelse(any(IMPACT == "HIGH") & (n() > 1), "PS5", ACMG_PS5)) %>%
# 	ungroup()
# # df |> filter(ACMG_PS5 == "PS5")
# 
# # PM2 ----
# # Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium
# 
# df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
# gnomad_max <- 1e-6 # round down to approx less than 1 on gnomad.
# df$ACMG_PM2 <- NA
# df <- df %>% dplyr::select(ACMG_PM2, everything())
# df$ACMG_PM2 <- ifelse(df$gnomAD_AF < gnomad_max, "PM2", NA)
# # df |> filter(ACMG_PM2 == "PM2")
# 
# # PM3 ----
# # For recessive disorders, detected in trans with a pathogenic variant
# # some redundancy with our PS5 since our rare disease cohort filtering call IMPACT==HIGH equates pathogenic
# df$ACMG_PM3 <- NA
# df <- df %>% dplyr::select(ACMG_PM3, everything())
# df <- df %>%
# 	group_by(sample, SYMBOL) %>%
# 	mutate(ACMG_PM3 = ifelse(comp_het_flag == 1 & (ACMG_PS1 == "PS1" | ACMG_PS5 == "PS5"), 
# 									 "PM3", ACMG_PM3)) %>%
# 	ungroup()
# # df |> filter(ACMG_PM3 == "PM3")
# 
# # PP3 in silico ----
# # Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.). CUSTOM: assigned if >=3 thresholds passed. 
# 
# # In-Silico Predictions: VarSome now implements the ClinGen recommendations from Evidence-based calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for clinical use of PP3/BP4 criteria: # Only one engine at a time is used, depending on availability of data, in order: MitoTip & MitImpact, MetaRNN, CADD (Premium only), DANN (if CADD is not available). # The maximum strength allowed for rules PP3 & BP4 is Strong, even if there may be evidence for Very Strong, with the exception of variants that are predicted splicing (ie: similar to PVS1). # The strength is limited to Supporting, if there's Moderate evidence from rules PM1 or PM5. # Splice prediction (scSNV) is given priority over the other in-silico predictions. # conservation is used for some low-sensitivity variant types, or if no other in-silico prediction is available. Please refer to PP3 and BP4 for more specific detail.
# 
# df <- df |> separate(SIFT, into = c("SIFT_label", "SIFT_score"), sep = "\\(", remove = TRUE) |>
# 	mutate(SIFT_score = str_replace(SIFT_score, "\\)", "")) 
# df <- df |> separate(PolyPhen, into = c("PolyPhen_label", "PolyPhen_score"), sep = "\\(", remove = TRUE) |>
# 	mutate(PolyPhen_score = str_replace(PolyPhen_score, "\\)", "")) 
# df$CADD_PHRED <- as.numeric(df$CADD_PHRED)
# df$REVEL_rankscore <- as.numeric(df$REVEL_rankscore)
# 
# # Define your conditions
# cond_CADD_PHRED <-            df$CADD_PHRED >= 30
# cond_REVEL_rankscore <-       df$REVEL_rankscore > .5
# cond_MetaLR_pred <-           df$MetaLR_pred == "D"
# cond_MutationAssessor_pred <- df$MutationAssessor_pred == "H"
# cond_SIFT_label <-            df$SIFT_label == "deleterious"
# cond_PolyPhen_label <-        df$PolyPhen_label == "probably_damaging"
# 
# # Initialize the ACMG_PP3 column with NA
# df$ACMG_PP3 <- NA
# df <- df %>% dplyr::select(ACMG_PP3, everything())
# 
# # Count the points and store them in ACMG_PP3
# df$ACMG_PP3_count <- rowSums(cbind(cond_CADD_PHRED, cond_REVEL_rankscore, cond_MetaLR_pred, 
# 									  cond_MutationAssessor_pred, cond_SIFT_label, cond_PolyPhen_label), na.rm = TRUE)
# 
# threshold <- 3
# df$ACMG_PP3 <- ifelse(df$ACMG_PP3_count >= threshold, "PP3", NA)
# df |> filter(ACMG_PP3 == "PP3")
# 
# # Remove temporary column
# df$ACMG_PP3_count <- NULL
# 
# rm(list=setdiff(ls(), c("df",  "df_acmg", "df_acmg_caveat", "file_suffix", "hold", "iuis", "varsome")))
# 
# # PP3 In silico: varsome ----
# # Varsome conditions
# varsome |>  dplyr::select(Engine) 
# 
# # Rename varsome to match our data
# names_to_replace <- list(
# 	c("BayesDel_addAF", "BayesDel_addAF_score"),
# 	c("BayesDel_noAF", "BayesDel_noAF_score"),
# 	c("CADD", "CADD_PHRED"),
# 	c("DANN", "DANN_score"),
# 	c("EIGEN", "Eigen.raw_coding"),
# 	c("EIGEN-PC", "Eigen.PC.phred_coding"),
# 	c("FATHMM", "FATHMM_score"),
# 	c("FATHMM-MKL", "fathmm.MKL_coding_score"),
# 	c("FATHMM-XF", "fathmm.XF_coding_score"),
# 	c("LRT", "LRT_score"),
# 	c("M-CAP", "M.CAP_score"),
# 	c("MetaLR", "MetaLR_score"),
# 	c("MetaSVM", "MetaSVM_score"),
# 	c("MetaRNN", "MetaRNN_score"),
# 	c("MutPred", "MutPred_score"),
# 	c("MutationAssessor", "MutationAssessor_score"),
# 	c("MutationTaster", "MutationTaster_score"),
# 	c("phastCons100way_vertebrate", "phastCons100way_vertebrate"),
# 	c("Polyphen2-HDIV", "Polyphen2_HDIV_score"),
# 	c("Polyphen2-HVAR", "Polyphen2_HVAR_score"),
# 	c("PROVEAN", "PROVEAN_score"),
# 	c("REVEL", "REVEL_score"),
# 	c("SIFT", "SIFT_score")
# )
# 
# # Loop over the list and replace the old names with the new names
# for (name_pair in names_to_replace) {
# 	varsome$Engine <- replace(varsome$Engine, varsome$Engine == name_pair[1], name_pair[2])
# }
# 
# # Not used: BLOSUM DANN DEOGEN2 EVE LIST-S2 M-CAP MVP MaxEntScan MitImpact MitoTip PrimateAI SIFT4G phyloP (PhyloP100Way) scSNV-ADA scSNV-RF
# 
# # varsome list of thresholds to tally conditions met
# calculate_varsome_score <- function(df, varsome, pathogenic_type) {
# 	varsome_list <- setNames(varsome[[pathogenic_type]], varsome$Engine)
# 	
# 	df[[pathogenic_type]] <- 0
# 	
# 	for (engine in names(varsome_list)) {
# 		if (!(engine %in% names(df))) {
# 			print(paste(engine, "not found in df, skipping..."))
# 			next
# 		}
# 		
# 		if (!is.numeric(df[[engine]])) {
# 			print(paste(engine, "is not numeric, converting..."))
# 			df[[engine]] <- as.numeric(df[[engine]])
# 		}
# 		
# 		condition <- df[[engine]] >= varsome_list[[engine]]
# 		
# 		condition <- tidyr::replace_na(condition, 0)
# 		
# 		print(paste(engine, ":", sum(is.na(condition)), "NAs.",
# 						pathogenic_type, ":", sum(condition, na.rm = TRUE)))
# 		
# 		df[[pathogenic_type]] <- df[[pathogenic_type]] + condition
# 	}
# 	
# 	return(df)
# }
# 
# df <- calculate_varsome_score(df, varsome, "Strong_pathogenic_GE")
# df <- calculate_varsome_score(df, varsome, "Moderate_pathogenic_GE")
# df <- calculate_varsome_score(df, varsome, "Supporting_pathogenic_GE")
# 
# df <- df |> dplyr::select(ends_with("_pathogenic_GE"), everything())
# 
# # distributions and thresholds 
# # library(tidyverse)
# varsome_thresholds <- varsome %>%
# 	dplyr::select(Engine, ends_with("_pathogenic_GE")) %>%
# 	pivot_longer(cols = -Engine,
# 					 names_to = "pathogenicity",
# 					 values_to = "threshold")
# 
# common_cols <- intersect(varsome_thresholds$Engine, names(df))
# 
# df_long <- df |>
# dplyr::select(all_of(common_cols)) |>
# 	pivot_longer(cols = all_of(common_cols),
# 					 names_to = "Engine",
# 					 values_to = "Score")
# 
# 
# # The Engine names are too long for our plot. Named vector where names are new (long) names and values are old (short) names
# name_mapping <- setNames(sapply(names_to_replace, `[[`, 1), sapply(names_to_replace, `[[`, 2))
# df_long$Engine_short <- name_mapping[df_long$Engine]
# 
# p.pathogenicity_distributions_engines <- df_long |>
# 	# ggplot(aes(x = NormScore, fill=..x..)) +
# 	ggplot(aes(x = Score, fill=..x..)) +
# 	geom_histogram(
# 		#color="black"
# 		) +
# 	facet_wrap(~Engine_short, scales = "free") +
# 	theme_minimal() +
# 	xlab("in silico prediction score") +
# 	ylab("No. qualifying variants")+ 
# 	guides(fill=FALSE) +
# 	scale_fill_scico(palette = 'bamako', direction = 1) # batlowK, acton, lajolla, lapaz, turku
# p.pathogenicity_distributions_engines
# ggsave(paste("../../images/singlecase/", file_suffix, "pathogenicity_distributions_engines.pdf", sep = "") ,plot = p.pathogenicity_distributions_engines)
# 
# # Append a suffix to the pathogenicity column in varsome_thresholds
# varsome_thresholds$pathogenicity <- paste0(varsome_thresholds$pathogenicity, "_threshold")
# 
# # Pivot varsome_thresholds to wide format
# varsome_thresholds_wide <- varsome_thresholds |>
# 	pivot_wider(names_from = pathogenicity, values_from = threshold)
# 
# # Join with df_long
# df_long <- left_join(df_long, varsome_thresholds_wide, by = "Engine")
# 
# 
# Strong_pathogenic_GE_threshold <- 3
# Moderate_pathogenic_GE_threshold <- 4
# Supporting_pathogenic_GE_threshold <- 10
# 
# threshold_results <- df |>
#   summarize(
#     total = n(),
#     strong_threshold = Supporting_pathogenic_GE_threshold,
#     strong_pass = sum(Strong_pathogenic_GE >= Strong_pathogenic_GE_threshold, na.rm = TRUE),
#     strong_percent = (strong_pass / total) * 100,
#     
#     total = n(),
#     moderate_threshold = Moderate_pathogenic_GE_threshold,
#     moderate_pass = sum(Moderate_pathogenic_GE >= Moderate_pathogenic_GE_threshold, na.rm = TRUE),
#     moderate_percent = (moderate_pass / total) * 100,
#     
#     total = n(),
#     supporting_threshold = Supporting_pathogenic_GE_threshold,
#     supporting_pass = sum(Supporting_pathogenic_GE >= Supporting_pathogenic_GE_threshold, na.rm = TRUE),
#     supporting_percent = (supporting_pass / total) * 100
#   )
# 
# threshold_results_long <- threshold_results |>
#   pivot_longer(everything(),
#                names_to = "Measurement",
#                values_to = "Value")
# 
# library(knitr)
# threshold_results_long %>%
#   mutate(Value = round(Value, 2)) %>% 
#   kable("latex", booktabs = TRUE)
# 
# # assign PP3 ----
# df$ACMG_PP3 <- ifelse(
#   df$Strong_pathogenic_GE >= Strong_pathogenic_GE_threshold |
#     df$Moderate_pathogenic_GE >= Moderate_pathogenic_GE_threshold |
#     df$Supporting_pathogenic_GE >= Supporting_pathogenic_GE_threshold, 
#   "PP3", 
#   df$ACMG_PP3
# )
# 
# df |> filter(ACMG_PP3 == "PP3")
# 
# 
# # independent fill scales ----
# # Preparing the data
# df_filtered <- df_long %>%
# 	group_by(Engine) %>%
# 	filter(!all(is.na(Score))) %>%
# 	ungroup()
# 
# # Creating a list of plots for each group
# p.list <- lapply(sort(unique(df_filtered$Engine_short)), function(i) {
# 	
# 	df_group <- df_filtered[df_filtered$Engine_short==i, ]
# 	
# 	df_group |>
# 		ggplot(aes(x = Score, fill=..x..)) +
# 		geom_histogram(bins = 30) +
# 		theme_minimal(base_size = 8) +
# 		labs(subtitle =i) +
# 		xlab("") +
# 		ylab("") +
# 		geom_vline(aes(xintercept = Supporting_pathogenic_GE_threshold),
# 					  linetype = "dashed", color = "#eeaf61") +
# 		geom_vline(aes(xintercept = Moderate_pathogenic_GE_threshold),
# 					  linetype = "dashed", color = "#ee5d6c") +
# 		geom_vline(aes(xintercept = Strong_pathogenic_GE_threshold),
# 					  linetype = "dashed", color = "#6a0d83")+ 
# 		guides(fill=FALSE) +
# 		scale_fill_scico(palette = 'bamako', direction = 1)
# 	
# })
# 
# # add legend
# df_empty <- data.frame()
# legend_only_plot <- 
# 	ggplot(df_empty) +
# 	geom_vline(aes(xintercept = 3, color = "#eeaf61")) +
# 	geom_vline(aes(xintercept = 2, color = "#ee5d6c")) +
# 	geom_vline(aes(xintercept = 1, color = "#6a0d83")) + 
# 	scale_color_identity("", 
# 								breaks = c("#eeaf61", "#ee5d6c", "#6a0d83"), 
# 								labels = c("Supporting pathogenic",
# 											  "Moderate pathogenic",
# 											  "Strong pathogenic"), 
# 								guide = "legend") +
# 	theme_void() +
# 	guides(color = guide_legend(reverse = TRUE))
# 
# legend <- get_legend(legend_only_plot)
# n_plots <- length(p.list)
# ncol <- 5
# nrow <- ceiling(n_plots / ncol) 
# plot_list <- c(p.list, 
# 					rep(list(NULL), nrow * ncol - n_plots - 1), 
# 					list(legend))
# 
# # Arrange all plots together
# p.pathogenicity_distributions_engines_threshold <-  
# 	annotate_figure(
# 		ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow),
# 		left = textGrob("No. qualifying variants", rot = 90, vjust = 1 ),
# 		bottom = textGrob("in silico prediction score" )
# 	)
# 
# p.pathogenicity_distributions_engines_threshold
# ggsave(paste("../../images/singlecase/", file_suffix, "pathogenicity_distributions_engines_threshold.pdf", sep = "") ,plot = p.pathogenicity_distributions_engines_threshold)
# 
# # thresholds passed 
# labels <- c( Strong_pathogenic_GE="Strong", Moderate_pathogenic_GE="Moderate", Supporting_pathogenic_GE="Supporting")
# 
# p.pathogenicity_distributions <- df |> 
# 	tidyr::pivot_longer(cols = ends_with("_pathogenic_GE"),
# 							  names_to = "pathogenicity",
# 							  values_to = "varsome_score") |> 
# 	mutate(pathogenicity = fct_relevel(pathogenicity, 
# 												  "Strong_pathogenic_GE", "Moderate_pathogenic_GE", "Supporting_pathogenic_GE")) |> 
# 	ggplot(aes(x = varsome_score, fill=..x..)) +
# 	geom_histogram(binwidth = 1, color="black") +
# 	geom_text_repel(stat='count', color = "black", 
# 						 box.padding = 0.5, max.overlaps = Inf,
# 						 # padding = unit(0.5, "lines"),
# 						 # nudge_y = 0.05,  
# 						 nudge_x = .4, 
# 						 direction = "y",
# 						 aes(label= ifelse(..count.. < 1000, ..count.., ''))
# 						 ) +
# 	facet_grid(pathogenicity ~ ., labeller=labeller(pathogenicity = labels)) +
# 	theme_minimal() +
# 	xlab("Pathogenicity\nthresholds passed") +
# 	ylab("No. variants")+ 
# 	guides(fill=FALSE) +
# 	scale_fill_scico(palette = 'bamako', direction = 1) # batlowK, acton, lajolla, lapaz, turku
# p.pathogenicity_distributions
# ggsave(paste("../../images/singlecase/", file_suffix, "pathogenicity_distributions.pdf", sep = "") ,plot = p.pathogenicity_distributions)
# 
# # acmg tally  ----
# # List of all ACMG labels
# acmg_labels <- c("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5", 
#                  "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6", 
#                  "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4")
# 
# # Check if each ACMG column exists, if not create it and fill with NA
# for (acmg_label in acmg_labels) {
#   if (!acmg_label %in% names(df)) {
#     df[[acmg_label]] <- NA
#   }
# }
# 
# # Then use coalesce to find the first non-NA ACMG label
# df$ACMG_highest <- dplyr::coalesce(!!!df[acmg_labels])
# df <- df %>% dplyr::select(ACMG_highest, everything())
# 
# # Count the number of non-NA values across the columns
# df$ACMG_count <- rowSums(!is.na(df[, acmg_labels ]))
# df <- df %>% dplyr::select(ACMG_count, everything())
# # df$ACMG_count[df$ACMG_count == 0] <- NA
# 
# p.criteria_count_each_gene <- df |> 
#   filter(ACMG_count > 1) |>
#   ggplot(aes(y = ACMG_count, x = SYMBOL)) +
#   geom_point() +
#   theme_minimal() +
#   theme(axis.text.x  = element_text(angle=45, hjust=1, vjust=1)) +
#   xlab("Gene symbol") +
#   ylab("ACMG criteria count (>1)")
# p.criteria_count_each_gene
# ggsave(paste("../../images/singlecase/", file_suffix, "criteria_count_each_gene.pdf", sep = "") ,plot = p.criteria_count_each_gene )
# 
# # as table
# df |> 
#   filter(ACMG_count > 1) |>
#   dplyr::select(sample, SYMBOL, ACMG_count) |>
#   arrange(desc(ACMG_count))
# 
# p.criteria_gene_total <- df %>%
#   group_by(SYMBOL) |>
#   summarise(acmg_count_per_symbol = sum(ACMG_count)) |>
#   na.omit() |>
#   ggplot(aes(x = acmg_count_per_symbol, fill=..x..) ) +
#   geom_histogram(stat="count", binwidth = 1, color="black"
#   ) +
#   theme_minimal() +
#   xlab("No. ACMG criteria (P) variants per gene") +
#   ylab("Number of genes") +
#   geom_text(stat='count', aes(label=..count.., y=..count..+20), color = "black") + 
#   guides(fill=FALSE) +
#   scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
# p.criteria_gene_total 
# ggsave(paste("../../images/singlecase/", file_suffix, "criteria_gene_total.pdf", sep = "") ,plot = p.criteria_gene_total )
# 
# # as table
# df |>
#   group_by(SYMBOL) |>
#   summarise(acmg_count_per_symbol = sum(ACMG_count)) |>
#   na.omit() |>
#   arrange(desc(acmg_count_per_symbol))
# 
# p.variants_per_criteria <- df |> 
#   ggplot(aes(x = ACMG_count, fill=..x..)) +
#   geom_histogram(binwidth = 1, color="black") +
#   xlab("No. ACMG criteria\nassigned (P)") +
#   ylab("No. variants") +
#   theme_minimal() +
#   geom_text(stat='count', aes(label=..count.., y=..count..+300), color = "black") + 
#   guides(fill=FALSE) +
#   scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
# p.variants_per_criteria
# ggsave(paste("../../images/singlecase/", file_suffix, "variants_per_criteria.pdf", sep = "") ,plot = p.variants_per_criteria , width = 9, height = 5)
# 
# # Check we only have approx. 1 "casual" variant per sample
# p.criteria_per_sample <- df %>%
#   group_by(sample) %>%
#   summarise(ACMG_count = max(ACMG_count, na.rm = TRUE))  %>%
#   ggplot(aes(x = ACMG_count, fill=..x..)) +
#   geom_histogram(binwidth = 1, color = "black") +
#   labs(x = "No. ACMG criteria\nassigned (P)", y = "No. samples") +
#   theme_minimal() +
#   geom_text(stat='count', aes(label=..count.., y=..count..+20), color = "black") + 
#   guides(fill=FALSE) +
#   scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
# p.criteria_per_sample
# ggsave(paste("../../images/singlecase/", file_suffix, "criteria_per_sample.pdf", sep = "") ,plot = p.criteria_per_sample )
# 
# # as table
# df |> 
#   group_by(sample, ACMG_count) |>
#   tally(n = "count_per_sample") |>
#   ungroup() |>
#   dplyr::select(-sample) |>
#   group_by(ACMG_count) |>
#   tally(n = "count_per_sample")
# 
# # ACMG Verdict----
# # Rules are combined using the point system described in PMID:32720330
# # Each rule triggered is assigned a number of points based on the strength of the evidence provided:
# # 
# # Supporting: 1 point
# # Moderate: 2 points
# # Strong: 4 points
# # Very Strong: 8 points
# # A total score is computed as the sum of the points from the pathogenic rules, minus the sum of the points from benign rules.
# # 
# # The total score is then compared to thresholds to assign the final verdict:
# # 	
# # Pathogenic if greater than or equal to 10,
# # Likely Pathogenic if between 6 and 9 inclusive,
# # Uncertain Significance if between 0 and 5,
# # Likely Benign if between -6 and -1,
# # Benign if less than or equal to -7.
# 
# # Benign if less than or equal to -7.
# 
# df <-  df |> dplyr::select("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5", 
#                            "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6", 
#                            "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4",
#                            everything())
# 
# 
# # Define scores for each ACMG label
# acmg_scores <- c("PVS1" = 8,
# 					  "PS1" = 4, "PS2" = 4, "PS3" = 4, "PS4" = 4, "PS5" = 4,
# 					  "PM1" = 2, "PM2" = 2, "PM3" = 2, "PM4" = 2, "PM5" = 2,
# 					  "PP3" = 1)
# 
# # Create ACMG_score column by looking up ACMG_highest in acmg_scores
# df$ACMG_score <- acmg_scores[df$ACMG_highest]
# 
# # If there are any ACMG labels that don't have a corresponding score, these will be NA. You may want to set these to 0.
# df$ACMG_score[is.na(df$ACMG_score)] <- 0
# df <- df |> dplyr::select(ACMG_score, everything())
# 
# # Total ACMG score ----
# # List of all ACMG columns
# acmg_columns <- grep("ACMG_P", colnames(df), value = TRUE)
# acmg_columns
# 
# # Mutate all ACMG columns
# df <- df %>% 
#   mutate_at(acmg_columns, function(x) acmg_scores[x])
# 
# # Replace NAs with 0 in ACMG columns only
# df[acmg_columns] <- lapply(df[acmg_columns], function(x) ifelse(is.na(x), 0, x))
# 
# # Calculate total ACMG score
# df$ACMG_total_score <- rowSums(df[acmg_columns])
# 
# df <- df |> dplyr::select(ACMG_total_score, everything())

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
ggsave(paste("../../images/singlecase/", file_suffix, "acmg_score.pdf", sep = "") ,plot = p.acmg_score )


# panel ----
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (
	(p.criteria_gene_total) / ( p.variants_per_criteria | p.criteria_per_sample ) / ( p.pathogenicity_distributions | p.acmg_score)
) + plot_annotation(tag_levels = 'A')
ggsave(paste("../../images/singlecase/", file_suffix, "patch1.pdf", sep = "") ,plot = patch1  + plot_annotation(tag_levels = 'A'), width = 8, height = 10 )

patch2 <- (
	(p.criteria_gene_total) / ( p.variants_per_criteria | p.criteria_per_sample ) / ( p.pathogenicity_distributions | p.acmg_score)
)  | (p.pathogenicity_distributions_engines_threshold) + plot_annotation(tag_levels = 'A')
# patch2
ggsave(paste("../../images/singlecase/", file_suffix, "patch2.pdf", sep = "") ,plot = patch2 + plot_annotation(tag_levels = 'A'), width = 16, height = 10 )
 
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
                ACMG_score, 
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
          desc(ACMG_score),
          sample.id)

colnames(df_report_main_text)[colnames(df_report_main_text) == 'Strong_pathogenic_GE'] <- 'Strong_patho'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'Moderate_pathogenic_GE'] <- 'Moder_patho'
colnames(df_report_main_text)[colnames(df_report_main_text) == 'Supporting_pathogenic_GE'] <- 'Suppor_patho'

saveRDS(df_report, file="../../data/", output_directory, "/df_report.Rds")

saveRDS(df_report_main_text, file="../../data/", output_directory, "/df_report_main_text.Rds")


geneset_MCL_ID <- "" #ignore pathway level info
write.csv(df_report_main_text,  paste0("../../data/", output_directory, "/AMCGuru_singlecase_df_report_main_text.csv"))

write.csv(df_report,  paste0("../../data/", output_directory, "/AMCGuru_singlecase_df_report.csv"))
