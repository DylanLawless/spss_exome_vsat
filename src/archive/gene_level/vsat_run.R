# Set paramenters in submission script

# read in current_MCLID from command line arguments. 
# this is required because I submit one job per MCLID (variant set ID / protein pathway ID). 
# the array allows this R code to run in parallel on each of the variant sets, 
# but this ID must be read from the scheduler 
# (SLRUM scheduler uisng SBATCH on EPLF's HPC: SCITAS).
# if we use "trailingOnly" we should be able to use only the input args 1,2,3 etc. 
# otherwise they would be args 6,7,8 etc.

# access the arguments
args <- commandArgs(trailingOnly = TRUE)

# print the args
cat("all args passed to R:")
print(args)

# set the variables read from the subssion script
current_MCLID <- args[1]
version_suffix <- args[2]
file_suffix <- args[3]
filter_method <- args[4]
gnomad_rare_thresh <- as.numeric(args[5])
cohort_max_carriers <- as.integer(args[6])
CRITICAL_n_samples <- as.integer(args[7])
chr_num <- as.integer(args[9])


cat("\nMCLID:", current_MCLID)
cat("\nversion:", version_suffix)
cat("\nfile_suffix:", file_suffix)
cat("\nfilter_method:", filter_method)
cat("\ngnomad_thresh:", gnomad_rare_thresh)
cat("\ncohort_thresh:", cohort_max_carriers)
cat("\nsamples:", CRITICAL_n_samples, "\n")
cat("\nchr num:", chr_num, "\n")

# redirect R output as a log file to save P-value from skat and other info
# this log will be converted to the final output from all VSAT.
sink(paste("log/log_", version_suffix, "/log_skat_", file_suffix, "_chr_", chr_num, "_MCL_", current_MCLID, ".log", sep = ""), append = TRUE)

# /////////////////////////////////////////////////////////////////////////////
# Required packages ----
# /////////////////////////////////////////////////////////////////////////////

# the "::" operator is used for namespace collisions, (dplyr).
# it is also used for several required packages in vcf_to_tables to avoiding loading the full packages.
library(dplyr)
library(tidyr)
library(SKAT)
# called in vcf_to_tables:
# ensemblVEP
# VariantAnnotation
# Rsamtools

# /////////////////////////////////////////////////////////////////////////////
# Project specifics ----
# /////////////////////////////////////////////////////////////////////////////

# VCFtool filter will show how many samples to use in CRITICAL_n_samples during submission

# frequency variables -----
gnomad_total <- 140000 
# gnomad_rare_thresh  set in submission script
gnomad_count_thresh <- gnomad_total*gnomad_rare_thresh
# cohort_max_carriers set in submission script
cohort_freq_max <- cohort_max_carriers/CRITICAL_n_samples

# confirm that objects do not already exist
suppressWarnings(rm(results_df,file_list))

# get PCs
X <- read.csv("../../data/joint/pca_output/bcftools_gatk_norm_pca.eigenvec",
				  sep = " ",
				  header = F)

# prep qualifying candidate variants
# generate file_list: read chromosome number from array task id
file_list <- paste0("../../data/joint/annotation/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad_chr_", chr_num, ".vcf.gz")

# file_list <- c("../../data/joint/annotation/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad_chr_21.vcf.gz")
# file_list <- c(
# 	paste0("../../data/joint/annotation/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad_chr_", 1:22, ".vcf.gz")
# 	# ,"../data/annotation/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad_chr_X.vcf.gz"
# 	# ,"../data/annotation/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad_chr_Y.vcf.gz"
# )


# keep gene in current chr if in pathway
# merge all to one df_pathway
# for (f in 6) { # do this for testing a specific chr. e.g. chr6
df_pathway_list <- list()
for (f in 1:length(file_list)) {
	source("gather.R") # note that this includes source("vcf_to_tables_2.R")
	source("genotype_clean.R")
	source("progress_bar.R")
	

# initialize an empty dataframe to store the results
df_all <- data.frame()
total_output <- data.frame()
df_pathway <- data.frame()

# clusters_HIGH_joint <- read.csv("../../data/ppi/mcl_clusters_recode_df_whole_genome_v5_c7_700.csv")
# clusters_HIGH_joint$SYMBOL <- clusters_HIGH_joint$Items

# # GENE_LEVEL
# # to run at gene-level instead of MCL_ID level we simply call each gene a 1-gene pathway and use the SYMBOL as an MCL_ID to keep the code consistent. Now MCL_ID == SYMBOL (gene). Note that we needd a gene list for all genes in the study
# # clusters_HIGH_joint$MCL_ID <- clusters_HIGH_joint$ID # pathway level
# clusters_HIGH_joint$MCL_ID <- clusters_HIGH_joint$Items

# clusters_HIGH_joint <- clusters_HIGH_joint |> dplyr::select(MCL_ID, SYMBOL)
# pathway_list <- clusters_HIGH_joint$MCL_ID |> unique()

# # /////////////////////////////////////////////////////////////////////////////
# # Gather gene per chr ----
# # /////////////////////////////////////////////////////////////////////////////

# # 1-947 pathways
# rm(df_pathway)	
# df_pathway <- data.frame()

# current_MCL_set <- clusters_HIGH_joint |> filter(MCL_ID == current_MCLID ) 
# current_SYMBOL <- current_MCL_set$SYMBOL |> unique() 
# current_SYMBOL_length <- current_SYMBOL |> length()
# cat("Now analysisng pathway ID:", current_MCLID, "\n")
# cat("Which consists of genes:", current_SYMBOL, "\n")	
# cat("Total number of genes:", current_SYMBOL_length, "\n")


# get the 90% genotype frequency based on AN
df_AN_GENOpc <- df %>%
	dplyr::select(SNP, AN) %>%
	summarise(max_AN = max(AN, na.rm = TRUE),
				 min_AN = min(AN, na.rm = TRUE),
				 AN_GENOpc = max_AN * 0.9) 

df_AN_GENOpc <- df_AN_GENOpc$AN_GENOpc

# Plot for testing, shows the AN for batch effect which is visible at around 44% and subsequently filtered out. 
# df |> dplyr::select(SNP, AN) |> group_by(SNP) |> unique() |> ggplot(aes(x=SNP, y=AN)) + geom_point() + geom_hline(yintercept = df_AN_GENOpc)

df <- df |> filter(AN > df_AN_GENOpc)
# rm(df_AN, df_AN_GENOpc)
rm(df_AN_GENOpc)

# Now return NA to 0 genotype for SKAT matrix
df <- df |> mutate(genotype = ifelse(is.na(genotype), 0, genotype))

	rm(df_main)

	# uses variables from submission script
	if(filter_method == "HIGH"){
	  df <- df |> filter(IMPACT == "HIGH")
	} else if (filter_method == "MODERATE-HIGH"){
	  df <- df |> filter(IMPACT %in% c("HIGH", "MODERATE"))
	}
	df <- df |> filter(gnomAD_AF < gnomad_rare_thresh)
	
	# hold the current info
	# df_hold <- merge(current_MCL_set, df )
	# GENE_LEVEL
	df_hold <- df



# current_MCL_set <- clusters_HIGH_joint |> filter(MCL_ID == current_MCLID ) 
# current_SYMBOL <- current_MCL_set$SYMBOL |> unique() 
current_MCL_set <- current_MCLID # GENE_LEVEL
current_SYMBOL <- current_MCLID |> unique() # GENE_LEVEL
current_SYMBOL_length <- current_SYMBOL |> length()
cat("Now analysisng pathway ID:", current_MCLID, "\n")
cat("Which consists of genes:", current_SYMBOL, "\n")	
cat("Total number of genes:", current_SYMBOL_length, "\n")
	
	# save to list for memory rather than a large df
	df_pathway_list[[f]] <- df_hold
	cat("Now analysisng chr", f, "\n")
	cat("which has genes:", unique(df_hold$SYMBOL), "\n")
	
}

df_pathway <- do.call(rbind, df_pathway_list)

# cleanup for memory
gc()

# GENE_LEVEL
# Use MCL_ID equal to gene SYMBOL
df_pathway$MCL_ID <- df_pathway$SYMBOL

cat("Fished pathway", unique(df_pathway$MCL_ID), "\n")
cat("Found genes:", unique(df_pathway$SYMBOL), "\n")
cat("\n")
cat("# \\\\\\\\\\\\\\\\\\\\", "\n")
cat("Summary:", "\n")
cat("# \\\\\\\\\\\\\\\\\\\\", "\n")
cat("Finshed analysisng MCL pathway ID:", current_MCLID, "\n")
# cat("Our data had:", length(unique(df_pathway$SYMBOL)), "\n")
# cat("genes:", unique(df_pathway$SYMBOL), "\n")
cat("total MCL genes:", current_SYMBOL_length, "\n")
cat("Number chromosomes:", length(file_list), "\n")
# cat("Our pathway data had rows:", length(df_pathway), "\n")
cat("\n")

# df_pathway[1,1]
# head(df_pathway)

# /////////////////////////////////////////////////////////////////////////////
# filter ----
# /////////////////////////////////////////////////////////////////////////////
# For gene and pathways use: SYMBOL or MCL_ID.

df <- df_pathway
df$cohort_pheno <- df$sample
df$cohort_pheno[grep("^setpt", df$sample)] <- "0"
df$cohort_pheno[!grepl("^setpt", df$sample)] <- "1"

OK <- df |> dplyr::select(sample, genotype, cohort_pheno, MCL_ID, SYMBOL, "rownames") 
OK$SNP <- df$"rownames"
OK <- OK |> na.omit(MCL_ID)

OK_filtered <- OK |>
	filter(cohort_pheno == 1) |>
	group_by(MCL_ID) |>
	summarize(total_genotype = sum(as.numeric(genotype))) 

cat("Our pathway data after filt for cases had rows:", length(df_pathway), "\n")
# OK_filtered[1,1]
head(OK_filtered)

# modify if you want variants which are missing in cases but present in control, such as a protective variant model.
symbols_to_keep <- OK_filtered |>
	filter(total_genotype > 0) |>
	dplyr::select(MCL_ID)

OK_filtered <- OK |>
	filter(MCL_ID %in% symbols_to_keep$MCL_ID)

# temporary - this one gene caused error in test data but is probably fine
OK_filtered <- OK_filtered |> filter(! SYMBOL == "RMRP") 

OK <- OK_filtered

# /////////////////////////////////////////////////////////////////////////////
# Skip if obj ----
# /////////////////////////////////////////////////////////////////////////////

if (!exists("obj")) {
	
	phenotype <- OK |>
		dplyr::select(sample, cohort_pheno) |>
		unique() |>
		dplyr::select(cohort_pheno)
	
	phenotype <- phenotype$cohort_pheno |> as.numeric()
	
	OK$genotype |> unique()
	OK$cohort_pheno |> unique()

	# check the number of unique samples in the filtered dataset
	n_pheno <- sort(phenotype, decreasing = TRUE) |> unique() |> length()
	n_pc <- X[,3] |> length()
	cat("Before obj, number of pheno:", n_pheno, "\n")
	cat("Before obj, number of PC1:", n_pc, "\n")

	# compare the number of unique samples to CRITICAL_n_samples
	if (n_pheno != 2) {
		warning("Number of phenotypes in the filtered dataset does not eqaul 2 (case/control).")
	}
	
	if (n_pc != CRITICAL_n_samples) {
		warning("Number of PC1 in the filtered dataset does not match CRITICAL_n_samples.")
	}
	
	# get unique SYMBOLs
	unique_symbols <- unique(OK$MCL_ID)
	
	# compute SKAT test statistic using bootstrap resampling and covariates
	obj <- SKAT_Null_Model(phenotype ~ X[,3] + X[,4] + X[,5] + X[,6],
							out_type="D", 
							type.Resampling="bootstrap", 
							n.Resampling = 100000,
							Adjustment=TRUE)
}

# initialize empty dataframe to store results
results_df <- data.frame(MCL_ID = character(), p_value = numeric(), stringsAsFactors = FALSE)

sink()

# /////////////////////////////////////////////////////////////////////////////
# SKAT RUN ----
# /////////////////////////////////////////////////////////////////////////////

# load the package
library(parallel)

skat_run_parallel <- function(i) {
	# set up log file
	sink(paste("log/log_", version_suffix, "/", "log_skat_", file_suffix, "_chr_", chr_num, "_SYMBOL_", i, ".log", sep = ""), append = TRUE)
	
	# filter for the current ~SYMBOL~ MCL_ID
	OK_filtered <- OK |> filter(MCL_ID == i)
	
	# check the number of unique samples in the filtered dataset
	n_unique_samples <- sort(OK_filtered$sample, decreasing = TRUE) |> unique() |> length()
	
	cat("Before skat, number of sample:", n_unique_samples, "\n")
	cat("Before skat, number of critial sample:", CRITICAL_n_samples, "\n")
	
	# compare the number of unique samples to CRITICAL_n_samples
	if (n_unique_samples != CRITICAL_n_samples) {
	  warning("Number of unique samples in the filtered dataset does not match CRITICAL_n_samples.")
	}
	
	# select relevant columns
	OK_filtered <- OK_filtered |> dplyr::select(sample, genotype, cohort_pheno, SNP) 
	
	# create phenotype vector
	phenotype <- OK_filtered |>
		dplyr::select(sample, cohort_pheno) |>
		unique() |>
		dplyr::select(cohort_pheno) 
	
	phenotype <- phenotype$cohort_pheno |> as.numeric()
	
	# prevent erorr with skat do to list when duplicates
	OK_filtered <- distinct(OK_filtered, sample, genotype, cohort_pheno, SNP) 

	# add the allele count for each genotype and phenotype
	OK_summary <- OK_filtered |> 
		group_by(cohort_pheno, genotype) |> 
		summarise(group_frequency = n(), .groups = "drop") |> 
		mutate(cohort_pheno = case_when(
			cohort_pheno == "0" ~ "control",
			cohort_pheno == "1" ~ "case",
			TRUE ~ cohort_pheno)) |>
		mutate(phenotype_genotype = paste(cohort_pheno, genotype, sep = "_")) |>
		ungroup() |>
		dplyr::select(phenotype_genotype, group_frequency) 

	OK_summary <- spread(OK_summary, phenotype_genotype, group_frequency)

	OK_var_count <- OK_filtered |> 
		dplyr::select(SNP) |> 
		unique() |>
		summarise(SNP = n())

	# this can be improved instead of using "OK" again, but testing to avoid error
	OK_gene_count <- OK |> 
		# filter(MCL_ID == unique_symbols[i]) |> 
		filter(MCL_ID == i) |> 
		dplyr::select(SYMBOL) |> 
		unique() |>
		summarise(SYMBOL = n())
	
	# pivot the data to create the Z matrix
	Z <- pivot_wider(OK_filtered, id_cols = sample, names_from = SNP, values_from = genotype)
	
	# remove the "sample" column from Z
	Z <- Z[, -1]
	
	# convert matrix elements to numeric
	Z <- apply(Z, 2, as.numeric)
	
	# compute SKAT test statistic using bootstrap resampling and covariates
	result_long_resamp_cov <- SKAT(Z, obj, method="SKATO")
	result <- Get_Resampling_Pvalue(result_long_resamp_cov)$p.value
	# result <- result_long_resamp_cov$p.value
	
	results_df <- rbind(results_df,
	data.frame(
		# MCL_ID = unique_symbols[i], # only at for pathway level
		MCL_ID = i, # NEEDS FIX AT for gene level
		p_value = result,
		OK_gene_count,
		OK_var_count,
		OK_summary,
		stringsAsFactors = FALSE))

	cat("Results for MCL, p_val, gene_count, var_count, ACs for case control:", paste0(results_df, collapse = " : "), "\n")
	# example output: this is correct as 6 unique variants gives (3280+8+2346+6)/6 = 940 samples.
	# MCL, p_val, var_count, ACs for case control: 600 : 0.576987424157896 : 6 : 3280 : 8 : 2346 : 6

	# workaround while awaiting pull request to SKAT package #17 3078cee, 
	# the SNP-ID from Z are only included in $test.snp.mac when there is more than 1 SNP in Z.
	# check if there is only one variant in the test set
	if (dim(Z)[2] == 1) {
		# if so, manually set the SNP ID in the output
		result_long_resamp_cov$test.snp.mac <- as.data.frame(result_long_resamp_cov$test.snp.mac)
		rownames(result_long_resamp_cov$test.snp.mac) <- colnames(Z)
	} else {
		# if not, proceed as before
		result_long_resamp_cov$test.snp.mac <- as.data.frame(result_long_resamp_cov$test.snp.mac)
		rownames(result_long_resamp_cov$test.snp.mac) <- rownames(result_long_resamp_cov$test.snp.mac)
	}

	# save the SNPs and count per MCL test
	test.snp.mac <- as.data.frame(result_long_resamp_cov$test.snp.mac) 
	test.snp.mac$SNP <- rownames(test.snp.mac)
	test.snp.mac <- test.snp.mac %>% setNames(c("test.snp.mac", "SNP"))
	rownames(test.snp.mac) <- NULL
	
	print(data.frame(MCL_ID = i, test.snp.mac, stringsAsFactors = FALSE))
	rm(test.snp.mac)
	
	# Close the log file
	sink()
	
	# return results_df at the end
	return(results_df)
	gc()
}

# determine the number of cores to use
# num_cores <- detectCores() - 1 # cores on node - note memory limit if use this
num_cores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1)) - 1


# Use mclapply to parallelize
results_list <- mclapply(unique_symbols, skat_run_parallel, mc.cores = num_cores)

# Combine results into one data frame
results_df <- do.call(rbind, results_list)
gc()

