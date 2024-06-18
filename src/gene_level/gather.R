# Used by vsat_run.R
# gathers genotypes from wide to long

# cat("\nVariables set :")
# cat(paste0("\nGnomad Freq< : ", gnomad_rare_thresh,
# 			  "\nGnomad count < : ", gnomad_count_thresh,
# 			  "\nCohort Freq< : ", cohort_freq_max,
# 			  "\nCohort count < : ", cohort_max_carriers))

vcfFile <- file_list[f]
source("vcf_to_tables.R")

df_main$AC<- as.numeric(df_main$AC)

# gather columns by CRITICAL_n_samples ----
# uses column number defined by
# n_sample_col_start, n_sample_col_end (see vcf_to_tables)
# n_sample_col_start <- 1+1
# n_sample_col_end <- 1+CRITICAL_n_samples

# cat("\nGather wide to long")
# cat("\nGathering: ", CRITICAL_n_samples, " sample columns." )
# cat("\nCRITICAL CHECK - Gathering: Columns ", n_sample_col_start, " to ", n_sample_col_end)
df <-
	tidyr::gather(df_main, n_sample_col_start:n_sample_col_end, key = "sample", value = "genotype_call") |>
	dplyr::select(sample, genotype_call, SYMBOL, HGVSp, HGVSc, Consequence,IMPACT, everything()) 
df$SNP <- df$"rownames"
rm(n_sample_col_start, n_sample_col_end)

# cat("\nFinished gathering sample: ", f)
