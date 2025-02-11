library(dplyr)
library(tidyr)

cat("\nVariables set :")
cat(paste0("\nGnomad Freq< : ", gnomad_rare_thresh,
			  "\nGnomad count < : ", gnomad_count_thresh,
			  "\nCohort Freq< : ", cohort_freq_max,
			  "\nCohort count < : ", cohort_max_carriers))

vcfFile <- file_list[f]
source("../stand_alone_vcf_to_table/vcf_to_tables.R")

df_main$AC<- as.numeric(df_main$AC)
df_main$AF.x<- as.numeric(df_main$AF.x)
df_main$MLEAC <- as.numeric(df_main$MLEAC)
df_main$MLEAF <- as.numeric(df_main$MLEAF)
df_main$PG<- as.character(df_main$PG)
df_main$RAW_MQandDP<- as.character(df_main$RAW_MQandDP)
df_main$OLD_VARIANT<- as.character(df_main$OLD_VARIANT)
df_main$ReadPosRankSum <- as.character(df_main$ReadPosRankSum)

# get CRITICAL_n_samples ----
n_sample_col_start <- 1+1
n_sample_col_end <- 1+CRITICAL_n_samples

cat("\nGather wide to long")
cat("\nGathering: ", CRITICAL_n_samples, " sample columns." )
cat("\nCRITICAL CHECK - Gathering: Columns ", n_sample_col_start, " to ", n_sample_col_end)

df <-
	tidyr::gather(df_main, n_sample_col_start:n_sample_col_end, key = "sample", value = "genotype_call") |>
	dplyr::select(sample, genotype_call, SYMBOL, HGVSp, HGVSc, Consequence,IMPACT, everything()) 
df$SNP <- df$"rownames"
rm(n_sample_col_start, n_sample_col_end)

cat("\nFinished gathering sample: ", f)
