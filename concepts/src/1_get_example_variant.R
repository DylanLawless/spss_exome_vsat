#!/bin/bash
# setwd("../../src/ACMGuru_post_ppi/")
# source("ACMGuru_post_ppi_vcurrent.R")
# setwd("../../concepts/src/")

geneset_MCL_ID <- c(22, 586)
df_report_main_text <- readRDS(file=paste0("../../data/ACMGuru_post_ppi/df_report_main_text_", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

names(df_report_main_text)

# select features
df_report_set <- df_report_main_text |> 
  dplyr::select(sample.id, rownames, 
                CHROM, REF, ALT, 
                POS, start, end, width, 
                Gene, SYMBOL, HGNC_ID, 
                HGVSp, HGVSc, Consequence,  
                IMPACT, genotype,
                Feature_type, Feature, BIOTYPE, VARIANT_CLASS, CANONICAL,
                    ) |> 
  filter(IMPACT == "HIGH") |>
  head(1)

df_report_set |> names()

# set example sample.id
df_report_set$sample.id <- "Canton_001"

# save example variant
output_directory <- "../../data/"
saveRDS(df_report_set, file=paste0(output_directory, "example_variant.Rds"))
