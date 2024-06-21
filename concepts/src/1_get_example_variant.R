#!/bin/bash
setwd("../../src/ACMGuru_post_ppi/")
source("ACMGuru_post_ppi_vcurrent.R")
setwd("../../concepts/src/")

# select features
df_report_set <- df |> 
  dplyr::select(sample, rownames, 
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

# save example variant
output_directory <- "../../data/"
saveRDS(df_report_set, file=paste0(output_directory, "example_variant.Rds"))
