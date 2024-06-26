# Used by vsat_run.R
# this file will take VEP annotated vcf.gz
# and output tables with
# all variants, comphet, and phenotype file.

# Require: index vcf
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz

# main ----
library(VariantAnnotation)

# cat("\nNow analysing file: ", vcfFile)

# vcf import ----
# cat("\nvcf import")

tabixVcf <- Rsamtools::TabixFile(file = vcfFile)
vcf <- VariantAnnotation::readVcf(file = tabixVcf)
evcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)
rm(vcfFile)
rm(tabixVcf)
rm(vcf)

# cat("\ncsq format")
# get what I want ----
csq <- ensemblVEP::parseCSQToGRanges(x = evcf)
df_csq <- as.data.frame(csq, row.names = NULL)
rm(csq)
df_csq$rownames <-
	ensemblVEP::parseCSQToGRanges(x = evcf) |> names()

# cat("\nget genotype")
df_geno <- as.data.frame(geno(evcf)[["GT"]])
df_info <- as.data.frame(info(evcf))
rm(evcf)

df_info$rownames <- df_info |> rownames()
df_geno$rownames <- df_geno |> rownames()

# cat("\nfilter canonical")
df_csq <- df_csq |> dplyr::filter(CANONICAL == "YES")
df_csq <- df_csq |> dplyr::select(rownames, SYMBOL, HGVSp, HGVSc, Consequence, IMPACT, gnomAD_AF)
df_info <- df_info |> dplyr::select(rownames, AC, AF, AN)

# cat("\nadd info")
df_merge <- merge(df_geno, df_info, by = "rownames")
rm(df_geno, df_info)

# get CRITICAL_n_samples ----
n_sample_col_start <- 1+1
n_sample_col_end <- 1+CRITICAL_n_samples
df_merge <- df_merge |> dplyr::select(1:n_sample_col_end,  AC, AF, AN)

# cat("\nadd csq")
df_main <- merge(df_merge, df_csq, by = "rownames")
rm(df_merge)
rm(df_csq)

