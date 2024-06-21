# This file will take VEP annotated vcf.gz
# and output tables with
# all variants, comphet, and phenotype file.

# Require: index vcf
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz

# main ----
# import ----
# if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")
# BiocManager::install("ensemblVEP")
library(VariantAnnotation)

cat("\nNow analysing file: ", vcfFile)

# vcf import ----
cat("\nvcf import")

tabixVcf <- Rsamtools::TabixFile(file = vcfFile)
vcf <- VariantAnnotation::readVcf(file = tabixVcf)
evcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)

cat("\ncsq format")

# Extract mandatory fields ----
# mandatory_fields <- data.frame(CHROM = vcf@rowRanges@seqnames,
#                                POS = start(vcf@rowRanges),
#                                ID = vcf@rowRanges@elementMetadata$ID,
#                                REF = vcf@rowRanges@elementMetadata$REF,
#                                ALT = vcf@rowRanges@elementMetadata$ALT)

# Extract mandatory fields into a dataframe and include rownames as a column
mandatory_fields <- data.frame(
  ROW_ID = names(rowRanges(vcf)),  # Use rownames from rowRanges
  CHROM = seqnames(rowRanges(vcf)),
  POS = start(rowRanges(vcf)),
  REF = as.character(ref(vcf)),  # Convert DNAStringSet to character
  ALT = sapply(alt(vcf), function(alts) {
    paste(sapply(alts, as.character), collapse = ",")
  }),
  stringsAsFactors = FALSE
)

# rm(vcfFile)
# rm(tabixVcf)
# rm(vcf)

# get vcf# get what I want ----
csq <- ensemblVEP::parseCSQToGRanges(x = evcf)
df_csq <- as.data.frame(csq, row.names = NULL)
# rm(csq)
df_csq$rownames <-
	ensemblVEP::parseCSQToGRanges(x = evcf) |> names()

cat("\nget genotype")
df_geno <- as.data.frame(geno(evcf)[["GT"]])
#evcf[csq$"VCFRowID"]
df_info <- as.data.frame(info(evcf))
rm(evcf)

df_info$rownames <- df_info |> rownames()
df_geno$rownames <- df_geno |> rownames()

cat("\nfilter canonical")
df_csq <-
	df_csq |> dplyr::filter(CANONICAL == "YES") # remove filter later, reduces size for testing.

cat("\nadd info")
df_merge <- merge(df_geno, df_info, by = "rownames")
rm(df_geno, df_info)
#names(df_merge)
df_merge <- df_merge |> dplyr::select(-CSQ) # now drop the original  list column.

cat("\nadd csq")
df_main <- merge(df_merge, df_csq, by = "rownames")
# rm(df_merge)
# rm(df_csq)

# Optional: Merge mandatory fields with the main dataframe
df_main_meta <- merge(df_main, mandatory_fields, by.x = "rownames", by.y = "ROW_ID")

df_main <- df_main_meta
