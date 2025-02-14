#!/bin/bash

# This version was run locally on macOS for some tests
# Confirmed that the included result is from this version with covariates and filters included
# Updating 3_assoc.sh to match

echo "START AT $(date)"
set -e
# module load intel # jed

cd ../../data/variant_level_extended

# plink \
~/tools/plink_mac_20241022/plink \
	--bfile output1 \
	--geno 0.05 \
	--hwe 1e-6 \
	--make-bed \
	--allow-no-sex \
	--out output3

# logistic model with confounders
# plink --bfile output3 --logistic --covar file.pca

# assoc model with control/case (1/2), PCA, covar
# plink 
~/tools/plink_mac_20241022/plink \
--bfile output3 --assoc --covar bcftools_gatk_norm_pca.eigenvec_covar --allow-no-sex

