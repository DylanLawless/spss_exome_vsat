#!/bin/bash

echo "START AT $(date)"
set -e
# module load intel # jed

# Master variables
source ./variables.txt
cd  ${VARLEVEL_DIR}

plink \
	--bfile output1 \
	--geno 0.05 \
	--hwe 1e-6 \
	--make-bed \
	--allow-no-sex \
	--out output3 

# logistic model with confounders
# plink --bfile output3 --logistic --covar file.pca

# assoc model with control/case (1/2), PCA, covar
plink --bfile output3 --assoc --covar bcftools_gatk_norm_pca.eigenvec_covar --allow-no-sex

