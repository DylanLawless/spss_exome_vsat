#!/bin/bash
echo "START AT $(date)"
set -e
# module load intel # jed

# Master variables
source ./variables.txt

cd  ${VARLEVEL_DIR}

# get PCA with covariates
# cp ${WORK_DIR}/pca_output/bcftools_gatk_norm_pca.eigenvec ./
cp ${WORK_DIR}/pca_output/bcftools_gatk_norm_pca.eigenvec_covar ./


