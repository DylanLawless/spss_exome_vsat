#!/bin/bash
echo "START AT $(date)"
set -e
# module load intel # jed

# Master variables
source ./variables.txt

cd  ${VARLEVEL_DIR}

# get PCA with covariates
# PC1-10, sex, age days^2, study site, and ICU stay
# "../../data/joint/covar/joint.csv"
# cp ${WORK_DIR}/pca_output/bcftools_gatk_norm_pca.eigenvec ./
cp ${WORK_DIR}/pca_output/bcftools_gatk_norm_pca.eigenvec_covar ./


