#!/bin/bash
echo "START AT $(date)"
set -e
# module load intel # jed

# Master variables
source ./variables.txt

cd  ${VARLEVEL_DIR}

# get PCA
cp ${WORK_DIR}/pca_output/bcftools_gatk_norm_pca.eigenvec ./


