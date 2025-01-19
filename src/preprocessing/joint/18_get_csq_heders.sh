#!/bin/bash
# This script gest the CSQ INFO headers added by VEP

echo "START AT $(date)"
set -e
module load intel # jed

# Master variables
source ./variables.txt

# Path
INPUT_DIR="${WORK_DIR}/pre_annotation_output"
OUTPUT_DIR="${WORK_DIR}/annotation"
BCFTOOLS="/home/lawless/bcftools/bcftools"

${BCFTOOLS} +split-vep -l \
	${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda.vcf |\
	cut -f2 \
	> ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_head_vep

${BCFTOOLS} +split-vep -l \
	${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small.vcf |\
	cut -f2 \
	> ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_head_vep
echo "END AT $(date)"
