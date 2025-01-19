#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 12G
#SBATCH --time 05:00:00
#SBATCH --job-name=filtervep_small
#SBATCH --output=./log/vep/filter_vep_impact_small_%J.out
#SBATCH --error=./log/vep/filter_vep_impact_small_%J.err
#SBATCH --comment="scitas_cost"

# CONDA ENV
# activate before submission
# conda activate vep

# 1 hour for impact filter 1000 samples with small vep
# 2 hours for 1000 smaples with big vep dbnsfp (~700 columns)

echo "START AT $(date)"
set -e

# Master variables
source ./variables.txt

# Path
INPUT_DIR="${WORK_DIR}/pre_annotation_output"
OUTPUT_DIR="${WORK_DIR}/annotation"

module load intel # jed

echo "filter_vep for either moderate or high impact variants."

# Keep either moderate or high impact variants
filter_vep \
-i ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small.vcf \
-o ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact.vcf \
--force_overwrite \
--filter "IMPACT in MODERATE,HIGH"

echo "wc -l before filter:"
wc -l ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small.vcf
echo "wc -l after filter:"
wc -l ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact.vcf

echo "END AT $(date)"
