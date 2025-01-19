#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 12G
#SBATCH --time 05:00:00
#SBATCH --job-name=filtergnomad_small
#SBATCH --output=./log/vep/filter_vep_gnomad_small_%J.out
#SBATCH --error=./log/vep/filter_vep_gnomad_small_%J.err
#SBATCH --comment="scitas_cost"

# CONDA ENV
# activate before submission
# conda activate vep

echo "START AT $(date)"
set -e

# Master variables
source ./variables.txt

# Path
INPUT_DIR="${WORK_DIR}/pre_annotation_output"
OUTPUT_DIR="${WORK_DIR}/annotation"

module load intel # jed

echo "filter_vep for gnomad_af freq."
filter_vep \
-i ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact.vcf \
-o ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad.vcf \
--force_overwrite \
--filter "gnomAD_AF < 0.2"

echo "wc -l before filter:"
wc -l ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact.vcf
echo "wc -l after filter:"
wc -l ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad.vcf

echo "END AT $(date)"
