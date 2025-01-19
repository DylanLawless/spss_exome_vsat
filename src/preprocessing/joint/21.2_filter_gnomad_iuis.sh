#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 12G
#SBATCH --time 05:00:00
#SBATCH --job-name=filtergnomad_big
#SBATCH --output=./log/vep/filter_vep_gnomad_big_%J.out
#SBATCH --error=./log/vep/filter_vep_gnomad_big_%J.err
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

# echo "filter_vep for gnomad_af freq."
# filter_vep \
# -i ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis.vcf \
# -o ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad.vcf \
# --force_overwrite \
# --filter "gnomAD_AF < 0.01"

echo "filter_vep for cohort_af freq."
filter_vep \
-i ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad.vcf \
-o ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1.vcf \
--force_overwrite \
--filter "AF < 0.01 or not AF"

echo "wc -l before filter:"
wc -l ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis.vcf
echo "wc -l after gnomad filter:"
wc -l ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad.vcf
echo "wc -l after AF0.1 filter:"
wc -l ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1.vcf \

echo "END AT $(date)"
