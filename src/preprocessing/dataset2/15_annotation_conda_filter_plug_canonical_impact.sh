#!/bin/bash
#SBATCH --array=0-24
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 12G
#SBATCH --time 02:00:00
#SBATCH --job-name=filterimpact
#SBATCH --output=./log/vep/filter_impact_%J.out
#SBATCH --error=./log/vep/filter_impact_%J.err
#SBATCH --comment="scitas_cost"

# CONDA ENV
# activate before submission
# conda activate vep

# 1 hour for impact filter 1000 samples

echo "START AT $(date)"
set -e
module load intel # jed

# Master variables
source ./variables.txt

# Path
OUTPUT_DIR="${WORK_DIR}/annotation"
FILENAME="bcftools_gatk_norm_vep_conda_plug_maf.recode"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Prepare input
declare -a NUMBER
for j in {1..22} X Y M; do NUMBER+=($j); done
INDEX=${NUMBER[$SLURM_ARRAY_TASK_ID]}

filter_vep \
-i ${OUTPUT_DIR}/${FILENAME}_canonical_chr_${INDEX}.vcf \
-o ${OUTPUT_DIR}/${FILENAME}_canonical_impact_chr_${INDEX}.vcf \
--force_overwrite \
--filter "IMPACT in HIGH,MODERATE"
# --filter "IMPACT = HIGH OR IMPACT = MODERATE" # This seems to only include MODERATE if also HIGH

# Check count
printf "[1] Variant transcripts with MODERATE but not HIGH: "
grep "MODERATE" \
${OUTPUT_DIR}/${FILENAME}_canonical_impact_chr_${INDEX}.vcf |\
grep -v "HIGH" | wc -l

printf "[2] Variant transcripts with HIGH but not MODERATE: "
grep "HIGH" \
${OUTPUT_DIR}/${FILENAME}_canonical_impact_chr_${INDEX}.vcf |\
grep -v "MODERATE" | wc -l

printf "[3] Variant transcripts with HIGH and MODERATE:     "
grep "HIGH" \
${OUTPUT_DIR}/${FILENAME}_canonical_impact_chr_${INDEX}.vcf |\
grep "MODERATE" | wc -l

echo "END AT $(date)"
