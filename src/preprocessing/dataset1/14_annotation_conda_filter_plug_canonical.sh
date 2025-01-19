#!/bin/bash
#SBATCH --array=0-24
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 12G
#SBATCH --time 02:00:00
#SBATCH --job-name=filtercanonical
#SBATCH --output=./log/vep/filter_canonical_%J.out
#SBATCH --error=./log/vep/filter_canonical_%J.err
#SBATCH --comment="scitas_cost"

# CONDA ENV
# activate before submission
# conda activate vep

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
-i ${OUTPUT_DIR}/${FILENAME}_chr_${INDEX}.vcf.gz \
-o ${OUTPUT_DIR}/${FILENAME}_canonical_chr_${INDEX}.vcf \
--force_overwrite \
--filter "CANONICAL is YES" 

echo "END AT $(date)"
