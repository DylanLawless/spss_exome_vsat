#!/bin/bash
#SBATCH --array=0-391
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 2G
#SBATCH --time 00:10:00
#SBATCH --job-name=tbi
#SBATCH --output=./log/hc/tbi_%J.out
#SBATCH --error=./log/hc/tbi_%J.err
#SBATCH --comment="scitas_cost"

# 0.08 min per file. 10 min per task should be OK. 

set -e
echo "START AT $(date)"

# Master variables
source variables.txt

GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# /////////////////////////////////////////////////////////////////////////////
# CHECK IF VCF RENAME WAS REQUIRED
# CHECK IF VCF RENAME WAS REQUIRED
# CHECK IF VCF RENAME WAS REQUIRED
# /////////////////////////////////////////////////////////////////////////////

OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output"
# OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output_rename"

# bgzip files instead of gzip
# then remake the tbi with IndexFeatureFile \

cd $OUTPUT_DIR

# Declare
declare -a total_file
for i in *.g.vcf.gz; do total_file+=($i); done
filebase=${total_file[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${filebase} .g.vcf.gz)

$GATK IndexFeatureFile \
	-I ${sample_id}.g.vcf.gz \
	-O ${sample_id}.g.vcf.gz.tbi

echo "END AT $(date)"
