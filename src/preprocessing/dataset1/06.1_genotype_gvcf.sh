#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 25
#SBATCH --cpus-per-task 2
#SBATCH --mem 20G
#SBATCH --time 01:00:00
#SBATCH --job-name=genotype_gvcf
#SBATCH --output=./log/genotype_gvcf/genotype_gvcf_%J.out
#SBATCH --error=./log/genotype_gvcf/genotype_gvcf_%J.err
#SBATCH --comment="scitas_cost"

# array for all intervals: 1-22 X Y M = 25
# sbatch --array=0-24 06.1_genotype_gvcf.sh
# Time: for 380 samples, chr 1: Processed 928140 total variants in 14.2 minutes.

echo "START AT $(date)"
set -e

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset1/data"

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
INPUT_DIR="${SCRATCH_DIR}/genomics_db_import_output"
OUTPUT_DIR="${SCRATCH_DIR}/genotype_gvcf_output"

# Path: ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Prepare input
declare -a NUMBER
for j in {1..22} X Y M; do NUMBER+=($j); done
INDEX=${NUMBER[$SLURM_ARRAY_TASK_ID]}
cd $INPUT

# Run GenotypeGVCFs
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms18G -Xmx18G -XX:ParallelGCThreads=2" GenotypeGVCFs \
	-R $REF \
	-V gendb://${INPUT_DIR}/chr${INDEX}_gdb \
	-O ${OUTPUT_DIR}/chr${INDEX}.g.vcf

echo "FINISH AT $(date)"

