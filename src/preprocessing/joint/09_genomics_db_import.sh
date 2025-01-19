#!/bin/bash 
#SBATCH --array=0-24
#SBATCH --nodes 1
#SBATCH --ntasks 25
#SBATCH --cpus-per-task 2
#SBATCH --mem 20G
#SBATCH --time 03:00:00
#SBATCH --job-name=genomics_db_import
#SBATCH --output=./log/genomics_db_import/genomics_db_import_%J.out
#SBATCH --error=./log/genomics_db_import/genomics_db_import_%J.err
#SBATCH --comment="scitas_cost"

# ~30 min per task
# array for all intervals: 1-22 X Y M = 25

echo "START AT $(date)"
set -e

# Master variables
source variables.txt

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Path: ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Path
# /////////////////////////////////////////////////////////////////////////////
# CHECK IF VCF RENAME WAS REQUIRED
# CHECK IF VCF RENAME WAS REQUIRED
# CHECK IF VCF RENAME WAS REQUIRED
# /////////////////////////////////////////////////////////////////////////////
INPUT_DIR1="${BASE_DIR}/data/control/haplotype_caller_output"
INPUT_DIR2="${BASE_DIR}/data/dataset1/haplotype_caller_output_rename"
INPUT_DIR3="${BASE_DIR}/data/dataset2/haplotype_caller_output_rename"
OUTPUT_DIR="${WORK_DIR}/genomics_db_import_output"

# Prepare input
declare -a NUMBER
for j in {1..22} X Y M; do NUMBER+=($j); done
INDEX=${NUMBER[$SLURM_ARRAY_TASK_ID]}
# For single dir:
# gvcf_files=$(for f in ${INPUT_DIR}/*.g.vcf.gz; do echo -n "-V $f " ; done)

# For multiple dir:
gvcf_files=$(find ${INPUT_DIR1} ${INPUT_DIR2} ${INPUT_DIR3} -name "*.g.vcf.gz" -print0 | xargs -0 -I {} echo -n "-V {} ")

# Run GenomicsDBImport
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms20G -Xmx20G -XX:ParallelGCThreads=2" GenomicsDBImport \
	--genomicsdb-workspace-path ${OUTPUT_DIR}/chr${INDEX}_gdb \
	-R $REF \
	$gvcf_files \
	--tmp-dir ${SCRATCH_DIR}/temp/${SLURM_JOBID} \
	--intervals chr${INDEX} \
	--genomicsdb-shared-posixfs-optimizations true

echo "FINISH AT $(date)"
