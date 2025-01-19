#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 25
#SBATCH --cpus-per-task 2
#SBATCH --mem 20G
#SBATCH --time 00:30:00
#SBATCH --job-name=genomics_db_import
#SBATCH --output=./log/genomics_db_import/genomics_db_import_%J.out
#SBATCH --error=./log/genomics_db_import/genomics_db_import_%J.err
#SBATCH --comment="scitas_cost"

# ~10 min per task
# array for all intervals: 1-22 X Y M = 25
# sbatch --array=0-24 05_genomics_db_import.sh

echo "START AT $(date)"
set -e

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset1/data"
INPUT_DIR="${WORK_DIR}/haplotype_caller_output_rename"

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
OUTPUT_DIR="${SCRATCH_DIR}/genomics_db_import_output"

# Path: ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Prepare input
declare -a NUMBER
for j in {1..22} X Y M; do NUMBER+=($j); done
INDEX=${NUMBER[$SLURM_ARRAY_TASK_ID]}
gvcf_files=$(for f in ${INPUT_DIR}/*.g.vcf.gz; do echo -n "-V $f " ; done)

# Run GenomicsDBImport
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms20G -Xmx20G -XX:ParallelGCThreads=2" GenomicsDBImport \
	--genomicsdb-workspace-path ${OUTPUT_DIR}/chr${INDEX}_gdb \
	-R $REF \
	$gvcf_files \
	--tmp-dir ${SCRATCH_DIR}/temp/${SLURM_JOBID} \
	--intervals chr${INDEX} \
	--genomicsdb-shared-posixfs-optimizations true

echo "FINISH AT $(date)"
