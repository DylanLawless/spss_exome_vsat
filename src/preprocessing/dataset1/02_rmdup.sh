#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 80G
#SBATCH --time 02:00:00
#SBATCH --job-name=rmdup
#SBATCH --output=./log/rmdup/rmdup_%J.out
#SBATCH --error=./log/rmdup/rmdup_%J.err

# each run takes ~10 minutes
# Submit the job as follow:
# sbatch --array=0-169 02_rmdup.sh # note 2 samples removed.
# sbatch --array=0-167 02_rmdup.sh

echo "START AT $(date)"
set -e

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Path
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset1/data"
SCRATCH_DIR="/scratch/lawless/spss_exome"
# INPUT_DIR="${WORK_DIR}/bam" # scratch copy data
INPUT_DIR="${SCRATCH_DIR}/bam" # scratch copy data
OUTPUT_DIR="${SCRATCH_DIR}/rmdup_output"

# Create directories
# mkdir -p ${OUTPUT_DIR}/metrics
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Save file list
cd $INPUT_DIR
ls *.bam | \
	cut -f 1 -d "." \
	> $WORK_DIR/bam/sample_bases.txt

declare -a total_file
for i in ${INPUT_DIR}/*bam; do total_file+=($i); done
filebase=${total_file[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${filebase} .bam)

# Run MarkDupSpark
${GATK} --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms70G -Xmx70G" MarkDuplicatesSpark \
	-I ${INPUT_DIR}/${sample_id}.bam \
	-O ${OUTPUT_DIR}/${sample_id}_rmdup.bam \
	--spark-master local[*] \
	-M ${OUTPUT_DIR}/${sample_id}_rmdup_metrics.txt \
	--remove-all-duplicates

echo "END AT $(date)"
