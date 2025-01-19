#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 12G
#SBATCH --time 04:00:00
#SBATCH --job-name=bqsr
#SBATCH --output=./log/bqsr/bqsr_%J.out
#SBATCH --error=./log/bqsr/bqsr_%J.err

# Modified from Ali: <https://github.com/AliSaadatV/Whole-Exome-Genome-Variant-Calling>
# Submit the job as follow:
# sbatch --array=0-379 03_bqsr.sh

echo "START AT $(date)"
set -e

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset2/data"

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
INPUT_DIR="${SCRATCH_DIR}/rmdup_output"
OUTPUT_DIR="${SCRATCH_DIR}/bqsr"
TEMP="${SCRATCH_DIR}/temp"

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Path
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
KNOWN_SITES="/work/gr-fe/saadat/pri/known_sites"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

declare -a total_file
for i in ${INPUT_DIR}/*bam; do total_file+=($i); done
filebase=${total_file[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${filebase} _rmdup.bam)

# Run and Apply BQSR
${GATK} --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms4G -Xmx4G -XX:ParallelGCThreads=2" BaseRecalibrator \
	-I ${INPUT_DIR}/${sample_id}_rmdup.bam \
	-R ${REF} \
	-O ${OUTPUT_DIR}/metrics/${sample_id}_bqsr.table \
	--known-sites ${KNOWN_SITES}/dbsnp_146.hg38.vcf.gz \
	--known-sites ${KNOWN_SITES}/Homo_sapiens_assembly38.known_indels.vcf.gz \
	--known-sites ${KNOWN_SITES}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz && \

${GATK} --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms4G -Xmx4G -XX:ParallelGCThreads=2" ApplyBQSR \
	-I ${INPUT_DIR}/${sample_id}_rmdup.bam \
	-R $REF \
	--bqsr-recal-file ${OUTPUT_DIR}/metrics/${sample_id}_bqsr.table \
	-O ${OUTPUT_DIR}/${sample_id}_rmdup_recal.bam

echo "END AT $(date)"

