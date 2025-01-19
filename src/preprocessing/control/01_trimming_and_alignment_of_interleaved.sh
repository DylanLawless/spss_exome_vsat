#!/bin/bash 
#SBATCH --array=0-391
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 50G
#SBATCH --time 02:00:00
#SBATCH --job-name=fastpbwa
#SBATCH --output=./log/fastpbwa/fastpbwa_%J.out
#SBATCH --error=./log/fastpbwa/fastpbwa_%J.err
#SBATCH --comment="scitas_cost"

# 392 interleaved fq:  --array=0-391
# source files see 00_set_up.sh
# RUN WITH sbatch since --array is set

set -e
echo "START AT $(date)"

# Master variables
source variables.txt

# Tools
BWA="/work/gr-fe/lawless/tool/bwa-0.7.17/bwa"
FASTP="/work/gr-fe/lawless/tool/fastp"

# Modules
module load intel
module load samtools/1.14

# Path: processing
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
INPUT_DIR="${SCRATCH_DIR}/fastq_raw" # main data
FQ_DIR="${SCRATCH_DIR}/fq" # 
FQ_METRICS="${SCRATCH_DIR}/fq_metrics" # 
OUTPUT_DIR="${SCRATCH_DIR}/bam"

# # rename file to match format
# cd $INPUT_DIR

# If you need to rename files
# for file in *.R1.fastq.gz
# do
#   mv "$file" "${file%.R1.fastq.gz}_R1.fastq.gz"
# done

# Extract interleaved fastq files. 
declare -a total_interleaved
for i in ${INPUT_DIR}/*interleaved.fq.gz; do total_interleaved+=($i); done
interleaved=${total_interleaved[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${interleaved} .interleaved.fq.gz)

printf "\n Running Fastp\n"
# Trimming, alignment, and converting to bam
$FASTP \
	--interleaved_in \
	--in1 ${INPUT_DIR}/${sample_id}.interleaved.fq.gz \
	-j "${FQ_METRICS}/${sample_id}_fastp.json" \
	-h "${FQ_METRICS}/${sample_id}_fastp.html" \
	--out1 ${FQ_DIR}/${sample_id}_R1.fq.gz \
	--out2 ${FQ_DIR}/${sample_id}_R2.fq.gz \
	--thread 28

printf "\n Running BWA\n"
# For BWA -p is used for interleaved
$BWA mem \
	${REF} \
	${FQ_DIR}/${sample_id}_R1.fq.gz \
	${FQ_DIR}/${sample_id}_R2.fq.gz \
	-p \
	-v 1 -M \
	-t 28 \
	-R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:unknownexome\tLB:${sample_id}" |\
	samtools view \
	--threads 28 \
	-O BAM -o ${OUTPUT_DIR}/${sample_id}.bam

echo "END AT $(date)"

