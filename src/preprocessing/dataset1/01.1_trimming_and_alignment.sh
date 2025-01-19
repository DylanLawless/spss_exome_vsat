#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 50G
#SBATCH --time 03:00:00
#SBATCH --job-name=fastpbwa
#SBATCH --output=./log/fastpbwa/fastpbwa_%J.out
#SBATCH --error=./log/fastpbwa/fastpbwa_%J.err
# #SBATCH --array=0-169

# source files
# /work/gr-fe/archive/sample_repository/Sepsis/FASTQ/raw

# RUN WITH:
#sbatch --array=0-169 01_trimming_and_alignment.sh

# For less than one node use SBATCH --partition=serial
# ls ./*fastq.gz | wc -l
# 760/2 = 380 dataset 2
# 340/2 = 170 datast 1

set -e

echo "START AT $(date)"

# Tools
BWA="/work/gr-fe/lawless/tool/bwa-0.7.17/bwa"
FASTP="/work/gr-fe/lawless/tool/fastp"

# Modules
# module load intel/19.0.5
module load intel
# module load samtools/1.10
module load samtools/1.14

# Path: processing
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
SCRATCH_DIR="/scratch/lawless/spss_exome"
INPUT_DIR="${SCRATCH_DIR}/fastqcat" # main data
OUTPUT_DIR="${SCRATCH_DIR}/bam"

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset1/data"
# METRICS_DIR="${WORK_DIR}/fastpmetrics"
BAM_DIR="${WORK_DIR}/bam"

# Create directories
mkdir -p $OUTPUT_DIR/metrics

# # rename file to match format
# cd $INPUT_DIR

# for file in *.R1.fastq.gz
# do
#   mv "$file" "${file%.R1.fastq.gz}_R1.fastq.gz"
# done

# for file in *.R2.fastq.gz
# do
#   mv "$file" "${file%.R2.fastq.gz}_R2.fastq.gz"
# done

# Extract fastq files
declare -a total_R1
for i in ${INPUT_DIR}/*R1.fastq.gz; do total_R1+=($i); done
R1=${total_R1[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${R1} _R1.fastq.gz)
R2=${INPUT_DIR}/${sample_id}_R2.fastq.gz

# Trimming, alignment, and converting to bam
$FASTP \
	--in1 ${INPUT_DIR}/${sample_id}_R1.fastq.gz \
	--in2 ${INPUT_DIR}/${sample_id}_R2.fastq.gz \
	-j "${OUTPUT_DIR}/metrics/${sample_id}_fastp.json" \
	-h "${OUTPUT_DIR}/metrics/${sample_id}_fastp.html" \
	--out1 ${OUTPUT_DIR}/${sample_id}_R1.fq.gz \
	--out2 ${OUTPUT_DIR}/${sample_id}_R2.fq.gz \
	--thread 28

$BWA mem \
	${REF} \
	${OUTPUT_DIR}/${sample_id}_R1.fq.gz \
	${OUTPUT_DIR}/${sample_id}_R2.fq.gz \
	-v 1 -M \
	-t 28 \
	-R '@RG\tID:$sample_id\tSM:$sample_id\tPL:NovaSeq6000twistexome\tLB:$sample_id' |\
	samtools view \
	--threads 28 \
	-O BAM -o ${OUTPUT_DIR}/${sample_id}.bam
# NOTE: Check files have taken the variable literally instead of printing sample name.
# can correct in vcf file with picard.
# I believe it should be:
# -R '@RG\tID:${sample_id} --- NO
# -R '@RG\tID:$sample_id --- Yes

# # move fq to their own dir (this should have been done originally instead of as clean up)
# mkdir $SCRATCH_DIR/fq
# mv ${OUTPUT_DIR}/*.fq.gz $SCRATCH_DIR/fq/

# cp ${OUTPUT_DIR}/metrics/${sample_id}_fastp.* ${METRICS_DIR}/
# cp ${OUTPUT_DIR}/${sample_id}.ba* ${BAM_DIR}/

echo "END AT $(date)"

# 2 samples failed:
# raw_BE059 - fastq are 14G each, 7x times normal size and probably error. Skip this file.
#	fastpbwa_1147183.err
# raw_CH064 - rerun single line command produced the same error: invalid compressed data--format violated.
#	fastpbwa_1147209.err
# Since the original file in archive has the same md5sum we must assume that the stored sample is currupted. 
# maybe there is an older backup available for this file.

