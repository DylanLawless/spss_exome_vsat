NOT IN USE
#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 50G
#SBATCH --time 03:00:00
#SBATCH --job-name=fastpbwa
#SBATCH --output=./log/fastpbwa/fastpbwa_%J.out
#SBATCH --error=./log/fastpbwa/fastpbwa_%J.err
# # #SBATCH --array=0-169 set number of samples

# source files
# /work/gr-fe/archive/sample_repository/?

# RUN WITH:
#sbatch --array=0-169 01_trimming_and_alignment.sh

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
	-j "${FQ_METRICS}/${sample_id}_fastp.json" \
	-h "${FQ_METRICS}/${sample_id}_fastp.html" \
	--out1 ${FQ_DIR}/${sample_id}_R1.fq.gz \
	--out2 ${FQ_DIR}/${sample_id}_R2.fq.gz \
	--thread 28

$BWA mem \
	${REF} \
	${FQ_DIR}/${sample_id}_R1.fq.gz \
	${FQ_DIR}/${sample_id}_R2.fq.gz \
	-v 1 -M \
	-t 28 \
	-R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:unknownexome\tLB:${sample_id}" |\
	samtools view \
	--threads 28 \
	-O BAM -o ${OUTPUT_DIR}/${sample_id}.bam

echo "END AT $(date)"

