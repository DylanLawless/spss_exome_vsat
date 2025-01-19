#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 20G
#SBATCH --time 02:30:00
#SBATCH --job-name=hc
#SBATCH --output=./log/hc/hc_%J.out
#SBATCH --error=./log/hc/hc_%J.err


# Modified from Ali: <https://github.com/AliSaadatV/Whole-Exome-Genome-Variant-Calling>
# Submit the job as follow:
# sbatch --array=0-379 04_haplotype_caller.sh
# 4GB-6GB bams requires 30-40 minutes. Max bam size is 11GB

echo "START AT $(date)"
set -e

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset2/data"

# Target list - the fastq files were acconpanied by the customer report which lists: twist_humancorePlusRefSeq_hg38. Therefore the panel bed file was downloaded from TwistBioscience to match. Including -L speeds up x5-10 times.
INTERVAL_LIST="$WORK_DIR/twist_humancorePlusRefSeq/Twist_Exome_RefSeq_targets_hg38.bed"

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
INPUT_DIR="${SCRATCH_DIR}/bqsr"
OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output"

# Path: ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
KNOWN_SITES="/work/gr-fe/saadat/pri/known_sites"

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Declare
declare -a total_file
for i in ${INPUT_DIR}/*bam; do total_file+=($i); done
filebase=${total_file[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${filebase} _rmdup_recal.bam)

# Run HC
${GATK} --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms18G -Xmx18G -XX:ParallelGCThreads=2" HaplotypeCaller \
	-I ${INPUT_DIR}/${sample_id}_rmdup_recal.bam \
	-R ${REF} \
	-ERC GVCF \
	-L ${INTERVAL_LIST} \
	-ip 30 \
	-O ${OUTPUT_DIR}/${sample_id}.g.vcf.gz \
	-D ${KNOWN_SITES}/dbsnp_146.hg38.vcf.gz

# Arguments
# -ERC = --emit-ref-confidence
# -L  = --intervals
# -ip = --interval-padding

echo "END AT $(date)"
