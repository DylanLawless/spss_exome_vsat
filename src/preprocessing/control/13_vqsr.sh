#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 12G
#SBATCH --time 01:00:00
#SBATCH --job-name=vqsr
#SBATCH --output=./log/vqsr/vqsr_%J.out
#SBATCH --error=./log/vqsr/vqsr_%J.err

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# NOTE
# For idels I required --max-gaussians 4 as the number is low.
# Controls: 287760 ./control/vqsr_output/merged_snp1.recal
# Controls:  18520 ./control/vqsr_output/merged_indel1.recal
# Dataset1: 277800 ./dataset1/vqsr_output/merged_snp1.recal
# Dataset1:  21095 ./dataset1/vqsr_output/merged_indel1.recal
# Therefore, SNP recal is normal, while indel recal is "low" but not concerning. 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# Runtine for 380 samples ~10 min

echo "START AT $(date)"
set -e

# Master variables
source variables.txt

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
KNOWN_SITES="/work/gr-fe/saadat/pri/known_sites"

# Path: temp work
INPUT_DIR="${SCRATCH_DIR}/genotype_gvcf_output"
OUTPUT_DIR="${SCRATCH_DIR}/vqsr_output"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Run VQSR
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms5G -Xmx5G -XX:ParallelGCThreads=4" VariantRecalibrator \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 \
	-tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
	-tranche 95.0 -tranche 94.0 \
	-tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-R $REF \
	-V ${INPUT_DIR}/merged.vcf.gz \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
	${KNOWN_SITES}/hapmap_3.3.hg38.vcf.gz \
	--resource:omni,known=false,training=true,truth=false,prior=12.0 \
	${KNOWN_SITES}/1000G_omni2.5.hg38.vcf.gz \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 \
	${KNOWN_SITES}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP -O ${OUTPUT_DIR}/merged_snp1.recal --tranches-file ${OUTPUT_DIR}/output_snp1.tranches && \

$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms5G -Xmx5G -XX:ParallelGCThreads=2" VariantRecalibrator \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 \
	-tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
	-tranche 95.0 -tranche 94.0 \
	-tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-R $REF \
	-V ${INPUT_DIR}/merged.vcf.gz \
	--resource:mills,known=false,training=true,truth=true,prior=12.0 \
	${KNOWN_SITES}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
	${KNOWN_SITES}/dbsnp_146.hg38.vcf.gz \
	-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode INDEL -O ${OUTPUT_DIR}/merged_indel1.recal --tranches-file ${OUTPUT_DIR}/output_indel1.tranches \
	--max-gaussians 4 && \

# Apply VQSR
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms5G -Xmx5G -XX:ParallelGCThreads=4" ApplyVQSR \
	-V ${INPUT_DIR}/merged.vcf.gz \
	--recal-file ${OUTPUT_DIR}/merged_snp1.recal \
	-mode SNP \
	--tranches-file ${OUTPUT_DIR}/output_snp1.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	-O ${OUTPUT_DIR}/snp_recal_99.7.vcf.gz && \

$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms5G -Xmx5G -XX:ParallelGCThreads=4" ApplyVQSR \
	-V ${OUTPUT_DIR}/snp_recal_99.7.vcf.gz \
	--recal-file ${OUTPUT_DIR}/merged_indel1.recal \
	-mode INDEL \
	--tranches-file ${OUTPUT_DIR}/output_indel1.tranches \
	--truth-sensitivity-filter-level 95 \
	--create-output-variant-index true \
	-O ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7.vcf.gz

echo "FINISH AT $(date)"
