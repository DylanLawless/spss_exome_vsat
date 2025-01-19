#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 04:00:00
#SBATCH --job-name=genotype_refine
#SBATCH --output=./log/genotype_refine/genotype_refine_%J.out
#SBATCH --error=./log/genotype_refine/genotype_refine_%J.err
#SBATCH --comment="scitas_cost"

# runtime for 1000 samples ~2 hours

echo "START AT $(date)"
set -e

# Master variables
source variables.txt

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Path: temp work
INPUT_DIR="${WORK_DIR}/vqsr_output"
OUTPUT_DIR="${WORK_DIR}/genotype_refinement_output"

# Ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
KNOWN_SITES="/work/gr-fe/saadat/pri/known_sites"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Run genotype refinement
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms24G -Xmx24G -XX:ParallelGCThreads=8" \
	CalculateGenotypePosteriors \
	-V ${INPUT_DIR}/indel_recal_95_snp_recal_99.7.vcf.gz \
	-O ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined.vcf.gz \
	--supporting-callsets ${KNOWN_SITES}/af-only-gnomad.hg38.vcf.gz \
	--num-reference-samples-if-no-call 20314 && \

$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms24G -Xmx24G -XX:ParallelGCThreads=8" \
        VariantFiltration \
        -V ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined.vcf.gz \
        -R $REF \
        --genotype-filter-expression "GQ < 20" \
        --genotype-filter-name "lowGQ" \
        -O ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined_GQ20.vcf.gz

echo "FINISH AT $(date)"
