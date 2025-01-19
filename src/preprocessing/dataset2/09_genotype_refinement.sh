#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 6G
#SBATCH --time 00:30:00
#SBATCH --job-name=genotype_refine
#SBATCH --output=./log/genotype_refine/genotype_refine_%J.out
#SBATCH --error=./log/genotype_refine/genotype_refine_%J.err

# runtime for 380 samples ~10 seconds

echo "START AT $(date)"
set -e

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
INPUT_DIR=${SCRATCH_DIR}/vqsr_output
OUTPUT_DIR=${SCRATCH_DIR}/genotype_refinement_output

# Ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
KNOWN_SITES="/work/gr-fe/saadat/pri/known_sites"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Run genotype refinement
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms5G -Xmx5G -XX:ParallelGCThreads=2" \
	CalculateGenotypePosteriors \
	-V ${INPUT_DIR}/indel_recal_95_snp_recal_99.7.vcf.gz \
	-O ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined.vcf.gz \
	--supporting-callsets ${KNOWN_SITES}/af-only-gnomad.hg38.vcf.gz \
	--num-reference-samples-if-no-call 20314 && \

$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms5G -Xmx5G -XX:ParallelGCThreads=2" \
        VariantFiltration \
        -V ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined.vcf.gz \
        -R $REF \
        --genotype-filter-expression "GQ < 20" \
        --genotype-filter-name "lowGQ" \
        -O ${OUTPUT_DIR}/indel_recal_95_snp_recal_99.7_refined_GQ20.vcf.gz

echo "FINISH AT $(date)"
