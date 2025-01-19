#!/bin/bash 
#SBATCH --array=0-391
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 12G
#SBATCH --time 04:00:00
#SBATCH --job-name=bqsr
#SBATCH --output=./log/bqsr/bqsr_%J.out
#SBATCH --error=./log/bqsr/bqsr_%J.err
#SBATCH --comment="scitas_cost"

# each run takes ~1.5 hours.
# the total time for ~400 controls was ~1 hour.
# (380 case samples was almost 7 hours but paramerts were not optmised).
# according to the docs from biowolf, increasing threads beyond 4 does not help much <https://hpc.nih.gov/training/gatk_tutorial/bqsr.html>
# try cpus = 4. on cpu = 20 it took 5.9 minutes to reach chr2 on log file. On the controls it seems to take half this time. 

echo "START AT $(date)"
set -e

# Master variables
source variables.txt

# Path: temp work
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
${GATK} --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms4G -Xmx4G -XX:ParallelGCThreads=4" BaseRecalibrator \
	-I ${INPUT_DIR}/${sample_id}_rmdup.bam \
	-R ${REF} \
	-O ${OUTPUT_DIR}/${sample_id}_bqsr.table \
	--known-sites ${KNOWN_SITES}/dbsnp_146.hg38.vcf.gz \
	--known-sites ${KNOWN_SITES}/Homo_sapiens_assembly38.known_indels.vcf.gz \
	--known-sites ${KNOWN_SITES}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz && \

${GATK} --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms4G -Xmx4G -XX:ParallelGCThreads=4" ApplyBQSR \
	-I ${INPUT_DIR}/${sample_id}_rmdup.bam \
	-R $REF \
	--bqsr-recal-file ${OUTPUT_DIR}/${sample_id}_bqsr.table \
	-O ${OUTPUT_DIR}/${sample_id}_rmdup_recal.bam

echo "END AT $(date)"

# For posterity
# lawless@jed:/scratch/lawless/spss_exome$ ls /work/gr-fe/saadat/pri/known_sites/
# 1000G_omni2.5.hg38.vcf.gz
# 1000G_omni2.5.hg38.vcf.gz.tbi
# 1000G_phase1.snps.high_confidence.hg38.vcf.gz
# 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
# af-only-gnomad.hg38.vcf.gz
# af-only-gnomad.hg38.vcf.gz.tbi
# Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
# Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
# dbsnp_146.hg38.vcf.gz
# dbsnp_146.hg38.vcf.gz.tbi
# hapmap_3.3.hg38.vcf.gz
# hapmap_3.3.hg38.vcf.gz.tbi
# HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
# HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi
# Homo_sapiens_assembly38.known_indels.vcf.gz
# Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
# known_sites_GRCH38
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
# README.txt
