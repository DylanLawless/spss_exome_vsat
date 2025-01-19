#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 20G
#SBATCH --time 01:00:00
#SBATCH --job-name=pre_annotation
#SBATCH --output=./log/pre_annotation/pre_annotation_%J.out
#SBATCH --error=./log/pre_annotation/pre_annotation_%J.err
#SBATCH --comment="scitas_cost"

set -e
echo "START AT $(date)"

# module load intel/19.0.5
# module load htslib/1.10.2
module load intel/2021.6.0 # jed
module load htslib/1.14 # jed

# bcftools; follow build from github:
#	http://samtools.github.io/bcftools/howtos/install.html
# vcftools; follow build from github:
#	https://github.com/vcftools/vcftools
#	but use:  ./configure prefix=$HOME
# vt; follow build from github:
#	https://genome.sph.umich.edu/wiki/Vt

# Tools
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"
BCFTOOLS="/home/lawless/bcftools/bcftools"
VCFTOOLS="/home/lawless/bin/vcftools"
VT="/work/gr-fe/lawless/tool/vt/vt"

# Path
SCRATCH_DIR="/scratch/lawless/spss_exome/dataset2"
INPUT_DIR=${SCRATCH_DIR}/genotype_refinement_output
OUTPUT_DIR=${SCRATCH_DIR}/pre_annotation_output

# ref
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

cd ${OUTPUT_DIR}

printf "bcftools filter QUAL DP GQ"
# Remove Batch Effect
$BCFTOOLS filter -i 'QUAL>=30 & INFO/DP>=20 & FORMAT/DP>=10 & FORMAT/GQ>=20' -S . -Ov \
	-o bcftools.vcf --threads 20 \
	${INPUT_DIR}/indel_recal_95_snp_recal_99.7_refined_GQ20.vcf.gz

printf "bgzip and index tbi"
bgzip -c bcftools.vcf > bcftools.vcf.gz
tabix -p vcf bcftools.vcf.gz
rm bcftools.vcf 

printf "GATK exclude filtered"
# Exclude filtered
$GATK --java-options "-Djava.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms18G -Xmx18G -XX:ParallelGCThreads=20" \
        SelectVariants \
        -V ${OUTPUT_DIR}/bcftools.vcf.gz \
        -R $REF \
        --exclude-filtered \
        -exclude-non-variants --remove-unused-alternates \
        -O ${OUTPUT_DIR}/bcftools_gatk.vcf.gz


printf "VT decompose and normalise"
# Normalize and break multiallelic
$VT decompose -s -o temp.vcf bcftools_gatk.vcf.gz
$VT normalize -r $REF -o bcftools_gatk_norm.vcf temp.vcf

# Uncomment the follwing line if you want to use GRCh99 build of SNPEFF. If you use GRCh38_no_alt, no need to change anything!
#cat bcftools_gatk_norm.vcf | sed "s/^chrM/MT/" > bcftools_gatk_norm_corrected.vcf

printf "final bgzip and index tbi"
bgzip -c bcftools_gatk_norm.vcf > bcftools_gatk_norm.vcf.gz
tabix -p vcf bcftools_gatk_norm.vcf.gz
rm temp.vcf
rm bcftools_gatk_norm.vcf


echo "FINISH AT $(date)"
