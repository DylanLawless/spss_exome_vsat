#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 32
#SBATCH --cpus-per-task 1
#SBATCH --mem 32G
#SBATCH --time 04:00:00
#SBATCH --job-name=pca
#SBATCH --output=./log/pca/pca_%J.out
#SBATCH --error=./log/pca/pca_%J.err
#SBATCH --comment="scitas_cost"

set -e
echo "START AT $(date)"

# Master variables
source ./variables.txt

module load intel/2021.6.0 # jed
module load htslib/1.14 # jed
# Path
INPUT_DIR="${WORK_DIR}/pre_annotation_output"
OUTPUT_DIR="${WORK_DIR}/pca_output"

#==============================================================================
# VCF to plink format
#==============================================================================

# double-id: causes both family and within-family IDs to be set to the sample ID.
# vcf-half-call m: Treat half-calls as missing.
# allow-extra-chr: PLINK reports an error if the input data contains unrecognized chromosome codes (such as hg19 haplotype chromosomes or unplaced contigs). 

plink2 \
	--vcf ${INPUT_DIR}/bcftools_gatk_norm.vcf \
	--make-bed \
	--double-id \
	--allow-extra-chr \
	--vcf-half-call m \
	--out ${OUTPUT_DIR}/bcftools_gatk_norm

#==============================================================================
# Make GRM and PCA
# Prune known and cohort-specific LD regions for the grm
#==============================================================================

# Removing long-range LD regions for PCA
plink2 \
	--bfile ${OUTPUT_DIR}/bcftools_gatk_norm \
	--exclude range ${OUTPUT_DIR}/high_LD_regions_build_grch38 \
	--make-bed \
	--out ${OUTPUT_DIR}/bcftools_gatk_norm_no_lrldr

# Produce a pruned subset of markers that are in approximate 
# linkage equilibrium with each other
# Writes the IDs to plink.prune.in
# (and the IDs of all excluded variants to plink.prune.out)
plink2 \
	--bfile ${OUTPUT_DIR}/bcftools_gatk_norm_no_lrldr \
	--indep-pairwise 50 5 0.5 \
	--out ${OUTPUT_DIR}/bcftools_gatk_norm_no_lrldr

# Extract the pruned subset of markers
plink2 \
	--bfile ${OUTPUT_DIR}/bcftools_gatk_norm_no_lrldr \
	--extract ${OUTPUT_DIR}/bcftools_gatk_norm_no_lrldr.prune.in \
	--make-bed --out ${OUTPUT_DIR}/bcftools_gatk_norm_no_lrldr.pruned

# PCA
# this takes apporx 10 minutes on 32 threads, 70 minutes on this dataset on ~18 thread

/work/gr-fe/lawless/tool/gcta64 \
	--bfile ${OUTPUT_DIR}/bcftools_gatk_norm_no_lrldr.pruned \
	--autosome --make-grm \
	--out ${OUTPUT_DIR}/bcftools_gatk_norm_grm\
	--thread-num 32

/work/gr-fe/lawless/tool/gcta64 \
	--grm ${OUTPUT_DIR}/bcftools_gatk_norm_grm \
	--pca 10 \
	--out ${OUTPUT_DIR}/bcftools_gatk_norm_pca \
	--thread-num 32

# have a look at the PCA in gwas_2/grm/plots_pca.pdf
module load gcc/11.3.0
module load r/4.1.3
Rscriptl pca.R


# # exclude 2 outliers
# plink \
#     --bfile $wdir/geno/cohort.filt \
#     --remove $wdir/grm/outliers \
#     --make-bed \
#     --out $wdir/geno/cohort.filtered \
#     --threads 12
# 
