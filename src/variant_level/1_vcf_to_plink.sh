#!/bin/bash

echo "START AT $(date)"

set -e

source ./variables.txt

INPUT_DIR="${WORK_DIR}/annotation"
cd  ${VARLEVEL_DIR}

# Unzip the VCF to "file.vcf" and keep the original.
gunzip -c \
	${INPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small_impact_gnomad.vcf.gz > file.vcf

wc -l file.vcf # count total lines
grep "^#" file.vcf | wc -l # count header lines

# ~/tools/plink1.9/
plink --vcf file.vcf --make-bed --out output1 --vcf-half-call m --double-id

# remove vcf
rm file.vcf

