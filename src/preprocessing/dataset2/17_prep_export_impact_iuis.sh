#!/bin/bash

set -e
source ./variables.txt

# Path
OUTPUT_DIR="${WORK_DIR}/annotation"

BCFTOOLS="/home/lawless/bcftools/bcftools"
VCFTOOLS="/home/lawless/bin/vcftools"

FILENAME="bcftools_gatk_norm_vep_conda_plug_maf.recode"
vcf_in=${OUTPUT_DIR}/${FILENAME}_canonical_impact_iuis_chr

printf "bgzip and index\n"
for i in {1..22}
do
	bgzip ${vcf_in}_${i}.vcf
	$BCFTOOLS index -t ${vcf_in}_${i}.vcf.gz
done 

for i in X
do
	bgzip ${vcf_in}_${i}.vcf
	$BCFTOOLS index -t ${vcf_in}_${i}.vcf.gz
done

for i in Y
do
	bgzip ${vcf_in}_${i}.vcf
	$BCFTOOLS index -t ${vcf_in}_${i}.vcf.gz
done

printf "done\n"
