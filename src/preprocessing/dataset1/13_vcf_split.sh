#!/bin/bash

set -e
source ./variables.txt

# Path
OUTPUT_DIR="${WORK_DIR}/annotation"
FILENAME="bcftools_gatk_norm_vep_conda_plug_maf.recode"

BCFTOOLS="/home/lawless/bcftools/bcftools"
VCFTOOLS="/home/lawless/bin/vcftools"

vcf_in=${OUTPUT_DIR}/${FILENAME}.vcf
vcf_out_stem=${OUTPUT_DIR}/${FILENAME}_chr
	
# Split by chr
printf "bgzip\n"
# compress
bgzip $vcf_in

printf "index\n"
# index (tbi)
$BCFTOOLS index -t ${vcf_in}.gz

printf "bcftools splitting vcf by chr\n"
# Split 
for i in {1..22}
do
	$BCFTOOLS view ${vcf_in}.gz --regions chr${i} -o ${vcf_out_stem}_${i}.vcf.gz -Oz
	$BCFTOOLS index -t ${vcf_out_stem}_${i}.vcf.gz
done

for i in X
do
	$BCFTOOLS view ${vcf_in}.gz --regions chr${i} -o ${vcf_out_stem}_${i}.vcf.gz -Oz
	$BCFTOOLS index -t ${vcf_out_stem}_${i}.vcf.gz
done

for i in Y
do
	$BCFTOOLS view ${vcf_in}.gz --regions chr${i} -o ${vcf_out_stem}_${i}.vcf.gz -Oz
	$BCFTOOLS index -t ${vcf_out_stem}_${i}.vcf.gz
done

printf "done\n"
