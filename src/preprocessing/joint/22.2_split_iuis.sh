#!/bin/bash

set -e

# Master variables
source ./variables.txt

# Path
OUTPUT_DIR="${WORK_DIR}/annotation"

BCFTOOLS="/home/lawless/bcftools/bcftools"
VCFTOOLS="/home/lawless/bin/vcftools"

# Add an array with the mainfile names
mainfiles=(
"bcftools_gatk_norm_maf01.recode_vep_conda_impact_iuis_gnomad_af1"
)

for mainfile in "${mainfiles[@]}"; do
    vcf_in=${OUTPUT_DIR}/${mainfile}.vcf
    vcf_out_stem=${OUTPUT_DIR}/${mainfile}_chr

    printf "Processing %s\n" "$mainfile"
    printf "bgzip\n"

    # compress
    bgzip $vcf_in

    # index (tbi)
    printf "index\n"
    $BCFTOOLS index -t ${vcf_in}.gz

    # Split
    printf "bcftools splitting vcf by chr\n"
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
done
printf "\nAll input files processed.\n"
