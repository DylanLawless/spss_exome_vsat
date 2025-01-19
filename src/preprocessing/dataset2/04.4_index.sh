#!/bin/bash
GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset2/data"
INPUT_DIR="${WORK_DIR}/haplotype_caller_output_rename"

# bgzip files instead of gzip
# then remake the tbi with IndexFeatureFile \

cd $INPUT_DIR
# gunzip ./*.g.vcf.gz
# for file in *_rename.g.vcf; do bgzip $file; done
for file in *_rename.g.vcf.gz; do $GATK IndexFeatureFile -I ${file} -O ${file}.tbi; done
