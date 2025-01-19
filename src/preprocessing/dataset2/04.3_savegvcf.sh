#!/bin/bash 

echo "START AT $(date)"
set -e

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset2/data"

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output_rename"

# Run HC
cp -r ${OUTPUT_DIR} ${WORK_DIR}

echo "END AT $(date)"
