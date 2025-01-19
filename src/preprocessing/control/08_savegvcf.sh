#!/bin/bash

echo "START AT $(date)"
set -e

# Master variables
source ./variables.txt

# /////////////////////////////////////////////////////////////////////////////
# CHECK IF VCF RENAME WAS REQUIRED
# CHECK IF VCF RENAME WAS REQUIRED
# CHECK IF VCF RENAME WAS REQUIRED
# /////////////////////////////////////////////////////////////////////////////

OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output"
# OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output_rename"

cp -r ${OUTPUT_DIR} ${WORK_DIR}/

echo "END AT $(date)"
