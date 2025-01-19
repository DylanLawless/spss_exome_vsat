#!/bin/bash

set -e

# Master variables
source ./variables.txt

OUTPUT_DIR="${SCRATCH_DIR}/genotype_gvcf_output"

# Save to main storage
cp -r $OUTPUT_DIR $WORK_DIR/
