#!/bin/bash
SCRATCH_DIR="/scratch/lawless/spss_exome"
OUTPUT_DIR="${SCRATCH_DIR}/genotype_gvcf_output"
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset1/data"

# Save to main storage
cp -r $OUTPUT_DIR $WORK_DIR/
