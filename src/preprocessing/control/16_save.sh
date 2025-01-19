#!/bin/bash
# Master variables
source ./variables.txt

DIR_1="${SCRATCH_DIR}/genotype_gvcf_output"
DIR_2="${SCRATCH_DIR}/genotype_refinement_output"
DIR_3="${SCRATCH_DIR}/pre_annotation_output"

# Save to main storage
cp -r $DIR_1 $WORK_DIR/
cp -r $DIR_2 $WORK_DIR/
cp -r $DIR_3 $WORK_DIR/
