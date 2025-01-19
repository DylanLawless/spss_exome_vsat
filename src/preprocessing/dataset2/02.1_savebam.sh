#!/bin/bash 
# save bams until the project is finished then delete

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
OUTPUT_DIR="${SCRATCH_DIR}/bam" # bams

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset2/data"
BAM_DIR="${WORK_DIR}/bam"

# copy to main storage
cp $OUTPUT_DIR/*  $BAM_DIR/
