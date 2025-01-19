#!/bin/bash 

SKIP FOR CONROL
# Path: processing
SCRATCH_DIR="/scratch/lawless/spss_exome"
INPUT_DIR="${SCRATCH_DIR}/fastqcat" # main data
OUTPUT_DIR="${SCRATCH_DIR}/bam"

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset1/data"
METRICS_DIR="${WORK_DIR}/fastpmetrics"
BAM_DIR="${WORK_DIR}/bam"

# move fq to their own dir (this should have been done originally instead of as clean up)
# mkdir $SCRATCH_DIR/fq
mv ${OUTPUT_DIR}/*.fq.gz $SCRATCH_DIR/fq/

cp ${OUTPUT_DIR}/metrics/*_fastp.* ${METRICS_DIR}/
# cp ${OUTPUT_DIR}/*.ba* ${BAM_DIR}/
