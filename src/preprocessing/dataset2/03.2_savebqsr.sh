#!/bin/bash 
# Usually I would save raw bams until the project is finished then delete, 
# as per 02.1_savebam.sh.
# I am saving the rmdup_output instead as I want the .bai (index) from bams. 
# I want to check bam/bai with IGV visually, which requires the bai file also.

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
OUTPUT_DIR="${SCRATCH_DIR}/rmdup_output"

# Path: storage
WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset2/data"

# copy to main storage
# cp -r $OUTPUT_DIR  $WORK_DIR/
printf("Skipping bam save. Storage requirements are high - only use if pausing pipeline with data on scratch.")
