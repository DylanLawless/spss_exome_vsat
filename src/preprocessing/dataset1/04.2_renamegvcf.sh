#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 24
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 00:02:00
#SBATCH --job-name=rename
#SBATCH --output=./log/hc/rename_%J.out
#SBATCH --error=./log/hc/rename_%J.err

# gunzip $INPUT_DIR/*vcf.gz
# sbatch --array=0-167 04.2_renamegvcf.sh

# module load picard/2.20.8 # fidis
module load picard/2.26.2 # jed

# Path: temp work
SCRATCH_DIR="/scratch/lawless/spss_exome"
INPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output"
OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output_rename"

# Declare
declare -a total_file
for i in ${INPUT_DIR}/*.g.vcf.gz; do total_file+=($i); done
filebase=${total_file[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${filebase} .g.vcf.gz)

picard RenameSampleInVcf \
-INPUT ${INPUT_DIR}/${sample_id}.g.vcf.gz \
-OUTPUT ${OUTPUT_DIR}/${sample_id}_rename.g.vcf.gz \
-NEW_SAMPLE_NAME ${sample_id}

# check names
# perl /work/gr-fe/lawless/tool/vcfhacks/getSampleNames.pl -i ${OUTPUT_DIR}/${sample_id}_rename.g.vcf.gz

# # save and rename the index filename to match
# cp $INPUT_DIR/*tbi $OUTPUT_DIR/
# cd  $OUTPUT_DIR
# for file in *.g.vcf.gz.tbi
# do
# 	cp "$file" "${file%.g.vcf.gz.tbi}_rename.g.vcf.gz.tbi"
# done
# Note that I later needed to rezip files with bgzip and index again, so I am
# note sure if the error was due to this tbi rename, or it was due to gzip 
# somewhere during testing.
