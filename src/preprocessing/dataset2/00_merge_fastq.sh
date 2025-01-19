#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 3
#SBATCH --cpus-per-task 3
#SBATCH --mem 60G
#SBATCH --time 04:00:00
#SBATCH --mail-user=dylan.lawless@epfl.ch
#SBATCH --mail-type=END 
#SBATCH --job-name=alignment_no_alt
#SBATCH --output=./log/fastqcat/fastqcat.out
#SBATCH --error=./log/fastqcat/fastqcat.err

set -e
# cut is used to split id1_L001_R1_001.fastq.gz by _ to get the sample ID.

WORK_DIR="/work/gr-fe/lawless/spss/exome/dataset2/data"
INPUT_DIR="${WORK_DIR}/fastq"
SCRATCH_DIR="/scratch/lawless/spss_exome"
OUTPUT_DIR="${SCRATCH_DIR}/fastqcat"
# INPUT_DIR="${SCRATCH_DIR}/testdata"

echo "START AT $(date)"

# Note that the cut options are specific to this dataset
# Since the data is provided in 3 dir, I will interate through them instead of rearranging.
# Source:
# /work/backup/gr-fe/lawless/spss/exome/dataset2/WES_Schlapbach_Oct21
# 220405_A00485_0265_AHC7F5DMXY
# 220407_A00485_0267_AHGHY3DSX3
# 220414_A00485_0270_AHH2H5DSX3

# Parallel runs
for dir in ${INPUT_DIR}/2204*; do cd $dir

	ls *_R1_* | cut -d _ -f 1-5 | sort | uniq \
	| while read id; do \
		cat $id*_R1_*.fastq.gz > $OUTPUT_DIR/$id\_R1.fastq.gz;
		cat $id*_R2_*.fastq.gz > $OUTPUT_DIR/$id\_R2.fastq.gz;
	done

done &
wait

echo "END AT $(date)"
# 380 samples present

# # Sanity check
# # Before 
# touch ${OUTPUT_DIR}/sanity_cat_before.txt
# wc -l ${INPUT_DIR}/2204*/*_R1_* >> ${OUTPUT_DIR}/sanity_cat_before.txt

# # After 
# touch ${OUTPUT_DIR}/sanity_cat_after.txt
# wc -l ${INPUT_DIR}/2204*/*_R1_* >> ${OUTPUT_DIR}/sanity_cat_after.txt

# # check: should be the same line count, if all samples IDs are cut -d -f1-5.

