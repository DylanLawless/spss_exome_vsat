#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 24
# #SBATCH --mem 50G
#SBATCH --mem-per-cpu=1G
#SBATCH --time 01:00:00
#SBATCH --job-name=setup
#SBATCH --output=./setup_log_%J.out
#SBATCH --error=./setup_log_%J.err
#SBATCH --comment="scitas_cost"

# Master variables
source ./variables.txt

# work
mkdir ${WORK_DIR}
cd ${WORK_DIR}
mkdir bam annotation

# src log
cd ${SRC_DIR}
mkdir log
cd log
mkdir bqsr fastpbwa gather_gvcf genomics_db_import genotype_gvcf genotype_refine hc pre_annotation rmdup vep vqsr

# scratch
mkdir ${SCRATCH_DIR}
cd ${SCRATCH_DIR}

mkdir bam bqsr fastq_raw fq fq_metrics genomics_db_import_output genotype_gvcf_output genotype_refinement_output haplotype_caller_output haplotype_caller_output_rename pre_annotation_output rmdup_output temp vqsr_output

# File prep
# file source
# ls /work/gr-fe/archive/sample_repository/SHCS392exomes/FASTQ/
# 392 interleaved fq
# IMPORTANT !!! For this set of exomes original FASTQ files were not available. Only BAMS were available. So the general steps were: ORIGINAL BAMs --> INTERLEAVED FASTQs --> UNMAPPED BAMs --> MAPPED BAMs --> GVCFs --> VCFs
# (when FASTQs are avilable the steps are: MultipleFASTQfiles --> rawFASTQs --> Filtered/Trimmed R1/R2 FASTQs --> UNMAPPED BAMs --> then as above) 

printf "Copying controls to scratch.\n"
# Set the source and destination directories
dir1="/work/gr-fe/archive/sample_repository/SHCS392exomes/FASTQ/"
dir2="${SCRATCH_DIR}/fastq_raw"

# Find all files in the source directory
find "$dir1" -maxdepth 1 -type f |

# Pass the file names to xargs, and run up to 24 copy processes at a time
xargs -P 24 -I {} cp {} "$dir2"

printf "Are fastq_raw now prepped?. Total cound and head:\n"
ls fastq_raw | wc -l
ls fastq_raw | head

printf "\nThe variables are set for src, work, scratch:\n"
ls ${SRC_DIR}
ls ${SCRATCH_DIR}
ls ${WORK_DIR}

cd ${SRC_DIR}
mv setup_log* ./log/
