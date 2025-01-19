#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 2G
#SBATCH --time 00:10:00
#SBATCH --job-name=tbi
#SBATCH --output=./log/hc/tbi_%J.out
#SBATCH --error=./log/hc/tbi_%J.err
#SBATCH --comment="scitas_cost"

# 0.08 min per file. 10 min per task should be OK. 
# Submit the job as follow:
# sbatch --array=0-167 04.3_index.sh

GATK="/work/gr-fe/saadat/tools/gatk/gatk-4.2.2.0/gatk"
SCRATCH_DIR="/scratch/lawless/spss_exome"
OUTPUT_DIR="${SCRATCH_DIR}/haplotype_caller_output_rename"

# bgzip files instead of gzip
# then remake the tbi with IndexFeatureFile \

cd $OUTPUT_DIR

# Declare
declare -a total_file
for i in *_rename.g.vcf.gz; do total_file+=($i); done
filebase=${total_file[$SLURM_ARRAY_TASK_ID]}
sample_id=$(basename ${filebase} _rename.g.vcf.gz)

$GATK IndexFeatureFile \
	-I ${sample_id}_rename.g.vcf.gz \
	-O ${sample_id}_rename.g.vcf.gz.tbi

# test: sequential method
# for file in *_rename.g.vcf; do bgzip $file; done
# for file in *_rename.g.vcf.gz; do $GATK IndexFeatureFile -I ${file} -O ${file}.tbi; done
