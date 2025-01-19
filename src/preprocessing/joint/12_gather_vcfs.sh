#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 36G
#SBATCH --time 03:00:00
#SBATCH --job-name=gather_gvcf
#SBATCH --output=./log/gather_gvcf/gather_gvcf_%J.out
#SBATCH --error=./log/gather_gvcf/gather_gvcf_%J.err
#SBATCH --comment="scitas_cost"

echo "START AT $(date)"
set -e

# Approx < 10 min 2 thread for ~400 samples

 # Master variables
source variables.txt

# Tools
PICARD="/work/gr-fe/saadat/tools/picard/picard.jar"

# Path
OUTPUT_DIR="${WORK_DIR}/genotype_gvcf_output"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Gather VCFs
input_list=$( for i in {1..22} X Y M; do echo -n "-I ${OUTPUT_DIR}/chr${i}.g.vcf " ; done)
java -Djva.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms32G -Xmx32G -XX:ParallelGCThreads=8 -jar $PICARD GatherVcfs \
        $input_list \
        -O ${OUTPUT_DIR}/merged.vcf.gz

echo "FINISH AT $(date)"
