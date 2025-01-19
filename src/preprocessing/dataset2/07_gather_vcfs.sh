#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 20G
#SBATCH --time 01:00:00
#SBATCH --job-name=gather_gvcf
#SBATCH --output=./log/gather_gvcf/gather_gvcf_%J.out
#SBATCH --error=./log/gather_gvcf/gather_gvcf_%J.err

echo "START AT $(date)"
set -e

# Approx < 10 min 2 thread.

 # Master variables
source variables.txt

# Tools
PICARD="/work/gr-fe/saadat/tools/picard/picard.jar"

# Path
OUTPUT_DIR="${WORK_DIR}/genotype_gvcf_output"

# Create directories
mkdir -p ${SCRATCH_DIR}/temp/${SLURM_JOBID}/io

# Gather VCFs
cd $OUTPUT
input_list=$( for i in {1..22} X Y M; do echo -n "-I ${OUTPUT_DIR}/chr${i}.g.vcf " ; done)
java -Djva.io.tmpdir=${SCRATCH_DIR}/temp/${SLURM_JOBID}/io -Xms18G -Xmx18G -XX:ParallelGCThreads=2 -jar $PICARD GatherVcfs \
        $input_list \
        -O ${OUTPUT_DIR}/merged.vcf.gz

echo "FINISH AT $(date)"
