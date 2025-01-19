#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G
#SBATCH --time 03:00:00
#SBATCH --job-name=maxmaf
#SBATCH --output=./log/vep/maxmaf_%J.out
#SBATCH --error=./log/vep/maxmaf_%J.err
#SBATCH --comment="scitas_cost"

set -e
echo "START AT $(date)"
source ./variables.txt
module load intel # jed
# module load htslib/1.14 # jed
VCFTOOLS="/home/lawless/bin/vcftools"

OUTPUT_DIR="${WORK_DIR}/annotation"
cd ${OUTPUT_DIR}
FILENAME="bcftools_gatk_norm_vep_conda_plug"

# Filter on MAF
printf "\nVCFtools filter max-maf\n"
echo "$(date)"
$VCFTOOLS --vcf ${FILENAME}.vcf \
	--max-maf 0.02 \
	--recode --recode-INFO-all \
	--out ${FILENAME}_maf
# Exclude on minor allele frequency. max-maf 0.02*168 samples = max 3.8 alleles (<2 hom or <4 het)


# printf "\nfinal bgzip and index tbi\n"
# echo "$(date)"
# bgzip -c ${FILENAME}_maf.recode.vcf > ${FILENAME}_maf.recode.vcf.gz
# tabix -p vcf ${FILENAME}_maf.recode.vcf.gz
# rm temp.vcf
# rm ${FILENAME}_maf

echo "END AT $(date)"
