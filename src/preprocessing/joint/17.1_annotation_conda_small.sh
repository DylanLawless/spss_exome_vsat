#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 30G
#SBATCH --time 03:00:00
#SBATCH --job-name=vep_small
#SBATCH --output=./log/vep/vep_small_%J.out
#SBATCH --error=./log/vep/vep_small_%J.err
#SBATCH --comment="scitas_cost"

# CONDA ENV
# activate before submission
# conda activate vep

echo "START AT $(date)"
set -e
source ./variables.txt
module load intel # jed

# Path
INPUT_DIR="${WORK_DIR}/pre_annotation_output"
OUTPUT_DIR="${WORK_DIR}/annotation"

# db
DATABASES="/work/gr-fe/databases"
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
GTF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/gtf_and_gff/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz"
GFF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/gtf_and_gff/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation_for_vep.gff.gz"
PLUGINS="${DATABASES}/vep_hg38/plugins"
LOFTEE="${DATABASES}/vep_hg38/loftee"
CACHE="${DATABASES}/vep_hg38/cache_hg38_ensembl"
CLINVAR="${DATABASES}/clinvar/hg38/20220730/clinvar_20220730.vcf.gz"
CADD_DIR="${DATABASES}/CADD_v1.6_GRCH38_hg38"
G2P_PANNEL="/work/gr-fe/saadat/panels/Primary_immunodeficiency.csv"
dbNSFP="/work/gr-fe/lawless/database/dbNSFP/data/dbNSFP4.2a_grch38.gz"
IntAct="/work/gr-fe/lawless/database/intact"
# REVEL="/work/gr-fe/lawless/ref/revel/new_tabbed_revel_grch38.tsv.gz"
# GNOMAD="/scratch/lawless/gnomad/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
# GNOMAD="${DATABASES}/gnomAD_GRCh38_sites/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
# GNOMAD="${DATABASES}/gnomAD_GRCh38_sites/gnomad.exomes.r2.1.1.sites.liftover_grch38_PASS.vcf.gz"

echo "Running VEP conda."
echo "Plugins limited for stat analysis etc."
echo "For extended version see downstream script."
echo "Stats file will be produced."

# VEP
vep \
--fasta $REF \
--offline \
--use_given_ref \
--everything \
--vcf \
--cache \
--dir_cache $CACHE \
--canonical \
--hgvs \
--symbol \
--mane \
--cache_version 104 \
--dir_plugins ${PLUGINS} \
--plugin SpliceConsensus \
--plugin CADD,${CADD_DIR}/whole_genome_SNVs.tsv.gz,${CADD_DIR}/gnomad.genomes.r3.0.indel.tsv.gz \
--custom ${CLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
-i ${INPUT_DIR}/bcftools_gatk_norm_maf01.recode.vcf.gz \
-o ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small.vcf \
--stats_file ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda_small.html \
--fork 24 \
--force_overwrite

echo "END AT $(date)"
