#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
# #SBATCH --ntasks 24
# #SBATCH --cpus-per-task 1
#SBATCH --mem 30G
#SBATCH --time 00:40:00
#SBATCH --job-name=vep
#SBATCH --output=./log/vep/vep_%J.out
#SBATCH --error=./log/vep/vep_%J.err
#SBATCH --comment="scitas_cost"

# CONDA ENV
# activate before submission
# conda activate vep

echo "START AT $(date)"
set -e

# Master variables
source variables.txt

# Path
INPUT_DIR="${WORK_DIR}/pre_annotation_output"
OUTPUT_DIR="${WORK_DIR}/annotation"

# db
REF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"
GTF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/gtf_and_gff/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz"
GFF="/work/gr-fe/saadat/Reference_Genome/GRCH38_no_alt/gtf_and_gff/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation_for_vep.gff.gz"
CACHE="/work/gr-fe/databases/vep_hg38/cache_hg38_ensembl"
PLUGINS="/work/gr-fe/databases/vep_hg38/plugins"
LOFTEE="/work/gr-fe/databases/vep_hg38/loftee"
G2P_PANNEL="/work/gr-fe/saadat/panels/Primary_immunodeficiency.csv"
CADD_DIR="/work/gr-fe/databases/CADD_v1.6_GRCH38_hg38"

# /////////////////////////////////////////////////////////////////////////////
# Load all modules required
# /////////////////////////////////////////////////////////////////////////////
module load intel # jed

# /////////////////////////////////////////////////////////////////////////////
# VEP
# /////////////////////////////////////////////////////////////////////////////
VEP_dir=/work/gr-fe/lawless/tool/ensembl-vep

# $VEP_dir/vep \

# Note with conda, I think if i use via tmux, only the initial tab has the 
# conda env active.

vep \
--fasta $REF \
--use_given_ref \
--everything --offline \
--vcf \
--cache \
--dir_cache $CACHE \
--canonical \
--hgvs \
--symbol \
--mane \
--cache_version 104 \
-i ${INPUT_DIR}/bcftools_gatk_norm.vcf.gz \
-o ${OUTPUT_DIR}/bcftools_gatk_norm_vep_conda.vcf \
--fork 24 \
--force_overwrite
# --flag_pick \
# --stats_file ${OUTPUT_DIR}/vep_cache_ensembl_flagpick.html \

# # tab output
# $VEP_dir/vep \
# --tab \
# --fasta $REF \
# --use_given_ref \
# --everything --offline \
# --cache \
# --dir_cache $CACHE \
# --canonical \
# --hgvs \
# --symbol \
# --mane \
# --flag_pick \
# --cache_version 104 \
# -i ${INPUT_DIR}/bcftools_gatk_norm.vcf.gz \
# -o ${OUTPUT_DIR}/vep_cache_ensembl_flagpick.tsv \
# --stats_file ${OUTPUT_DIR}/vep_cache_ensembl_flagpick.html \
# --fork 24 \
# --force_overwrite

# /////////////////////////////////////////////////////////////////////////////
# tests
# /////////////////////////////////////////////////////////////////////////////
#--gtf $GTF \
# vep \ # conda version

# # Activate vep from conda
# eval "$(conda shell.bash hook)"
# conda activate ensembl-vep


# Load all modules required
# Load all modules required
#--dir_plugins $PLUGINS \
#–transcript_version \
#--plugin CADD,${CADD_DIR}/whole_genome_SNVs.tsv.gz,${CADD_DOR}/gnomad.genomes.r3.0.indel.tsv.gz \
#vep \
#--offline --vcf --everything \
#--fasta $REF \
#--dir_cache $CACHE \
#--dir_plugins $PLUGINS \
#--plugin Condel,$PLUGINS/config/Condel/config/ \
#--plugin SpliceConsensus \
#--plugin Downstream \
#--terms SO \
#--af_gnomad \
#–coding_only \
#–transcript_version \
#--cache \
#--cache_version 104 \
#-i $INPUT \
#-o $OUTPUT \
#--fork 24 \
#--force_overwrite
#--assembly GRCh38 \

# module load samtools/1.10
# module load htslib/1.10.2

# I think the following may be incorrect - vep works on fidis but have perl trouble on helvetios
# While installing either vcftools, bcftools, or vt, a second perl install occurred
# making it difficult to see why I get a failure using vep due to perl.
# Instead use conda ensemble-vep:
#$ conda install -c bioconda ensembl-vep
#$ conda install ensembl-vep=105.0-0
echo "END AT $(date)"
