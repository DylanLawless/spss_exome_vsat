#!/bin/bash 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 30G
#SBATCH --time 03:00:00
#SBATCH --job-name=vep
#SBATCH --output=./log/vep/vep_%J.out
#SBATCH --error=./log/vep/vep_%J.err
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
REVEL="${DATABASES}/revel/new_tabbed_revel_grch38.tsv.gz"
GNOMAD="/scratch/lawless/gnomad/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
G2P_PANNEL="/work/gr-fe/saadat/panels/Primary_immunodeficiency.csv"
dbNSFP="/work/gr-fe/lawless/database/dbNSFP/data/dbNSFP4.2a_grch38.gz"
IntAct="/work/gr-fe/lawless/database/intact"
# GNOMAD="${DATABASES}/gnomAD_GRCh38_sites/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
# GNOMAD="${DATABASES}/gnomAD_GRCh38_sites/gnomad.exomes.r2.1.1.sites.liftover_grch38_PASS.vcf.gz"

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
--plugin dbNSFP,${dbNSFP},ALL \
--custom ${CLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
--plugin IntAct,mutation_file=${IntAct}/mutations.tsv,mapping_file=${IntAct}/mutation_gc_map.txt.gz \
--plugin IntAct,mutation_file=${IntAct}/mutations.tsv,mapping_file=${IntAct}/mutation_gc_map.txt.gz,minimal=1 \
-i ${INPUT_DIR}/bcftools_gatk_norm.vcf.gz \
-o ${OUTPUT_DIR}/bcftools_gatk_norm_vep_conda_plug.vcf \
--fork 24 \
--force_overwrite

# --plugin Downstream \
# --plugin Condel,${PLUGINS}/config/Condel/config/ \ 
# WARNING: Plugin 'Condel' went wrong:
 ### Phenotypes plugin: This will take some time but it will only run once per species, assembly and release
 # --plugin Phenotypes,dir=${OUTPUT_DIR},include_types=Gene \
# --plugin LoF,loftee_path:${LOFTEE} \
# --plugin gnomADc,${GNOMAD} \
# --plugin REVEL,${REVEL} \ # requires --assembly
# --plugin Condel,$VEPdb/VEP_plugins/config/Condel/config/ \
# --plugin LoF,loftee_path:${LOFTEE} --dir_plugins /path/to/loftee

# # vep custom annotation
# ./vep [...] --custom Filename , Short_name , File_type , Annotation_type , Force_report_coordinates , VCF_fields

# # For multiple custom files, use:
# ./vep [...] --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
#             --custom TOPMED_GRCh38_20180418.vcf.gz,topmed_20180418,vcf,exact,0,TOPMED \
#             --custom UK10K_COHORT.20160215.sites.GRCh38.vcf.gz,uk10k,vcf,exact,0,AF_ALSPAC

## Where the selected ClinVar INFO fields (from the ClinVar VCF file) are:
# - CLNSIG:     Clinical significance for this single variant
# - CLNREVSTAT: ClinVar review status for the Variation ID
# - CLNDN:      ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
# Of course you can select the INFO fields you want in the ClinVar VCF file

# --flag_pick \
# --stats_file ${OUTPUT_DIR}/vep_cache_ensembl_flagpick.html \

echo "END AT $(date)"
