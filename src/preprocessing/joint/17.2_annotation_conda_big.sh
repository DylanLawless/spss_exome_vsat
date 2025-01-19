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
G2P_PANNEL="/work/gr-fe/saadat/panels/Primary_immunodeficiency.csv"
dbNSFP="/work/gr-fe/lawless/database/dbNSFP/data/dbNSFP4.2a_grch38.gz"
IntAct="/work/gr-fe/lawless/database/intact"
# REVEL="/work/gr-fe/lawless/ref/revel/new_tabbed_revel_grch38.tsv.gz"
# GNOMAD="/scratch/lawless/gnomad/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
# GNOMAD="${DATABASES}/gnomAD_GRCh38_sites/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
# GNOMAD="${DATABASES}/gnomAD_GRCh38_sites/gnomad.exomes.r2.1.1.sites.liftover_grch38_PASS.vcf.gz"

echo "Running VEP conda."
echo "Plugins include clinvar, inact, cadd, dbnsfp, and output mane, cannonical, etc."
echo "dbnsfp custom selection."
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
--plugin dbNSFP,${dbNSFP},gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_AC,gnomAD_genomes_AF,gnomAD_genomes_AN,gnomAD_genomes_flag,chr,pos\(1-based\),ref,alt,aaref,aaalt,rs_dbSNP,hg19_chr,hg19_pos\(1-based\),hg18_chr,hg18_pos\(1-based\),aapos,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc,Uniprot_entry,HGVSc_ANNOVAR,HGVSp_ANNOVAR,HGVSc_snpEff,HGVSp_snpEff,HGVSc_VEP,HGVSp_VEP,APPRIS,GENCODE_basic,TSL,VEP_canonical,cds_strand,refcodon,codonpos,codon_degeneracy,Ancestral_allele,AltaiNeandertal,Denisova,VindijiaNeandertal,SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,Aloft_Fraction_transcripts_affected,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,CADD_raw,CADD_raw_rankscore,CADD_phred,CADD_raw_hg19,CADD_raw_rankscore_hg19,CADD_phred_hg19,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,fathmm-XF_coding_score,fathmm-XF_coding_rankscore,fathmm-XF_coding_pred,Eigen-raw_coding,Eigen-raw_coding_rankscore,Eigen-phred_coding,Eigen-PC-raw_coding,Eigen-PC-raw_coding_rankscore,Eigen-PC-phred_coding,GenoCanyon_score,GenoCanyon_rankscore,integrated_fitCons_score,integrated_fitCons_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_rankscore,HUVEC_confidence_value,LINSIGHT,LINSIGHT_rankscore,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian,phyloP30way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons30way_mammalian,phastCons30way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,bStatistic,bStatistic_converted_rankscore,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue,Geuvadis_eQTL_target_gene \
--custom ${CLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
--plugin IntAct,mutation_file=${IntAct}/mutations.tsv,mapping_file=${IntAct}/mutation_gc_map.txt.gz \
--plugin IntAct,mutation_file=${IntAct}/mutations.tsv,mapping_file=${IntAct}/mutation_gc_map.txt.gz,minimal=1 \
-i ${INPUT_DIR}/bcftools_gatk_norm_maf01.recode.vcf.gz \
-o ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda.vcf \
--stats_file ${OUTPUT_DIR}/bcftools_gatk_norm_maf01.recode_vep_conda.html \
--fork 24 \
--force_overwrite

# --plugin dbNSFP,${dbNSFP},ALL \
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
