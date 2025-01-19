#!/bin/bash

# All output files used for publication are listed in the two files `sync_list_figure` and `sync_list_table`.
# This study uses git to track code but the data dir is in gitignore due to privacy.
# This study contains many unused outputs such as figures (>1500) and numerous intermediate tables.
# Therefore, we only need to sync the selected content for a neat publication directory.


set -e

echo "First run the studyID indexing for publication tables"
Rscript tables_ID_Mapping.R

echo "Then move any images from data to images - some modules do not ouput to ./images by default"
cp ../data/variant_level/*.png ../images/variant_level/
cp ~/web/ProteoMCLustR/data/ppi_user_ouput/ppi_line_dist_whole_genome_v1_compared.pdf ../images/ProteoMCLustR
cp  ../data/archipelago/*pdf  ../images/archipelago/
cp  ../data/archipelago/*png  ../images/archipelago/
cp  ../data/archipelago/*jpg  ../images/archipelago/


# Image directories
src_img_dir="../images"
dest_img_base_dir="/Users/dylanlawless/Library/Mobile Documents/com~apple~CloudDocs/post/manuscripts_epfl/spss_exomes/lawless2023_spss_exome_vcurrent/publication_figure_sync"

# Table directories
src_tbl_dir="../data"
dest_tbl_base_dir="/Users/dylanlawless/Library/Mobile Documents/com~apple~CloudDocs/post/manuscripts_epfl/spss_exomes/lawless2023_spss_exome_vcurrent/publication_table_sync"

# Image files to sync
img_file_list="sync_list_figure"
# Table files to sync (create a similar list file for tables if needed)
tbl_file_list="sync_list_table"  # Specify the list file for table files

# Function to sync files from source to destination
sync_files () {
  local src_dir=$1
  local dest_base_dir=$2
  local file_list=$3

  while IFS= read -r line; do
    [[ "$line" =~ ^#.*$ ]] && continue

    dir_name=$(dirname "$line")
    file_name=$(basename "$line")
    # echo "$file_name"

    local src_file="$src_dir/$line"
    local dest_dir="$dest_base_dir/$dir_name"
    local dest_file="$dest_dir/$file_name"

    mkdir -p "$dest_dir"

    if [ -f "$src_file" ]; then
        cp "$src_file" "$dest_file"
       # echo "Copied $src_file to $dest_file"
    else
        echo "Source file does not exist: $src_file"
    fi
  done < "$file_list"
}

# Sync image files
sync_files "$src_img_dir" "$dest_img_base_dir" "$img_file_list"

# Sync table files
sync_files "$src_tbl_dir" "$dest_tbl_base_dir" "$tbl_file_list"

echo "All files synced from sync_list_figure and sync_list_table."


echo "Transfer standalone reference tables"
cp ../ref/annotation_datasets/data/datasources_compiled.tsv "${dest_tbl_base_dir}/ref/Table_S1_datasources_compiled.tsv"
cp ../ref/10875_2022_1289_MOESM2_ESM_DLcleaned.tsv          "${dest_tbl_base_dir}/ref/Table_S2_10875_2022_1289_MOESM2_ESM_DLcleaned.tsv"
cp ../ref/acmg_criteria_table.txt                           "${dest_tbl_base_dir}/ref/Table_S3_acmg_criteria_table.txt"
cp ../ref/acmg_criteria_table_caveats.txt                   "${dest_tbl_base_dir}/ref/Table_S4_acmg_criteria_table_caveats.txt"
cp ../ref/acmg_scoring_system_points.txt                    "${dest_tbl_base_dir}/ref/Table_S5_acmg_scoring_system_points.txt"
cp ../ref/acmg_scoring_system_for_tally.txt                 "${dest_tbl_base_dir}/ref/Table_S6_acmg_scoring_system_for_tally.txt"
cp ../ref/varsome_calibrated_insilico_thresholds.tsv        "${dest_tbl_base_dir}/ref/Table_S7_varsome_calibrated_insilico_thresholds.tsv"
cp "${dest_tbl_base_dir}/ACMGuru_post_ppi/df_report_discussion_22_586_836.csv" \
	"${dest_tbl_base_dir}/ACMGuru_post_ppi/Table_S8_df_report_discussion_22_586_836.csv"
cp "${dest_tbl_base_dir}/ACMGuru_post_ppi/ACMGuru_post_ppi_genetic_df_report_main_text_22_586_836_indexedstudyIDs.csv" \
	"${dest_tbl_base_dir}/ACMGuru_post_ppi/Table_S9_ACMGuru_post_ppi_genetic_df_report_main_text_22_586_836_indexedstudyIDs.csv"







echo "END"
