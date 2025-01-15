#!/bin/bash

# All output files used for publication are listed in the two files `sync_list_figure` and `sync_list_table`.
# This study uses git to track code but the data dir is in gitignore due to privacy.
# This study contains many unused outputs such as figures (>1500) and numerous intermediate tables.
# Therefore, we only need to sync the selected content for a neat publication directory.


set -e

# First move some images from data to images
cp ../data/variant_level/*.png ../images/variant_level/

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
    echo "$file_name"

    local src_file="$src_dir/$line"
    local dest_dir="$dest_base_dir/$dir_name"
    local dest_file="$dest_dir/$file_name"

    mkdir -p "$dest_dir"

    if [ -f "$src_file" ]; then
        cp "$src_file" "$dest_file"
        echo "Copied $src_file to $dest_file"
    else
        echo "Source file does not exist: $src_file"
    fi
  done < "$file_list"
}

# Sync image files
sync_files "$src_img_dir" "$dest_img_base_dir" "$img_file_list"

# Sync table files
sync_files "$src_tbl_dir" "$dest_tbl_base_dir" "$tbl_file_list"

