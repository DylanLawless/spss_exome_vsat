#!/bin/bash
set -e

# processes result from vsat_lauch.sh / vsat_run.R

read -p "Enter the version number (e.g., 5 for log_v5): " version

# define the variable
version_suffix="v$version"

# extract results file
cd ./log/log_${version_suffix}/
printf "\ngetting gene log result"
echo "log_file:Description:MCL_ID:p_val:gene_count:var_count:AC:AC:AC:AC:AC:AC:AC:AC:AC:AC" > vsat_log_${version_suffix}_result.txt
grep "Results for" log_skat_gene_whole_genome_interpret_skat_${version_suffix}_chr_*_SYMBOL_*.log >> vsat_log_${version_suffix}_result.txt
# printf "\n..................... done [1 of 3]\n"
printf "\n..................... done [1 of 1]\n"

# # extract genes file
# printf "\ngetting vsat log genes"
# grep "Found genes" log_skat_ppi_whole_genome_interpret_skat_${version_suffix}_*log | cut -f10-11 -d "_" > vsat_log_${version_suffix}_result_genes.txt 
# printf "\n.......................... done [2 of 3]\n"

# # extract variants and MAC per VSAT (MCL_ID) 
# printf "\ngetting vsat log variants"
# echo "MCL_ID test.snp.mac SNP" > vsat_log_${version_suffix}_result_snp_mac.txt
# for file in log_skat_ppi_whole_genome_interpret_skat_${version_suffix}_*log; do
# 	awk '/MCL_ID test.snp.mac/{flag=1;next} flag' "$file" \
# 		| tr -s ' ' \
# 		| cut -f2-4 -d" "
# done >> vsat_log_${version_suffix}_result_snp_mac.txt
# printf "\n............................... done [3 of 3]\n"
