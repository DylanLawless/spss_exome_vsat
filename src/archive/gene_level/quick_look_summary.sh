#!/bin/bash

read -p "Summarise MCL SKAT logs. Enter the version number (e.g., 5 for log_v5): " version
printf "\nDefault results:\n"

printf "number of MCL clusters: "
ls "./log/log_v$version"/log_skat_gene* | wc -l
printf "number of valid tests: "
grep "p_val" "./log/log_v$version"/log_skat_gene* | wc -l
printf "Top hit: "
grep "p_val" "./log/log_v$version"/log_skat_gene* | cut -f4 -d":" | sort -g | head

printf "\nSubset-resample results:\n"
printf "number of MCL clusters: "
ls "./log/log_v$version"/log_skat_subset* | wc -l
printf "number of valid tests: "
grep "p_val" "./log/log_v$version"/log_skat_subset* | wc -l
printf "Top hit: "
grep "p_val" "./log/log_v$version"/log_skat_subset* | cut -f4 -d":" | sort -g | head
