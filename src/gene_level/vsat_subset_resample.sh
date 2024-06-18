#!/bin/bash
printf "You have launched this script because"
printf "vsat_launch has complete but some vsat\ndid not have sufficient resampling.\n"
printf "Re-run subset with large resampling (20,000,000).\n\n"
read -p "Enter the version number (e.g., 5 for v5): " version
read -p "What is your p-value for resampling? (e.g., e-07): " crit

grep "p_val" "./log/log_v$version"/log* | grep "$crit" | cut -f3,4,5,6 -d":" \
	> subset_resampling_mclid_temp

printf "Your data for resampling (first 20):\n"
cat subset_resampling_mclid_temp | head -20

# note the logs currently start with a empty space
cat subset_resampling_mclid_temp \
	| tr -s ' ' \
	| cut -f2 -d" " \
	> subset_resampling_mclid.txt

printf "\nMCL_IDs (first 20):\n"
cat subset_resampling_mclid.txt | head -20
rm subset_resampling_mclid_temp

# Ask for the Resampling Mode
read -p "Enter the resampling mode (min [1]) or max [1,000,000,000]: " Resampling_Mode
while true; do
	read -p "Enter the resampling mode (min/max): " Resampling_Mode
	if [[ "$Resampling_Mode" == "max" || "$Resampling_Mode" == "min" ]]; then
		break
	else
		echo "Invalid input. Please enter 'max' or 'min'."
	fi
done

read -p "Continue? (y/n): " choice
if [ "$choice" == "y" ] || [ "$choice" == "Y" ]; then
	version_suffix="v${version}"
	printf "\nNow submitting:\nvsat_launch.sh subset_resampling_mclid.txt\nwith version ${version_suffix}\n"¬
	sh vsat_launch.sh "${version_suffix}" subset_resampling_mclid.txt "${Resampling_Mode}"
else
	printf "\nScript aborted.\n"
fi
