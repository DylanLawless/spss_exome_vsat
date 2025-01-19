#!/bin/bash
set -e

# how to run
#	1.  pick your version of parameters (v1, v2, v3, v4)
#	2A. run full analysis:                   sh vsat_launch.sh
#	2B. re-run subset with large resampling: sh vsat_subset_resample.sh

# about
#	- designed to take the variables for our multiple analysis options, 
#	- it passes them to a submission template 
#	- the template is used to create the job submission file: sbatch_script.sh.
#	- this sbatch_script.sh is submit and:
#	- saves a copy of the submission information
#	- the correct number of array jobs
#	- produces logs

# input
# - chr split vcf files from cohort case/control
# - MCL_ID list

# output
#	- results:	./log/log_v(version)/
#	- job info:	./log/log_v(version)/
#	- job logs:	./log/err_out_(version)/

# define the variables
CRITICAL_n_samples="940" # number of total samples
# n_variant_sets="947" # number of pathways to test from ProteoMCLusR
# n_variant_sets="1" # run a single pathway test
n_variant_sets="22" # run vsat_run.R once per chromosome

# read from command line argument if available, otherwise prompt
version_suffix="${1:-}"
if [[ -z "${version_suffix}" ]]; then
    read -p "Enter the version number (e.g., 5 for v5): " version
    version_suffix="v${version}"
fi

# pick your version based on: cohort MAC, gnomad MAC, and IMPACT
#	v1 50, 1e-2 ( .01), HIGH
#	v2 10, 1e-3 (.001), HIGH
#	v3 50, 1e-2 ( .01), HIGH-MODERATE
#	v4 10, 1e-3 (.001), HIGH-MODERATE

# example gnomad threshold frequencies
#	1e-1 = <14'000 on gnomad: 14000/140'000 = 0.1
#	1e-2 = <1'400  on gnomad: 1400/140'000 = 0.01
#	1e-3 = <140    on gnomad: 140/140'000 = 0.001
#	1e-4 = <14     on gnomad: 14/140'000 = 0.0001
#	1e-5 = <2      on gnomad: 2/140'000 = 0.00001

# Set variables according to version_suffix
case "${version_suffix}" in
	"v1")
		cohort_max_carriers="50"
		gnomad_rare_thresh="1e-2"
		filter_method="HIGH"
		;;
	"v2")
		cohort_max_carriers="10"
		gnomad_rare_thresh="1e-3"
		filter_method="HIGH"
		;;
	"v3")
		cohort_max_carriers="50"
		gnomad_rare_thresh="1e-2"
		filter_method="MODERATE-HIGH"
		;;
	"v4")
		cohort_max_carriers="10"
		gnomad_rare_thresh="1e-3"
		filter_method="MODERATE-HIGH"
		;;
	*)
		echo "Invalid version_suffix: ${version_suffix}"
		exit 1
		;;
esac

# default resampling value
n_Resampling=100000

# check for subset_resampling
Subset_resampling="${2:-}" # file path during submission for resampling a subset
Resampling_Mode="${3:-}" # resampling mode max/min

if [[ -n "${Subset_resampling}" ]]; then
    echo "Running subset resampling analysis with file: ${Subset_resampling}" # Debugging line
    file_suffix="subset_gene_whole_genome_interpret_skat_${version_suffix}"
    output_param="subset_"
    array_param=$(wc -l < "$Subset_resampling")
    if [[ "${Resampling_Mode}" == "max" ]]; then
        n_Resampling=1000000000
    else
        n_Resampling=1
    fi
    echo "file_suffix is: ${file_suffix}" # Debugging line
    echo "output_param is: ${output_param}" # Debugging line
    echo "array_param is: ${array_param}" # Debugging line
    echo "n_Resampling is: ${n_Resampling}" # Debugging line
else
    echo "Running full analysis" # Debugging line
    file_suffix="gene_whole_genome_interpret_skat_${version_suffix}"
    output_param=""
    array_param=$n_variant_sets
    echo "file_suffix is: ${file_suffix}" # Debugging line
    echo "output_param is: ${output_param}" # Debugging line
    echo "array_param is: ${array_param}" # Debugging line
    echo "n_Resampling is: ${n_Resampling}" # Debugging line
fi

cat > sbatch_script.sh << EOL
#!/bin/bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 50
#SBATCH --mem-per-cpu 9G
#SBATCH --time 06:00:00
#SBATCH --job-name=gene_${version_suffix}
#SBATCH --output=./log/out_err_${version_suffix}/${output_param}skat_gene_%J.out
#SBATCH --error=./log/out_err_${version_suffix}/${output_param}skat_gene_%J.err
#SBATCH --comment="scitas_cost"
#SBATCH --array=1-${array_param}
# #SBATCH --cpus-per-task 10

echo "START AT \$(date)"
set -e
module load gcc/11.3.0
module load r/4.1.3

# for subset_resampling use the correct line 
# from the Subset_resampling file as the current_MCLID
if [[ -n "${Subset_resampling}" ]]; then
	current_MCLID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$Subset_resampling")
else
	current_MCLID=\$(printf \$SLURM_ARRAY_TASK_ID)
fi

echo "current_MCLID is: \$current_MCLID"

# call Rscript with all variables passed as arguments
Rscript "vsat_run.R" \
	\$current_MCLID \
	${version_suffix} \
	${file_suffix} \
	${filter_method} \
	${gnomad_rare_thresh} \
	${cohort_max_carriers} \
	${CRITICAL_n_samples} \
	${n_Resampling} \
	\$SLURM_ARRAY_TASK_ID

echo "SCITAS: jed, 1 task per node 50G, 1 cpu, 4hours , array"
echo "END AT \$(date)"
EOL

# copy the new job script to the log directory
cp sbatch_script.sh ./log/log_${version_suffix}/${output_param}sbatch_script.sh

# submit the job and redirect submission output to your log directory
sbatch sbatch_script.sh > ./log/log_${version_suffix}/${output_param}sbatch_output.txt 2>&1

# clean up
rm sbatch_script.sh
