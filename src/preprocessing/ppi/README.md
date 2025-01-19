## Disease model versions
* To figure out the data format we test out different criteria;
* very-rare/rare, LoF/VUS.
* We called these "versions".
* Update this README on completion with final model definitions.

## Input
* Start with the output of GATK pipeline as found in `./data/joint`

## Script order
* The job is run by using: `vsat_launch.sh`

* `vsat_launch.sh` 
	* This is the only script that the user _should_ need to run. 
	* When run is asks for the disease model version (v1-4). 
	* This script takes the user variables and creates an SBATCH submission script for SLURM. 
	* It makes an array of jobs that run the `vsat_run.R` script for analysis.
	* It saves outputs and logs. 
	* If you get any results that require higher resampling, you will instead run `vsat_subset_resample.sh` which will handle the resubmission of the subset as well as re-submitting `vsat_launch.sh`. 
* `vsat_subset_resample.sh` 
	* When we have highly significant p-values, SKAT requires resampling at higher than default (e.g. 1,000,000 > 9,000,000).
	* To save costs and compute time, we initially run all pathways using the default `vsat_launch.sh`.
	* You can see a notice in results if any hits still require further resampling; we can run this script to re-run only the subset to save time/costs.
	* Like the default script, it will ask for your disease model version and it will ask for a p-value to identify pathways that need which further resampling (e.g. v3 with p-value  e-7).
* `vsat_run.R`	- This is main R code doing to work: importing VCF for analysis with SKAT based on MCL pathway IDs and saving output.
* after SkatterProt:
* `post_ppi/23_filter_MCL_genes.sh`: get pathway gene-set VCF for all in file `genelist` (MCL_ID list from output - manual)
* `post_ppi/vsat_log_result_vcurrent.R`
* `AMCGuru_post_ppi/ACMGuru_post_ppi_vcurrent.R`

* single case studies:
* `AMCGuru_singlecase/AMCGuru_singlecase_vcurrent.R`

* gene burden studies:

## Other notes
* Automation:
	* Ideally the last unfinished automation would be to replace `vsat_subset_resample.sh` such that we only need to ever run `vsat_launch.sh`. 
	* The problem is that `vsat_launch.sh` runs ~1000 job array of indeterminant time. 
	* On our dataset it requires ~2 hours to complete 1000 jobs for 1000 samples but the HPC scheduler may decide to delay job launches.
	* Therefore, we cannot know when to run the check if resampling is required (after all jobs complete). Note that you could set a flag for slurm to run afeter the last array finishes.
	* A possible easy fix is to allow SKAT to continue to run until it determines that resampling is sufficient.
	* However, this requires all jobs to have freer memory and times, thereby demoting our jobs' queue optimisation.
	* Each individual jobs can't decide after completing since it only knows its own resuts.
* Improve logic
	* Currently we take input data which is split per chromosome - a logical output for a large dataset from GATK pipeline
	* However, this complicates pathway analysis since: for every pathway we must import the chromosome, check for pathway genes, extract the variant set, store, and repeat until we have all chromosomes complete. Then we run skat on the stored dataset.
	* The logic would be easier to follow if we first pre-processed the input data such that it is already binned into protein pathways. Then we simply run skat on each pathway dataset.
* Mistake in output log: 
	* "Number of phenotypes in the filtered dataset does not match CRITICAL_n_samples."
	* This happened because we count to confirm 2 phenotypes for case/cotronl 0/1, but compared to sample count instead of "2".
	* Fixed with new warning: "Number of phenotypes in the filtered dataset does not eqaul 2 (case/control)."
