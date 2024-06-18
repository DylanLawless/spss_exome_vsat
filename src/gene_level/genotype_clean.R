# Used by vsat_run.R
# converts the multiple genotype styles from vcf into 0,1,2

# Create new column "genotype"
df$genotype <- df$genotype_call
df$genotype[df$genotype_call == "0/0"] <- "0"
df$genotype[df$genotype_call == "./0"] <- "0"
df$genotype[df$genotype_call == "0/."] <- "0"
df$genotype[df$genotype_call == "0|0"] <- "0"
df$genotype[df$genotype_call == "0/1"] <- "1"
df$genotype[df$genotype_call == "1/0"] <- "1"
df$genotype[df$genotype_call == "0|1"] <- "1"
df$genotype[df$genotype_call == "./1"] <- "1"
df$genotype[df$genotype_call == "1/."] <- "1"
df$genotype[df$genotype_call == ".|1"] <- "1"
df$genotype[df$genotype_call == "1/1"] <- "2"
df$genotype[df$genotype_call == "1|1"] <- "2"
df$genotype[df$genotype_call == "./."] <- "0"
df$genotype[df$genotype_call == "."] <- "0"

genotype_unique <- unique(df$genotype)

# print message and stop if anything expect 0,1,2
# if your data has another genotype style update the substitutions to include it
if (all(genotype_unique %in% c("0", "1", "2"))) {
} else {
	cat("\nGenotypes found: ", paste(genotype_unique))
	cat("\nError: Invalid genotype found.\nGenotypes must be format 0,1,2. See genotype_clean.R for details.\n")
	stop("\nStopping analysis.\n")
}

df$genotype <- as.numeric(df$genotype)

