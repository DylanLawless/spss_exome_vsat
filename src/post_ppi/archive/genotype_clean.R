# df$genotype_call |> head()
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
df$genotype[df$genotype_call == "./."] <- "NA"
df$genotype[df$genotype_call == "."] <- "NA"

genotype_unique <- unique(df$genotype)

if (all(genotype_unique %in% c("0", "1", "2", "NA"))) {
	# cat("\nGenotypes found: ", paste(genotype_unique), "\n")
} else {
	cat("\nGenotypes found: ", paste(genotype_unique))
	cat("\nError: Invalid genotype found.\nGenotypes must be format 0,1,2,NA. See genotype_clean.R for details.\n")
	stop("\nStopping analysis.\n")
}

# cat("Removing all genotype: 0.\n")
df$genotype <- as.numeric(df$genotype)

# Genotype missingness check
# df_AN <- df %>%
# 	group_by(SNP, genotype) %>%
# 	summarise(n_geno = n()) %>%
# 	mutate(NA_freq = sum(is.na(genotype))
# 								/(sum(n_geno)))



