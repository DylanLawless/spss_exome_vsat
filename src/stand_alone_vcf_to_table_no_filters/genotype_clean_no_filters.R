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
df$genotype[df$genotype_call == "./."] <- "0"
df$genotype[df$genotype_call == "."] <- "0"

genotype_unique <- unique(df$genotype)

if (all(genotype_unique %in% c("0", "1", "2"))) {
	# cat("\nGenotypes found: ", paste(genotype_unique), "\n")
} else {
	cat("\nGenotypes found: ", paste(genotype_unique))
	cat("\nError: Invalid genotype found.\nGenotypes must be format 0,1,2. See genotype_clean.R for details.\n")
	stop("\nStopping analysis.\n")
}

# cat("Removing all genotype: 0.\n")
df$genotype <- as.numeric(df$genotype)
# df <- df |> filter(genotype > 0)

