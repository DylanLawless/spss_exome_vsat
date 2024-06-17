# Load required packages
library(dplyr)
library(tidyr)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)  # for ggarrange() and annotate_figure()
library(scico) # devtools::install_github("thomasp85/scico")
# scico_palette_show()
library(stringr)
library(knitr)
library(patchwork)

# load single case report
df_report <- readRDS("../../data/singlecase/df_report.Rds")

df_report_main_text <- readRDS("../../data/singlecase/df_report_main_text.Rds")

# load clinical info
samples <- read.csv("../../data/cohort_summary_curated/SAMPLE_LIST", header = F)
samples$sample <- samples$V1
samples <- samples |> dplyr::select(-V1)

# Pheno ----
# Create new column "cohort_pheno"
samples$cohort_pheno <- samples$sample

# Replace any value that starts with "setpt" with "0" in the "cohort_pheno" column
samples$cohort_pheno[grep("^setpt", samples$sample)] <- "0"

# Replace any value that does not start with "setpt" with "1" in the "cohort_pheno" column
samples$cohort_pheno[!grepl("^setpt", samples$sample)] <- "1"

# clean IDs
samples <- separate(samples, sample, into = c("V1", "V2", "V3", "V4", "V5"))

samples <- samples |>
  mutate(V1 = ifelse(V1 == "raw", NA, V1))
  
samples <- samples |>
  unite(V1, V2, col = "sample.id", sep = "", na.rm = TRUE)

samples <- samples |> filter(cohort_pheno == 1)

# Clinical data
df <- read.csv("../../data/cohort_summary_curated/sepsis_v2.csv")
names(df)
df <- df |> dplyr::select(
  -exome_dataset_1,
  -exome_dataset_1.1,  
  -exome_dataset_2_path,
  -exome_dataset_1_path,
  -sqlpkey,
  -personal.id)

df$sample.id <- gsub("-", "", df$sample.id)

df <- merge(samples, df, by = "sample.id", all.x = TRUE)

missing_samples <- subset(df, is.na(study.site))
missing_sample_ids <- missing_samples$sample.id

# save for calling in other scripts ----
# df_cohort_clin_feat <-  df
saveRDS(df, file = "../../data/cohort_summary_curated/cohort_summary_curated_r_df_singlecase.Rds")

# text summary description Hmisc ----
df <- df |> dplyr::select(
  # -sample.id,
  -V3,
  -V4,
  -V5)

names(df)
names(df_report_main_text)
df <- merge(df_report_main_text, df, by="sample.id", all=T)

df$group <- "singlecase_damaging"
df$group <- ifelse(is.na(df$ACMG_total_score), "singlecase_NA", df$group)

# continuous: hist plots -----
# Convert the data from wide to long
# df_long <- df |>
#   dplyr::select(which(sapply(df, is.numeric))) |>  # Select numeric columns
#   gather(group, key = "variable", value = "value")  # Convert from wide to long format

# df_long <- df |>
  # dplyr::select(group, which(sapply(df, is.numeric))) |>  # Select numeric columns
  # gather(2:23, key = "variable", value = "value")  # Convert from wide to long format

df_long <- df |> 
  dplyr::select(episode.nr:group) |> 
  dplyr::select(c(group, where(is.numeric))) |>  
  pivot_longer(cols = -c(group), names_to = "variable", values_to = "value")  


# Define a helper function to get the next value in a vector
next_in_list <- function(lst, value) {
  ind <- which(lst == value)
  if (ind < length(lst)) {
    return(lst[ind + 1])
  } else {
    return(lst[ind])
  }
}

# Define a function to create a histogram
create_hist <- function(data, variable_name) {
  # Compute the number of bins based on the range of the data
  bins <- diff(range(data$value, na.rm = TRUE))
  
  # Set the bins to the next largest cap value from the vector
  cap_values <- c(5, 10, 20, 30)
  
  bins <- ifelse(bins <= min(cap_values), min(cap_values), bins)
  bins <- ifelse(bins > min(cap_values) & bins <= next_in_list(cap_values, min(cap_values)), 
                 next_in_list(cap_values, min(cap_values)), bins)
  bins <- ifelse(bins > next_in_list(cap_values, min(cap_values)) & bins <= next_in_list(cap_values, next_in_list(cap_values, min(cap_values))), 
                 next_in_list(cap_values, next_in_list(cap_values, min(cap_values))), bins)
  bins <- ifelse(bins > max(cap_values), max(cap_values), bins)
  
  # Create the plot
  data |>
    ggplot(aes(x = value, fill = group)) +
    geom_density(alpha=0.8, stat = "bin", bins = bins )+
    # geom_histogram(alpha=0.5, position = "dodge", color = "black") +
    theme_minimal() +
    guides(fill=FALSE) +
    scale_fill_scico_d(palette = 'berlin', direction = 1) +
    labs(subtitle = variable_name, 
         x = "", 
         y = "")
}

# Create a list of histograms, one for each variable
plot_list <- lapply(unique(df_long$variable), function(var) {
  create_hist(df_long |> filter(variable == var), var)
})

# Calculate the number of rows based on the number of plots and columns
n_plots <- length(plot_list)
ncol <- 5  # Adjust the number of columns as needed
nrow <- ceiling(n_plots / ncol)
padding <- list(NULL) # Create a list of NULL elements to pad the plot list

# Add padding to the plot list to make its length a multiple of the number of columns
plot_list <- c(plot_list, rep(padding, nrow * ncol - n_plots))

# Arrange all plots together
p_combined1 <-  
  annotate_figure(
    ggarrange(plotlist = plot_list, nrow = nrow, ncol = ncol),
    left = textGrob("No. of patients", rot = 90, vjust = 1),
    bottom = textGrob("Value")
  )
p_combined1

# Save combined plot to PDF
ggsave("../../images/cohort_summary_curated/cohort_plots_singlecase_continuous.pdf" ,plot = p_combined1, height = 10, width = 10)

# continuous stat ----
# Define a function to perform Kruskal-Wallis test and return p-value

# perform_kruskal_wallis_test <- function(data, variable_name) {
#   # Perform Kruskal-Wallis test
#   kw_test <- kruskal.test(data$value ~ data$group)
#   
#   # Return p-value
#   return(c(p_value = kw_test$p.value))
# }
# 
# # Apply the function to each variable
# test_results <- lapply(unique(df_long$variable), function(var) {
#   c(variable = var, perform_kruskal_wallis_test(df_long %>% filter(variable == var), var))
# })

# Drop NAs
# continuous stat ----
# Define a function to perform Kruskal-Wallis test and return p-value
perform_kruskal_wallis_test <- function(data, variable_name) {
  # Exclude NA values
  data <- data[!is.na(data$group) & !is.na(data$value), ]
  
  # Perform Kruskal-Wallis test
  kw_test <- kruskal.test(data$value ~ data$group)
  
  # Return p-value
  return(c(p_value = kw_test$p.value))
}

# Apply the function to each variable
test_results <- lapply(unique(df_long$variable), function(var) {
  c(variable = var, perform_kruskal_wallis_test(df_long %>% filter(variable == var), var))
})


# Unlist the results and bind them into a data frame
test_results_df <- do.call(rbind, test_results) |> as.data.frame()
# Reset row names
rownames(test_results_df) <- NULL

test_results_df$p_value <- as.numeric(test_results_df$p_value)

# Plot the -log10 of the p-values for each variable
stat_continuous <- 
  test_results_df |>
  ggplot(aes(x = variable, y=-log10(p_value))) +
  geom_point() +
  theme_minimal() +
  geom_hline(linetype="dotted", 
             yintercept=-log10( .05/nrow(test_results_df) ),
             color="red") + # Bonferroni correction threshold
  theme(axis.text.x  = element_text(angle=45, hjust=1, vjust=1)) # Rotate the x label

stat_continuous

# categorical: bar plots -----
# Convert the data from wide to long
# df_long <- df |>
#   dplyr::select(which(sapply(df, is.character))) |>  # Select character columns
#   gather(key = "variable", value = "value")  # Convert from wide to long format


df_long <- df |> 
  dplyr::select(episode.nr:group) |> 
  dplyr::select(c(group, where(is.character))) |>  
  pivot_longer(cols = -c(group), names_to = "variable", values_to = "value")  

# Filter out date variables
df_long <- df_long |> 
  filter(!(variable %in% c("hosp.adm", "hosp.dis", "picu.adm", "picu.dis", "bc.sampling", "death.date"))) 


df_long <- df_long |>
  mutate(value = str_replace_all(
    value, 
    c("abdominal infection" = "abdom. infec.",
      "haematologic or immunologic" = "haemat. or immun.",
      "third string" = "third replacement",
      "technology dependence" = "tech. depend.",
      "Group " = "Grp. ",
      "Streptococci" = "Strep.",
      "other middle eastern" = "other mid eastern",
      "congenital or genetic" = "congen. or genetic"
      )))

df_long <- df_long |>
  mutate(variable = str_replace_all(
    variable, 
    c("age.category2" = "age.category",
      "ethnicity" = "self.reported.ethnicity"
    )))

# drop pheno, because all cases
df_long <- df_long |> dplyr::filter(!variable == "cohort_pheno")

# # Define a function to create a bar plot
create_bar <- function(data, variable_name) {
  # Create the plot
  data |> 
  ggplot(aes(x = value, fill = group)) +
    geom_bar(color = "black", position="dodge") +
    # geom_density(alpha=0.8, stat = "count")+
    scale_fill_scico_d(palette = 'berlin', direction = 1) +
    guides(fill=FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(subtitle = variable_name,
         x = "",
         y = "")
}

  
# Create a list of bar plots, one for each variable
plot_list <- lapply(unique(df_long$variable), function(var) {
  create_bar(df_long %>% filter(variable == var), var)
})

# Calculate the number of rows based on the number of plots and columns
n_plots <- length(plot_list)
ncol <- 5
nrow <- ceiling(n_plots / ncol)
padding <- list(NULL) # Create a list of NULL elements to pad the plot list

# Add padding to the plot list to make its length a multiple of the number of columns
plot_list <- c(plot_list, rep(padding, nrow * ncol - n_plots))

# Arrange all plots together
p_combined2 <-  
  annotate_figure(
    ggarrange(plotlist = plot_list, nrow = nrow, ncol = ncol),
    left = textGrob("No. of patients", rot = 90, vjust = 1),
    bottom = textGrob("Category")
  )
p_combined2
# Save combined plot to PDF
ggsave("../../images/cohort_summary_curated/cohort_plots_singlecase_categorical.pdf", plot = p_combined2, height = 10, width = 10)

# categorical chi-sqr ----
# Define a function to perform chi-squared test and return p-value, X-squared and degrees of freedom
perform_chi_squared_test <- function(data, variable_name) {
  # Create a contingency table
  contingency_table <- table(data$group, data$value)
  
  # Perform chi-squared test
  chisq_test <- chisq.test(contingency_table)
  
  # Return p-value, X-squared, degrees of freedom
  return(c(p_value = chisq_test$p.value, X_squared = chisq_test$statistic, df = chisq_test$parameter))
}


# Apply the function to each variable
test_results <- lapply(unique(df_long$variable), function(var) {
  c(variable = var, perform_chi_squared_test(df_long %>% filter(variable == var), var))
})

# Unlist the results and bind them into a data frame
test_results_df <- do.call(rbind, test_results) |> as.data.frame()
# Reset row names
rownames(test_results_df) <- NULL
class(test_results_df)
names(test_results_df)

test_results_df$p_value <- as.numeric(test_results_df$p_value)

stat_categorical <- test_results_df |>
  ggplot(aes(x = variable, y=-log10(p_value))) +
  geom_point() +
  theme_minimal() +
  geom_hline(linetype="dotted", 
             yintercept=-log10( .05/nrow(test_results_df) ),
             color="red") + # Bonferroni correction threshold
  theme(axis.text.x  = element_text(angle=45, hjust=1, vjust=1)) # Rotate the x label

stat_continuous
  
# Now, test_results_df contains the variable name, p-value, X-squared value, and degrees of freedom for each variable.


# patchwork ----
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (p_combined2 | p_combined1) + plot_annotation(tag_levels = 'A')

patch2 <- (stat_continuous / stat_categorical) + plot_annotation(tag_levels = 'A')
  

# red is VSAT contributer
ggsave("../../images/cohort_summary_curated/cohort_plots_singlecase_cat_con_bluePATHO.pdf", plot = patch1, height = 10, width = 20)

ggsave("../../images/cohort_summary_curated/cohort_plots_singlecase_cat_con_statistic.pdf", plot = patch2, height = 6, width = 6)

# kable latex tables ----
ACMG_total_score_cutoff_pathogenic <- 6
# df_summaries <- df |> filter(ACMG_score >= 4)
df_summaries <- df |> filter(ACMG_total_score >= ACMG_total_score_cutoff_pathogenic)

df_summaries |> 
  group_by(genotype) |>
  summarise(variants = n()) |>
  kable("latex", booktabs = TRUE)

df_summaries |> 
  group_by(Inheritance) |>
  summarise(n())|>
  kable("latex", booktabs = TRUE)

df_summaries |> 
  ungroup() |>
  group_by(ACMG_total_score, Consequence, IMPACT) |>
  summarise(unique_variants = n()) |>
  arrange(desc(ACMG_total_score), IMPACT, unique_variants) |>
  kable("latex", booktabs = TRUE)


# final table ----
df_summaries <- df_summaries |> dplyr::select(-Strong_patho, -Moder_patho, -Suppor_patho)

saveRDS(df_summaries, file="../../data/singlecase/ACMGuru_singlecase_df_report_cohort_data.Rds")
write.csv(df_summaries,  paste0("../../data/singlecase/ACMGuru_singlecase_df_report_cohort_data.csv"))


df_dedup <- df_summaries |> dplyr::select(sample.id, study.site: psofa.hem) |> unique()
write.csv(df_dedup,  paste0("../../data/singlecase/ACMGuru_singlecase_df_report_dedup.csv"))

# df_report_main_text |>
#   filter(ACMG_total_score >= 6) |> 
#   dplyr::select(SYMBOL, 
#                 Consequence, 
#                 sample.id,
#                 genotype,
#                 Inheritance) |>
#   kable("latex", booktabs = TRUE)

# collapse Inheritance column 
df_report_main_text %>%
  filter(ACMG_total_score >= 6) %>%
  group_by(SYMBOL, Consequence, sample.id, genotype) %>%
  summarise(Inheritance = paste(Inheritance, collapse = ", ")) %>%
  ungroup() %>%
  kable("latex", booktabs = TRUE)



# Textual clinical report ----
# Load necessary libraries
library(dplyr)
library(stringr)

# Assume df_summaries is pre-loaded and appropriately structured
df_summaries <- df_summaries %>%
  mutate(across(where(is.character), as.factor))

df_summaries_length <- length(df_summaries)

# Summarize each column with numeric summaries for numeric data and frequency counts for factors
summarize_data <- function(data) {
  data %>%
    summarise(
      across(where(is.numeric), 
             ~ paste(median(., na.rm = TRUE),
                     "[", min(., na.rm = TRUE), 
                     "-", max(., na.rm = TRUE), "]"),
             .names = "{.col}.stats"),
      across(where(is.factor),
             ~ paste(names(table(.)), table(.), sep=": ", collapse=", "),
             .names = "{.col}.distribution"))
}

# Apply the function to summarize the dataset
cohort_features <- summarize_data(df_summaries)

print(cohort_features)

# Define variable descriptions based on provided definitions
# Define variable descriptions and group them by category
variable_descriptions <- list(
  "Demographic Information" = list(
    episode.nr.stats = "Number of previous sepsis episodes registered in the same child",
    age.days.stats = "Patient age at blood culture sampling in days",
    gender.distribution = "Gender distribution of the cohort",
    age.category2.distribution = "Age categories based on the time of blood culture sampling",
    ethnicity.distribution = "Ethnic background of the cohort"
  ),
  "Hospitalization Data" = list(
    hosp.dur.stats = "Total length of hospital stay in days",
    hosp.dur.post.bc.stats = "Hospital stay length after blood culture sampling in days",
    picu.dur.stats = "Total length of stay in the Pediatric Intensive Care Unit (PICU) in days",
    picu.dur.post.bc.stats = "PICU stay length after blood culture sampling in days",
    hospital.adm.delay.stats = "Delay in hospital admission from time of initial presentation in days"
    # hosp.dis.distribution = "Dates of hospital discharges",
    # picu.dis.distribution = "Dates of PICU discharges"
  ),
  "Clinical Outcomes" = list(
    cons05.score.stats = "Total number of organ failures as defined by the 2005 consensus",
    pelod.score.stats = "Total number of organ failures as defined by the Pediatric Logistic Organ Dysfunction Score (PELOD-2)",
    psofa.score.stats = "Total number of organ failures as defined by the 2017 pSOFA",
    outcome.death.distribution = "Mortality outcomes within 30 days post-admission",
    outcome.picu.los.distribution = "Impact of PICU length of stay on outcomes",
    outcome.death.picu.distribution = "Mortality outcomes specifically within the PICU settings"
  ),
  "Organ Failures and Sepsis Details" = list(
    cons05.cvs.distribution  = "Cardiovascular failure score under 2005 consensus definitions",
    cons05.resp.distribution = "Respiratory failure score under 2005 consensus definitions",
    cons05.cns.distribution = "Central nervous system failure score under 2005 consensus definitions",
    cons05.ren.distribution = "Renal failure score under 2005 consensus definitions",
    cons05.hep.distribution = "Hepatic failure score under 2005 consensus definitions",
    cons05.hem.distribution = "Hematological failure score under 2005 consensus definitions"
  ),
  "Pathogen Information" = list(
    clin.focus.distribution = "Primary clinical focus or reasons for medical intervention",
    pathogen.grp.distribution = "Types of pathogens identified in blood cultures",
    cvc.clabsi.distribution = "Incidence of central venous catheter-associated bloodstream infections"
  )
)

# sample.id.distribution = "Unique identifier for each episode, used to track blood samples",
# ACMG_highest.distribution = "Highest classification of genetic variants per the ACMG standards",
# SYMBOL.distribution = "Gene symbols associated with observed genetic variants",
# rownames.distribution = "Row identifiers based on the structure of the dataset",
# chr.distribution = "Chromosomal locations for the identified genetic variants",
# HGVSp.distribution = "Protein changes caused by genetic variants",
# HGVSc.distribution = "Coding DNA sequence changes due to genetic variants",
# Consequence.distribution = "Functional consequences of genetic variants on protein",
# IMPACT.distribution = "Impact level of genetic variants according to their severity",
# Inheritance.distribution = "Inheritance patterns of identified genetic conditions",
# CLIN_SIG.distribution = "Clinical significance of identified genetic variants",
# bc.sampling.distribution = "Dates when blood cultures were sampled",
# hosp.adm.distribution = "Dates of hospital admissions",
# picu.adm.distribution = "Dates of PICU admissions",
# death.date.distribution = "Dates of death within the cohort",
# hosp.dis.distribution = "Dates of hospital discharges",
# picu.dis.distribution = "Dates of PICU discharges",


# Function to produce the full clinical report with subheadings and missing data check
create_full_report <- function(descriptions, feature_list) {
  report <- ""
  for (category in names(descriptions)) {
    report <- paste(report, sprintf("### %s\n", category), sep="\n")  # Add category heading
    for (variable in names(descriptions[[category]])) {
      description <- descriptions[[category]][[variable]]
      # Check if the variable exists in the dataset
      if (!is.null(feature_list[[variable]])) {
        summary <- feature_list[[variable]]
        report <- paste(report, sprintf("%s (%s): %s.", description, variable, summary), sep=" ")
      } else {
        # If the variable is missing, provide a note in the report
        report <- paste(report, sprintf("%s (%s): \n !!! Data not available: check the variable name in variable_descriptions !!! \n ", description, variable), sep=" ")
      }
    }
    report <- paste(report, "\n", sep="\n")  # Add spacing between groups for readability
  }
  return(report)
}

# Example usage
full_report <- create_full_report(variable_descriptions, cohort_features)
cat(full_report)



# Define the text to prepend with dynamic values included
prepend_text <- paste0(
  "ACMGuru also performed a cohort summary for single case analysis using the study handbook variables. ",
  "Variants which were considered pathogenic or likely pathogenic (ACMG total score = ", 
  ACMG_total_score_cutoff_pathogenic, 
  ") were retained. After analysis, ", 
  df_summaries_length, 
  " patients were identified. Values are indicated as median [min - max]."
)

# Print to check
cat(prepend_text)

# 
# The total score is then compared to thresholds to assign the final verdict:
# 	
# Pathogenic if greater than or equal to 10,
# Likely Pathogenic if between 6 and 9 inclusive,
# Uncertain Significance if between 0 and 5,
# Likely Benign if between -6 and -1,
# Benign if less than or equal to -7.

# Combine the prepend text with the full report
full_report_with_intro <- paste(prepend_text, full_report, sep="")

# Print the full report with intro to the console
cat(full_report_with_intro)

# Define the file path (adjust the path as needed for your directory structure)
text_path <- "../../data/singlecase/ACMGuru_singlecase_df_report_clinical_text.txt"

# Save the report as a text file
writeLines(full_report_with_intro, text_path)
