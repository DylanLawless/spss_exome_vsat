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
library(patchwork)

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
df <- df |>dplyr::select(
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
saveRDS(df, file = "../../data/cohort_summary_curated/cohort_summary_curated_r_df.Rds")

# text summary description Hmisc ----
df <- df |>dplyr::select(
  -sample.id,
  -V3,
  -V4,
  -V5)

# Use the describe function
df_desc <- describe(df)

# Print the result
df_desc

# Generate a LaTeX file from the describe output
# latex(df_desc, file = "../latex/df_desc.tex")


# continuous: hist plots -----
# Convert the data from wide to long
df_long <- df |>
 dplyr::select(which(sapply(df, is.numeric))) |>  #dplyr::select numeric columns
  gather(key = "variable", value = "value")  # Convert from wide to long format

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
    ggplot(aes(x = value, fill = ..x..)) +
    geom_histogram(bins = bins, color = "black") +
    theme_minimal() +
    guides(fill=FALSE) +
    scale_fill_scico(palette = 'nuuk', direction = 1) +
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
ggsave("../../images/cohort_summary_curated/cohort_plots_continuous.pdf" ,plot = p_combined1, height = 10, width = 10)

# categorical: bar plots -----
# Convert the data from wide to long
df_long <- df |>
 dplyr::select(which(sapply(df, is.character))) |>  #dplyr::select character columns
  gather(key = "variable", value = "value")  # Convert from wide to long format

# Filter out date variables
df_long <- df_long |> 
  filter(!(variable %in% c("hosp.adm", "hosp.dis", "picu.adm", "picu.dis", "bc.sampling", "death.date"))) 

library(stringr)
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

# Define a function to create a bar plot
create_bar <- function(data, variable_name) {
  # Create the plot
  data |>
  ggplot(aes(x = value, fill = value)) +
    geom_bar(color = "black") +
    scale_fill_scico_d(palette = 'nuuk', direction = 1) +
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
ggsave("../../images/cohort_summary_curated/cohort_plots_categorical.pdf", plot = p_combined2, height = 10, width = 10)

# patchwork ----
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (p_combined2 | p_combined1) + plot_annotation(tag_levels = 'A')

ggsave("../../images/cohort_summary_curated/cohort_plots_cat_con.pdf", plot = patch1, height = 10, width = 20)

# Summary stats ----
# table categorical ----
# Function to calculate frequency count and percentage for a categorical variable
calculate_category_counts <- function(data, variable_name) {
  # Count the frequency of each category
  category_counts <- data %>%
   dplyr::select(all_of(variable_name)) %>%
    table() %>%
    as.data.frame()
  
  # Rename columns for clarity and add a variable name column
  names(category_counts) <- c("Category", "Count")
  category_counts$Variable <- variable_name
  
  # Calculate the percentage of each category
  total_count <- sum(category_counts$Count)
  category_counts$Percentage <- (category_counts$Count / total_count) * 100
  
  # Round the percentage values to two decimal places
  category_counts$Percentage <- round(category_counts$Percentage, digits = 2)
  
  return(category_counts)
}

# Identify character variables and exclude date variables
categorical_vars <- names(df)[sapply(df, is.character)]
date_vars <- c("hosp.adm", "hosp.dis", "picu.adm", "picu.dis", "bc.sampling", "death.date")
categorical_vars <- setdiff(categorical_vars, date_vars)

# Apply the function to all categorical variables and combine the results
all_category_counts <- lapply(categorical_vars, function(var) {
  calculate_category_counts(df, var)
}) %>% bind_rows()

all_category_counts <- 
  all_category_counts %>% 
  dplyr::select(Variable, Category, Count, Percentage) %>%
  filter(Category !="no")


# Print the combined frequency table
print(all_category_counts)

# Save the combined frequency table to a CSV file
write.csv(all_category_counts, file = "../../data/cohort_summary_curated/all_category_counts.csv", row.names = FALSE)

# table continuous ----
# Convert the data from wide to long (similar to the histogram plot preparation)
df_long <- df |>
 dplyr::select(which(sapply(df, is.numeric))) |>  #dplyr::select numeric columns
  gather(key = "variable", value = "value")  # Convert from wide to long format

# Function to calculate summary statistics for each variable
calculate_continuous_summary_stats <- function(data, variable_name) {
  variable_data <- data %>% filter(variable == variable_name)
  
  summary_stats <- variable_data %>%
    summarise(
      Mean = mean(value, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      Min = min(value, na.rm = TRUE),
      Max = max(value, na.rm = TRUE),
      N = n()
    )
  
  # Round the statistics to two decimal places
  summary_stats$Mean <- round(summary_stats$Mean, digits = 2)
  summary_stats$Median <- round(summary_stats$Median, digits = 2)
  summary_stats$SD <- round(summary_stats$SD, digits = 2)
  summary_stats$Min <- round(summary_stats$Min, digits = 2)
  summary_stats$Max <- round(summary_stats$Max, digits = 2)
  
  
  # Add a variable name column
  summary_stats$Variable <- variable_name
  
  return(summary_stats)
}

# Apply the function to each unique variable in the df_long dataframe
all_continuous_summary_stats <- lapply(unique(df_long$variable), function(var) {
  calculate_continuous_summary_stats(df_long, var)
}) %>% bind_rows()

# Rearrange columns for better readability
all_continuous_summary_stats <- all_continuous_summary_stats %>%
  dplyr::select(Variable, Mean, Median, SD, Min, Max, N)

# Print the combined summary statistics table
print(all_continuous_summary_stats)

# Save the combined summary statistics table to a CSV file
write.csv(all_continuous_summary_stats, file = "../../data/cohort_summary_curated/all_continuous_stats.csv", row.names = FALSE)


# publication table ----
print(all_category_counts)
print(all_continuous_summary_stats)

## cohort ----
select_cohort <- all_category_counts |> 
  filter(Variable == "cohort_pheno") |>
  mutate(Variable = case_when(
  Variable == "cohort_pheno" ~ "Cohort",
  TRUE ~ Variable)) |>
  mutate(Category = case_when(
    Category == "1" ~ "Sepsis patient",
    TRUE ~ Category))

select_cohort

## sex ----
select_sex <- all_category_counts |> 
  filter(Variable == "gender") |>
  mutate(Variable = case_when(
    Variable == "gender" ~ "Sex",
    TRUE ~ Variable))

select_sex

## age ----
select_age <- all_category_counts |> filter(Variable == "age.category2")

select_age <- select_age |>
  mutate(Category = case_when(
  Category == "Preterm newborn" ~ "1.Preterm newborn",
  Category == "Term newborn" ~ "2.Term newborn",
  Category == "Infant (1mt - 1y)" ~ "3.Infant (1mt - 1y)",
  Category == "Toddler (2y - 5y)" ~ "4.Toddler (2y - 5y)",
  Category == "School age (6y - 12y)" ~ "5.School age (6y - 12y)",
  Category == "Adolescent (13y - 17y)" ~ "6.Adolescent (13y - 17y)",
  TRUE ~ Category  # Default 
))  %>%
  arrange(Category)  |>
  mutate(Category = sub("^\\d+\\.", "", Category)) |>
  mutate(Variable = case_when(
    Variable == "age.category2" ~ "Agr group",
    TRUE ~ Variable  # Keeps any other variable unchanged
  ))

select_age

## hospital.acp ----
select_hospaq <- all_category_counts |> 
  filter(Variable == "hospital.acquired") |>
  mutate(Variable = case_when(
    Variable == "hospital.acquired" ~ "Hospital acquired",
    TRUE ~ Variable))

select_hospaq

## icu ----
select_icu <- all_category_counts |> 
  filter(Variable == "icu") |>
  mutate(Variable = case_when(
    Variable == "icu" ~ "ICU",
    TRUE ~ Variable))

select_icu

## comorbid ----
select_comorbidity <- all_category_counts |> 
  filter(Variable == "comorbidity") |>
  mutate(Variable = case_when(
    Variable == "comorbidity" ~ "Comorbidities",
    TRUE ~ Variable))

select_comorbidity

## ccc.final ----
select_comorbidgrp <- all_category_counts |> 
  filter(Variable == "ccc.final") |>
  mutate(Variable = case_when(
    Variable == "ccc.final" ~ "Comorbidity group",
    TRUE ~ Variable))

select_comorbidgrp

## selected_dur_hosp ----
selected_dur_hosp <- 
  all_continuous_summary_stats |> 
  select(Variable, Median, Min, Max, SD) |>
  filter(
    Variable == "hosp.dur.post.bc" # Length of stay after sepsis onset (days)
  )

selected_dur_hosp

## selected_dur_picu ----
selected_dur_picu <- 
  all_continuous_summary_stats |> 
  select(Variable, Median, Min, Max, SD) |>
  filter(
    Variable == "picu.dur.post.bc" # Length of ICU stay after sepsis onset
    # Variable == "picu.dur" |  # Length of ICU stay
  )


# rbind ----
print("Demographics and clinical characteristics")

selec_cat_bind <- 
  bind_rows(
select_cohort, 
select_sex,
select_age,
select_hospaq,
select_icu,
select_comorbidity,
select_comorbidgrp,
)

selec_cont_bind <- 
  bind_rows(
    selected_dur_hosp,
    selected_dur_picu
  )

# publication table Reformat  ----
# Format categorical data
formatted_categorical <- selec_cat_bind %>%
  mutate(Value = paste0(Count, " (", Percentage, "%)")) %>%
  select(Variable, Category, Value)

# Format continuous data
formatted_continuous <- selec_cont_bind %>%
  mutate(Value = paste0("Median: ", Median, ", Range: [", Min, "-", Max, "], SD: ", SD)) %>%
  select(Variable, Value)

# Combine both formatted data into a single dataframe
combined_data_clinical_summary <- bind_rows(formatted_categorical, formatted_continuous)

# Print or view the combined table
print(combined_data_clinical_summary)

# Save the combined summary statistics table to a CSV file
write.table(combined_data_clinical_summary, file = "../../data/cohort_summary_curated/combined_data_clinical_summary.tsv", quote = FALSE, row.names = FALSE, sep = "\t")


# date ----
# df$hosp.adm
# library(lubridate)
# 
# # Define the date columns
# date_columns <- c("hosp.adm", "bc.sampling", "picu.adm", "death.date", "hosp.dis", "picu.dis")  # Replace with the actual column names containing the dates
# 
# # Convert date columns to Date class
# df[date_columns] <- lapply(df[date_columns], function(x) as.Date(x, format = "%d.%m.%y"))
# 
# # Create a list of ggplot objects for date columns
# plots_date <- lapply(date_columns, function(column) {
#   ggplot(df, aes_string(x = column)) +
#     geom_bar(fill = "blue", color = "black") +
#     labs(x = column, y = "Count") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_x_date(date_labels = "%Y-%m-%d")  # Customize date labels if needed
# })
# 
# # Print plots for date columns
# for (i in seq_along(plots_date)) {
#   print(plots_date[[i]])
# }
# 
# 
# # Save the arranged plots to a PDF
# ggsave("plots_dates.pdf", plot = gridExtra::marrangeGrob(grobs = plots_date, ncol = 2, nrow = 2))



