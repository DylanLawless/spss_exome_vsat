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
library(patchwork)

# load VSAT report
df_report_sample_vsat <- readRDS("../../data/ACMGuru_post_ppi/df_report_sample_vsat_22_586.Rds")

geneset_MCL_ID <- c(22, 586)

df_report_main_text <- readRDS(paste0("../../data/ACMGuru_post_ppi/df_report_main_text_", paste(geneset_MCL_ID, collapse="_"), ".Rds"))

# load single case report
# df_report_main_text <- readRDS(file="../../data/singlecase/df_report_main_text.Rds")

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
saveRDS(df, file = "../../data/cohort_summary_curated/cohort_summary_curated_r_df.Rds")

# text summary description Hmisc ----
df <- df |> dplyr::select(
  # -sample.id,
  -V3,
  -V4,
  -V5)

df <- merge(df_report_sample_vsat, df, by="sample.id", all=T)

df$group <- ifelse(is.na(df$group), "not_contributer", df$group)

# continuous: hist plots -----
# Convert the data from wide to long
# df_long <- df |>
#   dplyr::select(which(sapply(df, is.numeric))) |>  # Select numeric columns
#   gather(group, key = "variable", value = "value")  # Convert from wide to long format

df_long <- df |>
  dplyr::select(group, which(sapply(df, is.numeric))) |>  # Select numeric columns
  gather(2:23, key = "variable", value = "value")  # Convert from wide to long format


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
    # geom_histogram(bins = bins, alpha=0.5, position = "dodge", color = "black") +
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
ggsave("../../images/cohort_summary_curated/cohort_plots_post_ppi_continuous.pdf" ,plot = p_combined1, height = 10, width = 10)


# categorical: bar plots -----
# Convert the data from wide to long
# df_long <- df |>
#   dplyr::select(which(sapply(df, is.character))) |>  # Select character columns
#   gather(key = "variable", value = "value")  # Convert from wide to long format

df_long <- df |>
  dplyr::select(group, which(sapply(df, is.character))) |>  # Select numeric columns
  gather(3:31, key = "variable", value = "value")  # Convert from wide to long format

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

# Define a function to create a bar plot
create_bar <- function(data, variable_name) {
  # Create the plot
  data |>
  ggplot(aes(x = value, fill = group)) +
    geom_bar(color = "black", position="dodge") +
    # geom_density(alpha=0.8, stat = "bin", bins = bins )+ 
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
ggsave("../../images/cohort_summary_curated/cohort_plots_post_ppi_categorical.pdf", plot = p_combined2, height = 10, width = 10)

# patchwork ----
# plot1 + (plot2 + plot3) + plot_layout(ncol = 1)
patch1 <- (p_combined2 | p_combined1) + plot_annotation(tag_levels = 'A')

# red is VSAT contributer
ggsave("../../images/cohort_summary_curated/cohort_plots_post_ppi_cat_con_redVSAT.pdf", plot = patch1, height = 10, width = 20)

df |> count(group)

# kable latex tables ----
library(knitr)

df_report_main_text |>
# df_summaries |> 
  group_by(genotype) |>
  summarise(variants = n()) |>
  kable("latex", booktabs = TRUE)

df_report_main_text |>
  # df_summaries |> 
  group_by(Inheritance) |>
  summarise(n())|>
  kable("latex", booktabs = TRUE)

df_report_main_text |>
  # df_summaries |> 
  ungroup() |>
  group_by(ACMG_total_score, Consequence, IMPACT) |>
  summarise(unique_variants = n()) |>
  arrange(desc(ACMG_total_score), IMPACT, unique_variants) |>
  kable("latex", booktabs = TRUE)


# df_report_main_text |>
#   filter(ACMG_total_score >= 6) |> 
#   dplyr::select(SYMBOL, 
#                 Consequence, 
#                 sample.id,
#                 genotype,
#                 Inheritance) |>
#   kable("latex", booktabs = TRUE)


# collapse Inheritance column 
df_report_main_text |>
  filter(ACMG_total_score >= 1) |>
  group_by(SYMBOL, Consequence, sample.id, genotype) |>
  summarise(Inheritance = paste(Inheritance, collapse = ", ")) |>
  dplyr::select(-Inheritance) |> # drop inheritance since we have none
  ungroup() |>
  kable("latex", booktabs = TRUE)

# collapse Inheritance column 
df_report_main_text |>
  filter(ACMG_total_score >= 6) |>
  group_by(SYMBOL, HGVSp, Consequence, sample.id, genotype) |>
  summarise(Inheritance = paste(Inheritance, collapse = ", ")) |>
  dplyr::select(-Inheritance) |> # drop inheritance since we have none
  ungroup() |>
  kable("latex", booktabs = TRUE)

# gene list
unique_symbols <- df_report_main_text %>%
  dplyr::pull(SYMBOL) %>% 
  unique() %>%
  paste(collapse = ", ")

unique_symbols

nrow(df_report_main_text)
nrow(df_report_main_text |> dplyr::select(rownames) |> unique())  
nrow(df_report_main_text |> dplyr::select(rownames, sample.id) |> unique())  
nrow(df_report_main_text |> dplyr::select(rownames, SYMBOL) |> unique())  

