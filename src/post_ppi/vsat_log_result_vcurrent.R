# Load the necessary library
library(dplyr)
library(stringr)
library(ggplot2)

# set version 1-4
# version <- "v1"
# version <- "v2"
# version <- "v3"
version <- "v4"

file <- c(paste("../../data/post_ppi/vsat_log_", version, "_result", sep = ""))
file_suffix <- paste("vsat_log_result_", version, sep = "")

skat_log_read <- read.csv(paste(file, ".txt", sep = ""), sep = ":", header = TRUE)
skat_log_read$rowid <-  rownames(skat_log_read)

skat_log_read_genes <- read.csv(paste(file, "_genes.txt", sep = ""),
                                sep = ":",
                                header = FALSE)
skat_log_read_genes$rowid <-  rownames(skat_log_read_genes)
skat_log_read_genes$Genes <- skat_log_read_genes$V3
skat_log_read_genes <- skat_log_read_genes |> dplyr::select(-V3)
skat_log_read_genes <- skat_log_read_genes |>
  mutate(V1 = str_replace(V1, "\\.log$", "")) |>
  rename(MCL_ID = V1)

skat_log <- merge(skat_log_read_genes, skat_log_read, by = "MCL_ID")
skat_log <- skat_log |> select(MCL_ID, p_val, gene_count, var_count, Genes, everything())
skat_log$MCL_ID <- as.numeric(skat_log$MCL_ID)

# clean up
rm(list = setdiff(ls(), c("skat_log", "file_suffix", "version")))

skat_log <- skat_log |>
  filter(gene_count > 5) |> # v4
  filter(var_count < 100) #v4

# psig threshold
MCL_ID_count <- unique(skat_log$MCL_ID)
psig <- .05 / length(MCL_ID_count)

# Remove quotes from "Name" column
skat_log$Genes <- gsub('"', '', skat_log$Genes)

# remove leading/trailing spaces and split by spaces
skat_log$Genes_clean <- str_trim(skat_log$Genes)
skat_log$gene_count <- sapply(str_split(skat_log$Genes_clean, "\\s+"), length)

skat_log |>
  filter(p_val < psig) |>
  select(Genes, gene_count, MCL_ID, p_val) |>
  arrange(p_val)

# drop qc faile sets
qc <- read.table("../post_ppi/vsat_qc_id.txt")
skat_log <- skat_log %>% filter(!MCL_ID %in% qc$V1)

# Plot the data with highlighted color
pskat <- ggplot(skat_log, aes(x = MCL_ID, y = -log10(p_val))) +
  geom_hline(yintercept = -log10(psig)) +
  geom_point()

# Show the plot
pskat
ggsave(paste("../../data/post_ppi/", file_suffix, ".plot1.pdf", sep = "") ,
       plot = pskat)


pskat2 <- pskat  +
  geom_text(aes(label = ifelse(psig > p_val, MCL_ID, "")), vjust = 1.5, size = 5)

pskat2
ggsave(paste("../../data/post_ppi/", file_suffix, ".plot2.pdf", sep = "") ,
       plot = pskat2)


# Define the colors for highlight and non-highlight
highlight_color <- "red"
default_color <- "black"
dim_color <- "grey"

# Add a new column for color based on the condition
skat_log$color <- ifelse(-log10(skat_log$p_val) < -log10(psig),
                         default_color,
                         highlight_color)
skat_log$color[skat_log$gene_count < 5] <- dim_color

pskat3 <- ggplot(skat_log, aes(
  x = MCL_ID,
  y = -log10(p_val),
  color = color
)) +
  geom_hline(yintercept = -log10(psig)) +
  geom_point() +
  # scale_color_manual(values = c(default_color, highlight_color)) + # v4
  scale_color_manual(values = c(dim_color, default_color, highlight_color)) + # v1
  guides(color = "none")

pskat3 <-
  pskat3 +
  geom_text(aes(label = ifelse(psig > p_val, MCL_ID, "")), vjust = 1.5, size = 5)
pskat3
ggsave(paste("../../data/post_ppi/", file_suffix, ".plot3.pdf", sep = "") ,
       plot = pskat3)

library(ggrepel)
skat_log2_wrapped <- skat_log |>
  filter(psig > p_val) |>
  mutate(Genes_clean = str_wrap(Genes_clean, width = 20))

pskat4 <- pskat3 +
  geom_label_repel(
    data = skat_log2_wrapped |> filter(psig > p_val) |> filter(gene_count > 4),
    aes(
      label = Genes_clean,
      x = MCL_ID,
      y = -log10(p_val)
    ),
    label.padding = unit(0.15, "lines"),
    label.size = 0,
    direction = "y",
    force = 20,
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf,
    nudge_y       = 3
  ) + coord_cartesian(clip = "off") +
  guides(color = "none")

pskat4
ggsave(paste("../../data/post_ppi/", file_suffix, ".plot4.pdf", sep = "") ,
       plot = pskat4)

pskat5 <- skat_log |>
  ggplot(aes(x = AC, y = AC)) +
  geom_point()
pskat5
ggsave(paste("../../data/post_ppi/", file_suffix, ".plot5.pdf", sep = "") ,
       plot = pskat5)

pskat6 <- skat_log |>
  ggplot(aes(x = AC)) +
  geom_histogram()
pskat6
ggsave(paste("../../data/post_ppi/", file_suffix, ".plot6.pdf", sep = "") ,
       plot = pskat6)

# save sig gene list ----
skat_log_psig <-
  skat_log |>
  filter(p_val < psig)

skat_log_psig$Genes_clean
# Perform the replacement operation
skat_log_psig$Genes_clean_csv <- str_replace_all(skat_log_psig$Genes_clean, " ", ",")
skat_log_psig_csv <- skat_log_psig |> dplyr::select(MCL_ID, Genes_clean_csv)

# Write dataframe to a text file
write.table(
  skat_log_psig,
  "../../data/post_ppi/skat_log_psig.txt",
  sep = "\t",
  row.names = FALSE
)
write.table(
  skat_log_psig_csv,
  "../../data/post_ppi/skat_log_psig_csv.txt",
  sep = "\t",
  row.names = FALSE
)

# Not useds ----

# Calculate IQR
Q1 <- quantile(skat_log$AC, 0.25)
Q3 <- quantile(skat_log$AC, 0.75)
IQR <- Q3 - Q1

# Filter out the outliers
# skat_log <- skat_log |> filter(AC >= (Q1 - 1.5 * IQR) & AC <= (Q3 + 1.5 * IQR))
# skat_log <- skat_log |> filter(AC >= (Q1 - 1 * IQR) & AC <= (Q3 + 1 * IQR))


# direction ----
# Let's see if we have any strange direction of effect.
# Clean data to do a simple fisher test on heterozygous cariage.


# This function counts how many columns are used per row, then populate new columns
library(tidyr)
df_clean <- skat_log |>
  pivot_longer(cols = starts_with("AC"),
               names_to = "AC",
               values_to = "value") |>
  group_by(MCL_ID) |>
  na.omit() |>
  mutate(row_number = row_number(), # Add row number within each group
         count = n()) |>
  ungroup()

df_clean <- df_clean |>
  mutate(
    AC = case_when(
      count == 4 & row_number <= 2 ~ "case",
      count == 4 & row_number > 2 ~ "control",
      count == 6 & row_number <= 3 ~ "case",
      count == 6 & row_number > 3 ~ "control",
      count == 5 & row_number <= 3 ~ "case",
      count == 5 & row_number > 3 ~ "control",
      TRUE ~ AC
    ),
    AC = case_when(
      (AC == "case" |
         AC == "control") &
        (count == 4 |
           count == 6) &
        row_number %% 2 == 1 ~ paste(AC, "0", sep = ""),
      (AC == "case" |
         AC == "control") &
        (count == 4 |
           count == 6) &
        row_number %% 2 == 0 ~ paste(AC, "1", sep = ""),
      (AC == "case" |
         AC == "control") &
        count == 5 &
        row_number <= 3 ~ paste(AC, row_number - 1, sep = ""),
      (AC == "case" |
         AC == "control") &
        count == 5 &
        row_number > 3 ~ paste(AC, row_number - 4, sep = ""),
      TRUE ~ AC
    )
  )

df_clean <- df_clean |>
  separate(
    AC,
    into = c("group", "allele"),
    sep = "(?<=[a-z])(?=[0-9])",
    remove = TRUE
  ) |>
  pivot_wider(names_from = "group", values_from = "value")

df_clean <-
  df_clean |> select(MCL_ID, allele, row_number, case, control)

# Reshape the data
df_clean <- df_clean |>
  ungroup() |>
  pivot_longer(
    cols = c(case, control),
    names_to = "group",
    values_to = "count"
  ) |>
  drop_na(count)

# Fill the data
df_clean <- df_clean |>
  group_by(MCL_ID, allele) |>
  complete(group = c("case", "control"))

# for simplicity lets just use heterozygous
df_clean <- df_clean |> filter(row_number < 5)

# Perform Fisher's Exact Test for each MCL_ID
library(purrr)
library(broom)

results <- df_clean |>
  group_by(MCL_ID) |>
  do(tidy(fisher.test(matrix(
    .$count, nrow = 2, byrow = TRUE
  )))) |>
  as.data.frame()

head(results)

# The Fisher's Exact Test in R returns an odds ratio by default, which is stored in the estimate column
results <- results |>
  mutate(direction_of_effect = if_else(estimate > 1, "Increase in controls", "Increase in cases"))

skat_log <- merge(skat_log, (results |> select(direction_of_effect, MCL_ID)), by = "MCL_ID")
skat_log <- skat_log |> select(direction_of_effect, everything())

skat_log <-
  skat_log |> filter(direction_of_effect == "Increase in cases")

ggplot(results, aes(x = MCL_ID, y = -log10(p.value))) +
  geom_point()

ggplot(results, aes(x = -log10(p.value), y = -log10(p.value))) +
  geom_point()

results |>
  ggplot(aes(x = direction_of_effect)) +
  geom_bar(stat = "count")

