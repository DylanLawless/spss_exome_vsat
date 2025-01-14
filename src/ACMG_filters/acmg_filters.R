# acmg_filters ----
print("We are adding PS3 now")

# BA1 ----
df_acmg |> filter(ACMG_label == "BA1") |> dplyr::select(Criteria)
# Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium
df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
gnomad_BA1 <- .05 
df$ACMG_BA1 <- NA
df$ACMG_BA1 <- ifelse(df$gnomAD_AF > gnomad_BA1, "BA1", NA)
df <- df %>% dplyr::select(ACMG_BA1, everything())
# df |> filter(ACMG_BA1 == "BA1") |> dplyr::select(gnomAD_AF)

# BS1 ----
df_acmg |> filter(ACMG_label == "BS1") |> dplyr::select(Criteria)
# Allele frequency is greater than expected for disorder
# duplicate previous frequency unless other info
df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
gnomad_BS1 <- .05 
df$ACMG_BS1 <- NA
df$ACMG_BS1 <- ifelse(df$gnomAD_AF > gnomad_BS1, "BS1", NA)
df <- df %>% dplyr::select(ACMG_BS1, everything())
# df |> filter(ACMG_BS1 == "BS1") |> dplyr::select(gnomAD_AF)

# BS2 ----
df_acmg |> filter(ACMG_label == "BS2") |> dplyr::select(Criteria)
# Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an early age
# already removed
df$ACMG_BS2 <- NA
df <- df %>% dplyr::select(ACMG_BS2, everything())
# df |> filter(ACMG_BS2 == "BS2") |> dplyr::select(gnomAD_AF) 

# BS3 ----
df_acmg |> filter(ACMG_label == "BS3") |> dplyr::select(Criteria)
# Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing
# already removed
df$ACMG_BS3 <- NA
df <- df %>% dplyr::select(ACMG_BS3, everything())
# df |> filter(ACMG_BS3 == "BS3") |> dplyr::select(gnomAD_AF) 

# BS4 ----
df_acmg |> filter(ACMG_label == "BS4") |> dplyr::select(Criteria)
# Lack of segregation in affected members of a family
# not available
df$ACMG_BS4 <- NA
df <- df %>% dplyr::select(ACMG_BS4, everything())
# df |> filter(ACMG_BS4 == "BS4") |> dplyr::select(gnomAD_AF) 

# BS5 ----
df_acmg |> filter(ACMG_label == "BS5") |> dplyr::select(Criteria)
# The user has additional (value) strong benign evidence
# not available
df$ACMG_BS5 <- NA
df <- df %>% dplyr::select(ACMG_BS5, everything())
# df |> filter(ACMG_BS5 == "BS5") |> dplyr::select(gnomAD_AF) 

# BP1 ----
df_acmg |> filter(ACMG_label == "BP1") |> dplyr::select(Criteria)
# Missense variant in a gene for which primarily truncating variants are known to cause disease
# not available
df$ACMG_BP1 <- NA
df <- df %>% dplyr::select(ACMG_BP1, everything())
# df |> filter(ACMG_BP1 == "BP1") |> dplyr::select(gnomAD_AF) 

# BP2 ----
df_acmg |> filter(ACMG_label == "BP2") |> dplyr::select(Criteria)
# Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern
df$ACMG_BP2 <- NA
df <- df %>% dplyr::select(ACMG_BP2, everything())

# BP3 ----
df_acmg |> filter(ACMG_label == "BP3") |> dplyr::select(Criteria)
# In-frame deletions/insertions in a repetitive region without a known function
df$ACMG_BP3 <- NA
df <- df %>% dplyr::select(ACMG_BP3, everything())
# Custom logic for BP3

# BP4 ----
df_acmg |> filter(ACMG_label == "BP4") |> dplyr::select(Criteria)
# Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary, splicing impact, etc.)
df$ACMG_BP4 <- NA
df <- df %>% dplyr::select(ACMG_BP4, everything())
# Custom logic for BP4

# BP5 ----
df_acmg |> filter(ACMG_label == "BP5") |> dplyr::select(Criteria)
# Variant found in a case with an alternate molecular basis for disease
df$ACMG_BP5 <- NA
df <- df %>% dplyr::select(ACMG_BP5, everything())
# Custom logic for BP5

# BP6 ----
df_acmg |> filter(ACMG_label == "BP6") |> dplyr::select(Criteria)
# Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation
df$ACMG_BP6 <- NA
df <- df %>% dplyr::select(ACMG_BP6, everything())
# Custom logic for BP6

# BP7 ----
df_acmg |> filter(ACMG_label == "BP7") |> dplyr::select(Criteria)
df$ACMG_BP7 <- NA
df <- df %>% dplyr::select(ACMG_BP7, everything())
# Custom logic for BP7

# BP8 ----
df_acmg |> filter(ACMG_label == "BP8") |> dplyr::select(Criteria)
# A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved
# already removed
df$ACMG_BP8 <- NA
df <- df %>% dplyr::select(ACMG_BP8, everything())
# Custom logic for BP8

# PVS1 ----
# PVS1 are null variants where IMPACT=="HIGH" and inheritance match, in gene where LoF cause disease.
df$ACMG_PVS1 <- NA
df <- df %>% dplyr::select(ACMG_PVS1, everything())
df$ACMG_PVS1 <- ifelse(df$IMPACT == "HIGH" & df$genotype == 2, "PVS1", NA) # homozygous
df$ACMG_PVS1 <- ifelse(df$IMPACT == "HIGH" & df$Inheritance == "AD", "PVS1", df$ACMG_PVS1) # dominant
# df |> filter(ACMG_PVS1 == "PVS1")

# include comp_het if both HIGH impact. WARNING NOT PHASE CHECKED
df <- df |>
  group_by(sample, SYMBOL) |>
  mutate(ACMG_PVS1 = ifelse(ACMG_PVS1 == "PVS1", "PVS1", 
                            ifelse(sum(IMPACT == "HIGH" & comp_het_flag == 1) >= 2 & IMPACT == "HIGH", "PVS1", ACMG_PVS1))) %>%
  ungroup() 
# df |> filter(ACMG_PVS1 == "PVS1")

# PS1 ----
# df$ACMG_PS1 <- NA
# df <- df %>% dplyr::select(ACMG_PS1, everything())
# df$ACMG_PS1 <- ifelse(df$CLIN_SIG == "pathogenic", "PS1", NA)

# Updated
# PS1 Same amino acid change as a previously established pathogenic variant regardless of nucleotide change. Note to keep splicing variant as PSV1 (these are covered by IPACT HIGH).
# Keep any thing from CLIN_SIG or anything that has a ClinVar_CLNDN expect for NA/benign.
# ? Allow CLIN_SIG %in% Conflicting_interpretations_of_pathogenicity for more modest screening
df$ACMG_PS1 <- NA
df <- df %>% dplyr::select(ACMG_PS1, everything())

# This code essentially labels entries as "PS1" where the clinical significance is marked as "pathogenic", the ClinVar_CLNDN.y is not missing, and it does not include undesirable terms unless explicitly containing "not_specified&not_provided".

# df$ACMG_PS1 <- ifelse(df$CLIN_SIG %in% c("pathogenic") & 
												# !is.na(df$ClinVar_CLNDN.y) & 
												# !grepl("benign", df$CLIN_SIG) &
												# !df$ClinVar_CLNDN.y %in% c("not_provided", "not_specified", "not_specified&not_provided") |
												# grepl("not_specified&not_provided", df$ClinVar_CLNDN.y), "PS1", NA)

df$ACMG_PS1 <- ifelse(
  # df$CLIN_SIG %in% c("pathogenic") &
    grepl("pathogenic", df$CLIN_SIG) #&
    # !is.na(df$ClinVar_CLNDN.y) &
    # !grepl("benign", df$CLIN_SIG) &
    # !grepl("not_specified|not_provided", df$ClinVar_CLNDN.y),
    , "PS1", NA)

# print("!!! PS1 test !!\n\n")
# print(df$CLIN_SIG |> unique())
# print(df$ClinVar_CLNDN.y |> unique())
# print(df |> dplyr::select(CLIN_SIG, ClinVar_CLNDN.y) |> unique())

# df |> filter(ACMG_PS1 == "PS1")



# PS2 skip ----
# PS2 De novo (both maternity and paternity confirmed) in a patient with the disease and no family history
# Skip due to no parental genetics.

print("We are adding PS3 now")
# PS3 ----
# PS3 Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.
# df$ACMG <- ifelse(uniprot == "pathogenic" |
# 							pubmed == "pathogenic",
# 						"PS3")
# df$ACMG_PS3 <- ifelse(! is.na(df$Inheritance), "PS3", NA) 

# Create a new column 'ACMG_PS3' and initialise with NA
df$ACMG_PS3 <- NA

# Set the label 'PS3' based on multiple criteria
df$ACMG_PS3[(
  (df$genotype >= 1 & df$Inheritance == "AD") |
    (df$genotype >= 1 & df$Inheritance == "AD/AR") |
    (df$genotype == 2 & df$Inheritance == "AR") |
    (df$genotype == 2 & df$Inheritance == "XL") |
    (df$genotype == 2 & df$Inheritance == "XLR")
)] <- "PS3"

# PS4 skip ----
# The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls
# Skip do to statistical analysis separately

# PS5 ----
# The user has additional (value) strong pathogenic evidence
df$ACMG_PS5 <- NA
df <- df |> dplyr::select(ACMG_PS5, everything())

# comp_het with at least 1 HIGH impact. WARNING NOT PHASE CHECKED
df <- df %>%
  group_by(sample, SYMBOL) %>%
  mutate(ACMG_PS5 = ifelse(any(IMPACT == "HIGH") & (n() > 1), "PS5", ACMG_PS5)) %>%
  ungroup()
# df |> filter(ACMG_PS5 == "PS5")

# PM2 ----
# Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium

df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
gnomad_max <- 1e-6 # round down to approx less than 1 on gnomad.
df$ACMG_PM2 <- NA
df <- df %>% dplyr::select(ACMG_PM2, everything())
df$ACMG_PM2 <- ifelse(df$gnomAD_AF < gnomad_max, "PM2", NA)
# df |> filter(ACMG_PM2 == "PM2")

# PM3 ----
# For recessive disorders, detected in trans with a pathogenic variant
# some redundancy with our PS5 since our rare disease cohort filtering call IMPACT==HIGH equates pathogenic
df$ACMG_PM3 <- NA
df <- df %>% dplyr::select(ACMG_PM3, everything())
df <- df %>%
  group_by(sample, SYMBOL) %>%
  mutate(ACMG_PM3 = ifelse(comp_het_flag == 1 & (ACMG_PS1 == "PS1" | ACMG_PS5 == "PS5"), 
                           "PM3", ACMG_PM3)) %>%
  ungroup()
# df |> filter(ACMG_PM3 == "PM3")

# PP3 in silico ----
# Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.). CUSTOM: assigned if >=3 thresholds passed. 

# In-Silico Predictions: VarSome now implements the ClinGen recommendations from Evidence-based calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for clinical use of PP3/BP4 criteria: # Only one engine at a time is used, depending on availability of data, in order: MitoTip & MitImpact, MetaRNN, CADD (Premium only), DANN (if CADD is not available). # The maximum strength allowed for rules PP3 & BP4 is Strong, even if there may be evidence for Very Strong, with the exception of variants that are predicted splicing (ie: similar to PVS1). # The strength is limited to Supporting, if there's Moderate evidence from rules PM1 or PM5. # Splice prediction (scSNV) is given priority over the other in-silico predictions. # conservation is used for some low-sensitivity variant types, or if no other in-silico prediction is available. Please refer to PP3 and BP4 for more specific detail.
# 
# df <- df |> separate(SIFT, into = c("SIFT_label", "SIFT_score"), sep = "\\(", remove = TRUE) |>
#   mutate(SIFT_score = str_replace(SIFT_score, "\\)", "")) 
# 
# df <- df |> separate(PolyPhen, into = c("PolyPhen_label", "PolyPhen_score"), sep = "\\(", remove = TRUE) |>
#   mutate(PolyPhen_score = str_replace(PolyPhen_score, "\\)", "")) 


# Check if the SIFT columns have already been separated, if not then separate them
if("SIFT" %in% colnames(df) && !("SIFT_label" %in% colnames(df) && "SIFT_score" %in% colnames(df))) {
  df <- df %>%
    separate(SIFT, into = c("SIFT_label", "SIFT_score"), sep = "\\(", remove = TRUE) %>%
    mutate(SIFT_score = str_replace(SIFT_score, "\\)", ""))
}

# Check if the PolyPhen columns have already been separated, if not then separate them
if("PolyPhen" %in% colnames(df) && !("PolyPhen_label" %in% colnames(df) && "PolyPhen_score" %in% colnames(df))) {
  df <- df %>%
    separate(PolyPhen, into = c("PolyPhen_label", "PolyPhen_score"), sep = "\\(", remove = TRUE) %>%
    mutate(PolyPhen_score = str_replace(PolyPhen_score, "\\)", ""))
}


df$CADD_PHRED <- as.numeric(df$CADD_PHRED)
df$REVEL_rankscore <- as.numeric(df$REVEL_rankscore)

# Define your conditions
cond_CADD_PHRED <-            df$CADD_PHRED >= 30
cond_REVEL_rankscore <-       df$REVEL_rankscore > .5
cond_MetaLR_pred <-           df$MetaLR_pred == "D"
cond_MutationAssessor_pred <- df$MutationAssessor_pred == "H"
cond_SIFT_label <-            df$SIFT_label == "deleterious"
cond_PolyPhen_label <-        df$PolyPhen_label == "probably_damaging"

# Initialize the ACMG_PP3 column with NA
df$ACMG_PP3 <- NA
df <- df %>% dplyr::select(ACMG_PP3, everything())

# Count the points and store them in ACMG_PP3
df$ACMG_PP3_count <- rowSums(cbind(cond_CADD_PHRED, cond_REVEL_rankscore, cond_MetaLR_pred, 
                                   cond_MutationAssessor_pred, cond_SIFT_label, cond_PolyPhen_label), na.rm = TRUE)

threshold <- 3
df$ACMG_PP3 <- ifelse(df$ACMG_PP3_count >= threshold, "PP3", NA)
df |> filter(ACMG_PP3 == "PP3")

# Remove temporary column
df$ACMG_PP3_count <- NULL

rm(list=setdiff(ls(), c("df",  "df_acmg", "df_acmg_caveat", "geneset_MCL_ID", "file_suffix", "output_directory", "hold", "iuis", "varsome", "ac_count_per_var")))

# PP3 In silico: varsome ----
# Varsome conditions
varsome |>  dplyr::select(Engine) 

# Rename varsome to match our data
names_to_replace <- list(
  c("BayesDel_addAF", "BayesDel_addAF_score"),
  c("BayesDel_noAF", "BayesDel_noAF_score"),
  c("CADD", "CADD_PHRED"),
  c("DANN", "DANN_score"),
  c("EIGEN", "Eigen.raw_coding"),
  c("EIGEN-PC", "Eigen.PC.phred_coding"),
  c("FATHMM", "FATHMM_score"),
  c("FATHMM-MKL", "fathmm.MKL_coding_score"),
  c("FATHMM-XF", "fathmm.XF_coding_score"),
  c("LRT", "LRT_score"),
  c("M-CAP", "M.CAP_score"),
  c("MetaLR", "MetaLR_score"),
  c("MetaSVM", "MetaSVM_score"),
  c("MetaRNN", "MetaRNN_score"),
  c("MutPred", "MutPred_score"),
  c("MutationAssessor", "MutationAssessor_score"),
  c("MutationTaster", "MutationTaster_score"),
  c("phastCons100way_vertebrate", "phastCons100way_vertebrate"),
  c("Polyphen2-HDIV", "Polyphen2_HDIV_score"),
  c("Polyphen2-HVAR", "Polyphen2_HVAR_score"),
  c("PROVEAN", "PROVEAN_score"),
  c("REVEL", "REVEL_score"),
  c("SIFT", "SIFT_score")
)

# Loop over the list and replace the old names with the new names
for (name_pair in names_to_replace) {
  varsome$Engine <- replace(varsome$Engine, varsome$Engine == name_pair[1], name_pair[2])
}

# Not used: BLOSUM DANN DEOGEN2 EVE LIST-S2 M-CAP MVP MaxEntScan MitImpact MitoTip PrimateAI SIFT4G phyloP (PhyloP100Way) scSNV-ADA scSNV-RF

# varsome list of thresholds to tally conditions met
calculate_varsome_score <- function(df, varsome, pathogenic_type) {
  varsome_list <- setNames(varsome[[pathogenic_type]], varsome$Engine)
  
  df[[pathogenic_type]] <- 0
  
  for (engine in names(varsome_list)) {
    if (!(engine %in% names(df))) {
      print(paste(engine, "not found in df, skipping..."))
      next
    }
    
    if (!is.numeric(df[[engine]])) {
      print(paste(engine, "is not numeric, converting..."))
      df[[engine]] <- as.numeric(df[[engine]])
    }
    
    condition <- df[[engine]] >= varsome_list[[engine]]
    
    condition <- tidyr::replace_na(condition, 0)
    
    print(paste(engine, ":", sum(is.na(condition)), "NAs.",
                pathogenic_type, ":", sum(condition, na.rm = TRUE)))
    
    df[[pathogenic_type]] <- df[[pathogenic_type]] + condition
  }
  
  return(df)
}

df <- calculate_varsome_score(df, varsome, "Strong_pathogenic_GE")
df <- calculate_varsome_score(df, varsome, "Moderate_pathogenic_GE")
df <- calculate_varsome_score(df, varsome, "Supporting_pathogenic_GE")

df <- df |> dplyr::select(ends_with("_pathogenic_GE"), everything())

# distributions and thresholds 
# library(tidyverse)
varsome_thresholds <- varsome %>%
  dplyr::select(Engine, ends_with("_pathogenic_GE")) %>%
  pivot_longer(cols = -Engine,
               names_to = "pathogenicity",
               values_to = "threshold")

common_cols <- intersect(varsome_thresholds$Engine, names(df))

df_long <- df |>
  dplyr::select(all_of(common_cols)) |>
  pivot_longer(cols = all_of(common_cols),
               names_to = "Engine",
               values_to = "Score")


# The Engine names are too long for our plot. Named vector where names are new (long) names and values are old (short) names
name_mapping <- setNames(sapply(names_to_replace, `[[`, 1), sapply(names_to_replace, `[[`, 2))
df_long$Engine_short <- name_mapping[df_long$Engine]

p.pathogenicity_distributions_engines <- df_long |>
  # ggplot(aes(x = NormScore, fill=..x..)) +
  ggplot(aes(x = Score, fill=..x..)) +
  geom_histogram(
    #color="black"
  ) +
  facet_wrap(~Engine_short, scales = "free") +
  theme_minimal() +
  xlab("in silico prediction score") +
  ylab("No. qualifying variants")+ 
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'bamako', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.pathogenicity_distributions_engines
ggsave(paste("../../images/", output_directory ,file_suffix, "pathogenicity_distributions_engines.pdf", sep = "") ,plot = p.pathogenicity_distributions_engines, width = 9, height = 5)

# Append a suffix to the pathogenicity column in varsome_thresholds
varsome_thresholds$pathogenicity <- paste0(varsome_thresholds$pathogenicity, "_threshold")

# Pivot varsome_thresholds to wide format
varsome_thresholds_wide <- varsome_thresholds |>
  pivot_wider(names_from = pathogenicity, values_from = threshold)

# Join with df_long
df_long <- left_join(df_long, varsome_thresholds_wide, by = "Engine")

Strong_pathogenic_GE_threshold <- 3
Moderate_pathogenic_GE_threshold <- 4
Supporting_pathogenic_GE_threshold <- 10

threshold_results <- df |>
  summarize(
    total = n(),
    strong_threshold = Supporting_pathogenic_GE_threshold,
    strong_pass = sum(Strong_pathogenic_GE >= Strong_pathogenic_GE_threshold, na.rm = TRUE),
    strong_percent = (strong_pass / total) * 100,
    
    total = n(),
    moderate_threshold = Moderate_pathogenic_GE_threshold,
    moderate_pass = sum(Moderate_pathogenic_GE >= Moderate_pathogenic_GE_threshold, na.rm = TRUE),
    moderate_percent = (moderate_pass / total) * 100,
    
    total = n(),
    supporting_threshold = Supporting_pathogenic_GE_threshold,
    supporting_pass = sum(Supporting_pathogenic_GE >= Supporting_pathogenic_GE_threshold, na.rm = TRUE),
    supporting_percent = (supporting_pass / total) * 100
  )

threshold_results_long <- threshold_results |>
  pivot_longer(everything(),
               names_to = "Measurement",
               values_to = "Value")

threshold_results_long %>%
  mutate(Value = round(Value, 2)) %>% 
  kable("latex", booktabs = TRUE)

# assign PP3 ----
df$ACMG_PP3 <- ifelse(
  df$Strong_pathogenic_GE >= Strong_pathogenic_GE_threshold |
    df$Moderate_pathogenic_GE >= Moderate_pathogenic_GE_threshold |
    df$Supporting_pathogenic_GE >= Supporting_pathogenic_GE_threshold, 
  "PP3", 
  df$ACMG_PP3
)

df |> filter(ACMG_PP3 == "PP3")

# independent fill scales ----
# Preparing the data
df_filtered <- df_long %>%
  group_by(Engine) %>%
  filter(!all(is.na(Score))) %>%
  ungroup()

# Creating a list of plots for each group
p.list <- lapply(sort(unique(df_filtered$Engine_short)), function(i) {
  
  df_group <- df_filtered[df_filtered$Engine_short==i, ]
  
  df_group |>
    ggplot(aes(x = Score, fill=..x..)) +
    geom_histogram(bins = 30) +
    theme_minimal(base_size = 8) +
    labs(subtitle =i) +
    xlab("") +
    ylab("") +
    geom_vline(aes(xintercept = Supporting_pathogenic_GE_threshold),
               linetype = "dashed", color = "#eeaf61") +
    geom_vline(aes(xintercept = Moderate_pathogenic_GE_threshold),
               linetype = "dashed", color = "#ee5d6c") +
    geom_vline(aes(xintercept = Strong_pathogenic_GE_threshold),
               linetype = "dashed", color = "#6a0d83")+ 
    guides(fill=FALSE) +
    scale_fill_scico(palette = 'bamako', direction = 1)
  
})

# add legend
df_empty <- data.frame()
legend_only_plot <- 
  ggplot(df_empty) +
  geom_vline(aes(xintercept = 3, color = "#eeaf61")) +
  geom_vline(aes(xintercept = 2, color = "#ee5d6c")) +
  geom_vline(aes(xintercept = 1, color = "#6a0d83")) + 
  scale_color_identity("", 
                       breaks = c("#eeaf61", "#ee5d6c", "#6a0d83"), 
                       labels = c("Supporting pathogenic",
                                  "Moderate pathogenic",
                                  "Strong pathogenic"), 
                       guide = "legend") +
  theme_void() +
  guides(color = guide_legend(reverse = TRUE))

legend <- get_legend(legend_only_plot)
n_plots <- length(p.list)
ncol <- 5
nrow <- ceiling(n_plots / ncol) 
plot_list <- c(p.list, 
               rep(list(NULL), nrow * ncol - n_plots - 1), 
               list(legend))

# Arrange all plots together
p.pathogenicity_distributions_engines_threshold <-  
  annotate_figure(
    ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow),
    left = textGrob("No. qualifying variants", rot = 90, vjust = 1 ),
    bottom = textGrob("in silico prediction score" )
  )
p.pathogenicity_distributions_engines_threshold
ggsave(paste("../../images/", output_directory ,file_suffix, "pathogenicity_distributions_engines_threshold.pdf", sep = "") ,plot = p.pathogenicity_distributions_engines_threshold)


ggsave(paste("../../images/", output_directory,  file_suffix, "pathogenicity_distributions_engines_threshold.pdf", sep = "") ,plot = p.pathogenicity_distributions_engines_threshold)
 
# 
# thresholds passed 
labels <- c( Strong_pathogenic_GE="Strong", Moderate_pathogenic_GE="Moderate", Supporting_pathogenic_GE="Supporting")

p.pathogenicity_distributions <- df |> 
  tidyr::pivot_longer(cols = ends_with("_pathogenic_GE"),
                      names_to = "pathogenicity",
                      values_to = "varsome_score") |> 
  mutate(pathogenicity = fct_relevel(pathogenicity, 
                                     "Strong_pathogenic_GE", "Moderate_pathogenic_GE", "Supporting_pathogenic_GE")) |> 
  ggplot(aes(x = varsome_score, fill=..x..)) +
  geom_histogram(binwidth = 1, color="black") +
  geom_text_repel(stat='count', color = "black", 
                  box.padding = 0.5, max.overlaps = Inf,
                  # padding = unit(0.5, "lines"),
                  # nudge_y = 0.05,  
                  nudge_x = .0,
                  nudge_y = .1,
                  direction = "y",
                  aes(label= ifelse(..count.. < 500, ..count.., ''))
  ) +
  facet_grid(pathogenicity ~ ., labeller=labeller(pathogenicity = labels)) +
  theme_minimal() +
  xlab("Pathogenicity\nthresholds passed") +
  ylab("No. variants")+ 
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'bamako', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.pathogenicity_distributions
ggsave(paste("../../images/", output_directory ,file_suffix, "pathogenicity_distributions.pdf", sep = "") ,plot = p.pathogenicity_distributions)


# acmg tally  ----
# List of all ACMG labels
# acmg_labels <- c("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5", "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6", "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4")

# Transform 'Evidence_type' to 'P' for pathogenic and 'B' for benign
df_acmg$code_prefix <- ifelse(df_acmg$Evidence_type == "pathogenicity", "P", "B")

# Create the ACMG code by combining the new prefix and the label, prepending 'ACMG_'
df_acmg$ACMG_code <- paste0("ACMG_", df_acmg$code_prefix, df_acmg$label)

acmg_labels <- df_acmg$ACMG_code
print(acmg_labels)

# df_acmg$ACMG_label
names(df) |> head(30) |> as.character()

# Check if each ACMG column exists, if not create it and fill with NA
for (acmg_label in acmg_labels) {
  if (!acmg_label %in% names(df)) {
    print("missing label")
    df[[acmg_label]] <- NA
  }
}

# Then use coalesce to find the first non-NA ACMG label
df$ACMG_highest <- dplyr::coalesce(!!!df[acmg_labels])
df <- df %>% dplyr::select(ACMG_highest, everything())

# Count the number of non-NA values across the columns
df$ACMG_count <- rowSums(!is.na(df[, acmg_labels ]))
df <- df %>% dplyr::select(ACMG_count, everything())
# df$ACMG_count[df$ACMG_count == 0] <- NA

p.criteria_count_each_gene <- df |> 
  filter(ACMG_count > 1) |>
  ggplot(aes(y = ACMG_count, x = SYMBOL)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle=45, hjust=1, vjust=1)) +
  xlab("\nGene symbol") +
  ylab("ACMG criteria count (>1)")
p.criteria_count_each_gene
ggsave(paste("../../images/", output_directory ,file_suffix, "criteria_count_each_gene.pdf", sep = "") ,plot = p.criteria_count_each_gene )

# as table
df |> 
  filter(ACMG_count > 1) |>
  dplyr::select(sample, SYMBOL, ACMG_count) |>
  arrange(desc(ACMG_count))

p.criteria_gene_total <- df %>%
  group_by(SYMBOL) |>
  summarise(acmg_count_per_symbol = sum(ACMG_count)) |>
  na.omit() |>
  ggplot(aes(x = acmg_count_per_symbol, fill=..x..) ) +
  geom_histogram(stat="count", binwidth = 1, color="black"
  ) +
  theme_minimal() +
  xlab("No. ACMG criteria (P) variants per gene") +
  ylab("Number of genes") +
  geom_text(stat='count', aes(label=..count.., y=..count..+1), color = "black") + 
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.criteria_gene_total 
ggsave(paste("../../images/", output_directory ,file_suffix, "criteria_gene_total.pdf", sep = "") ,plot = p.criteria_gene_total )

# as table
df |>
  group_by(SYMBOL) |>
  summarise(acmg_count_per_symbol = sum(ACMG_count)) |>
  na.omit() |>
  arrange(desc(acmg_count_per_symbol))

p.variants_per_criteria <- df |> 
  ggplot(aes(x = ACMG_count, fill=..x..)) +
  geom_histogram(binwidth = 1, color="black") +
  xlab("No. ACMG criteria\nassigned (P)") +
  ylab("No. variants") +
  theme_minimal() +
  geom_text(stat='count', aes(label=..count.., y=..count..+20), color = "black") + 
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.variants_per_criteria
ggsave(paste("../../images/", output_directory ,file_suffix, "variants_per_criteria.pdf", sep = "") ,plot = p.variants_per_criteria , width = 9, height = 5)

# Check we only have approx. 1 "casual" variant per sample
p.criteria_per_sample <- df %>%
  group_by(sample) %>%
  summarise(ACMG_count = max(ACMG_count, na.rm = TRUE))  %>%
  ggplot(aes(x = ACMG_count, fill=..x..)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(x = "No. ACMG criteria\nassigned (P)", y = "No. samples") +
  theme_minimal() +
  geom_text(stat='count', aes(label=..count.., y=..count..+10), color = "black") + 
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.criteria_per_sample
ggsave(paste("../../images/", output_directory ,file_suffix, "criteria_per_sample.pdf", sep = "") ,plot = p.criteria_per_sample, width = 9, height = 5)

# as table
df |> 
  group_by(sample, ACMG_count) |>
  tally(n = "count_per_sample") |>
  ungroup() |>
  dplyr::select(-sample) |>
  group_by(ACMG_count) |>
  tally(n = "count_per_sample")

# ACMG Verdict----
# Rules are combined using the point system described in PMID:32720330
# Each rule triggered is assigned a number of points based on the strength of the evidence provided:
# 
# Supporting: 1 point
# Moderate: 2 points
# Strong: 4 points
# Very Strong: 8 points
# A total score is computed as the sum of the points from the pathogenic rules, minus the sum of the points from benign rules.
# 
# The total score is then compared to thresholds to assign the final verdict:
# 	
# Pathogenic if greater than or equal to 10,
# Likely Pathogenic if between 6 and 9 inclusive,
# Uncertain Significance if between 0 and 5,
# Likely Benign if between -6 and -1,
# Benign if less than or equal to -7.

# df <-  df |> dplyr::select("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5", 
#                  "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6", 
#                  "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4",
#                  everything())

df <-  df |> dplyr::select(acmg_labels, everything())

# df <- df %>% dplyr::select(ACMG_highest, ACMG_score, ACMG_count,ACMG_PP3:ACMG_PVS1, ACMG_PS2:ACMG_PP4, everything())

# Define scores for each ACMG label
# acmg_scores <- c("PVS1" = 8,
# 					  "PS1" = 4, "PS2" = 4, "PS3" = 4, "PS4" = 4, "PS5" = 4,
# 					  "PM1" = 2, "PM2" = 2, "PM3" = 2, "PM4" = 2, "PM5" = 2,
# 					  "PP3" = 1)

# Define scores for Pathogenic criteria
pathogenic_scores <- c(
  "PVS1" = 8,
  setNames(rep(4, 5), paste0("PS", 1:5)),
  setNames(rep(2, 7), paste0("PM", 1:7)),
  "PP3" = 1
)

# Define scores for Benign criteria
benign_scores <- c(
  "BA1" = -8,
  setNames(rep(-4, 5), paste0("BS", 1:5)),
  setNames(rep(-1, 8), paste0("BP", 1:8))
)

# Combine both scoring systems into one vector
acmg_scores <- c(pathogenic_scores, benign_scores)

# Print the complete ACMG scoring system
print(acmg_scores)

# Create ACMG_score column by looking up ACMG_highest in acmg_scores
df$ACMG_score <- acmg_scores[df$ACMG_highest]

# If there are any ACMG labels that don't have a corresponding score, these will be NA. You may want to set these to 0.
df$ACMG_score[is.na(df$ACMG_score)] <- 0
df <- df |> dplyr::select(ACMG_score, everything())

# Total ACMG score ----
# List of all ACMG columns
# acmg_columns <- grep("ACMG_P", colnames(df), value = TRUE)
# Update the way columns are selected based on the specific acmg_labels. this is now redundant, switch to "acmg_labels" insted of acmg_columns
# acmg_labels <- colnames(df)[colnames(df) %in% acmg_labels]

# Mutate all ACMG columns
df <- df %>% 
  mutate_at(acmg_labels, function(x) acmg_scores[x])

# Replace NAs with 0 in ACMG columns only
df[acmg_labels] <- lapply(df[acmg_labels], function(x) ifelse(is.na(x), 0, x))

# Calculate total ACMG score
df$ACMG_total_score <- rowSums(df[acmg_labels])

df <- df |> dplyr::select(ACMG_total_score, everything())



# PS4 method ----

# df$ACMG_PS4 <- NA
# df <- df %>% dplyr::select(ACMG_PS4, everything())

# temp <- df %>% dplyr::select(sample, rownames, genotype, cohort_pheno)
# 
# temp <- temp %>% filter(rownames == "chr21:10485736_T/C")
# temp <- temp %>% unique()
# 
# # Create a subset of the data with only cases (cohort_pheno == 1)
# cases <- temp %>%
# 	filter(cohort_pheno == "1")
# 
# # Create a subset of the data with only controls (cohort_pheno == 0)
# controls <- temp %>%
# 	filter(cohort_pheno == "0")
# 
# # Define a function to perform the test for a given variant
# test_variant <- function(variant_name) {
# 	# Calculate the contingency table for the variant
# 	table_var <- table(cases$genotype[cases$rownames == variant_name],
# 							 controls$genotype[controls$rownames == variant_name])
# 	
# 	# Perform a Fisher's exact test for the difference between cases and controls
# 	fisher.test(table_var)$p.value
# }
# 
# # Apply the test_variant function to all variants and store the p-values in a list
# p_values <- lapply(unique(df$rownames), test_variant)
# 
# # Convert the list of p-values to a data frame and add the variant names as a column
# results <- data.frame(variant_name = unique(df$rownames),
# 							 p_value = unlist(p_values))

# Notes ----
# CADD_PHRED >30 likely deleterious. Variants with scores over 30 are predicted to be the 0.1% most deleterious possible substitutions in the human genome. We strongly recommend the actual score is used when assessing a variant and a cut-off appropriate to your requirements is chosen.

# REVEL  It integrates scores from MutPred, FATHMM v2.3, VEST 3.0, PolyPhen-2, SIFT, PROVEAN, MutationAssessor, MutationTaster, LRT, GERP++, SiPhy, phyloP, and phastCons. Score range from 0 to 1 and variants with higher scores are predicted to be more likely to be pathogenic.
# REVEL does not provide a descriptive prediction but for convenience, we display scores above 0.5, as 'likely disease causing' and display scores below 0.5 as 'likely benign'. REVEL_rankscore, REVEL_score

# MetaLR uses logistic regression to integrate nine independent variant deleteriousness scores and allele frequency information to predict the deleteriousness of missense variants. Variants are classified as 'tolerated' or 'damaging'; a score between 0 and 1 is also provided and variants with higher scores are more likely to be deleterious.

# MutationAssessor predicts the functional impact of amino-acid substitutions in proteins using the evolutionary conservation of the affected amino acid in protein homologs. We display the prediction, which is one of 'neutral', 'low', 'medium' and 'high', and the rank score, which is between 0 and 1 where variants with higher scores are more likely to be deleterious. 

# PolyPhen and SIFT results are heavily dependent on sequence conservation estimates derived from protein sequence alignments and using different versions of the protein databases can result in substantial variance in the predictions and scores obtained.
# Polyphen greater than 0.908	"Probably Damaging"
# SIFT a score < 0.05 are called 'deleterious' and all others are called 'tolerated'.

# GERP conservation scores as computed with the Genomic Evolutionary Rate Profiling GERP software on Multiple Sequence Alignments of whole-genomes. GERP identifies constrained loci in multiple sequence alignments by comparing the level of substitution observed to that expected if there was no functional constraint. Positive scores represent highly-conserved positions while negative scores represent highly-variable positions. the highest score of any base in a multi-base deletion is displayed. the mean of the scores of the two flanking bases is shown for an insertion

# GERP.._NR
# GERP.._RS_rankscore
# GERP.._RS