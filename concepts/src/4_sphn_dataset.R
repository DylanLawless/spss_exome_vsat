
df_sphn <- readxl::read_excel(path = "../data/SPHN_dataset_release_2024_2_20240502.xlsx",
                           sheet="Concepts") # what kind of idiot would release this as excel sheets

df_sub <- df_sphn |> dplyr::filter(`concept reference` == "Single Nucleotide Variation")

names(df_sub)

df_sub <- 
  df_sub |> 
  dplyr::select("concept reference",
                "general concept name",
                "cardinality for composedOf",
                "concept or concept compositions or inherited",
                "type",
                "parent",
                everything())
