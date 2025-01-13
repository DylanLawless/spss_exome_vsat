## ggraph ----
library(ggplot2)
library(ggraph)
library(dplyr)
library(tidyr)
library(igraph)
set.seed(666)
# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("STRINGdb")
# <https://bioinformatics.stackexchange.com/questions/282/how-to-manipulate-protein-interaction-network-from-string-database-in-r>

# Import and prep ----
df_report_sample_vsat <- readRDS(file = "../../data/ACMGuru_post_ppi/df_report_sample_vsat_22_586_836.Rds")

# Read enriched pathway
mcl_enriched <-
  data.table::fread('../../data/post_ppi/skat_log_psig_csv.txt') 

mcl_enriched_long <- mcl_enriched  |>
  separate_rows(Genes_clean_csv, sep = ",") |>
  rename(Gene_ID = Genes_clean_csv)

head(mcl_enriched_long)

# Load STRING data from file
string_db <- readRDS(file = "../ProteoMCLustR_github_clone/data/ppi/string_data_700.rds")

# Read protein info file into data frame and select relevant columns
string_id_df <-
  data.table::fread('../ProteoMCLustR_github_clone/data/ppi/9606.protein.info.v11.5.txt.gz') |>
  dplyr::select(#STRING_ID = protein_external_id,
    STRING_ID = "#string_protein_id", Gene_ID = preferred_name)

# function plot ----
plot_network_for_mcl_id <- function(current_subset_MCL_ID, mcl_enriched_long, string_db, string_id_df) {
  
  # Network plot ----

  # get gene list as a vector
  selected_genes <- mcl_enriched_long |> 
    dplyr::filter(MCL_ID %in% current_subset_MCL_ID) |>
    dplyr::pull(Gene_ID) 
  
  # Create a label string from current_subset_MCL_ID
  label_id <- paste(current_subset_MCL_ID, collapse = ", ")
  
  # subset ----
  # Get the STRINGdb subset for gene list
  string_id_df <- string_id_df %>%
    dplyr::filter(Gene_ID %in% selected_genes)
  
  # Select the relevant interactions (edges) from the network
  selected_interactions <- string_db$get_interactions(string_id_df$STRING_ID)
  
  # Create an igraph object from the selected edges
  string_igraph <- graph_from_data_frame(selected_interactions[, c("from", "to")])
  
  # Assign the group attribute to the vertices
  V(string_igraph)$group <- "groupX"
  
  # Create a color vector based on the degree of each vertex
  degree_colors <- colorRampPalette(c("red", "yellow"))(vcount(string_igraph))
  
  # Assign the group attribute to the vertices
  V(string_igraph)$group <- degree_colors
  
  # Set vertex attributes (protein names)
  protein_names <- unique(c(selected_interactions$from, 
                            selected_interactions$to))
  
  # gets the ID as used in plot, etc.
  V(string_igraph)$name <- protein_names[match(V(string_igraph)$name, protein_names)] # ENSP....
  
  # gets the symbol as used in plot, etc.
  V(string_igraph)$gene_symbol <- string_id_df$Gene_ID[match(V(string_igraph)$name, string_id_df$STRING_ID)] # gene name
  
  # Set edge weights using the combined score from STRINGdb
  E(string_igraph)$weight <- selected_interactions$combined_score
  rm(selected_interactions)

  # scale ggraph ----
  g <- string_igraph
  
  # weights <- normalise(E(string_igraph)$weight, to = c(.01, 1))
  weights <- (E(string_igraph)$weight)
  
  # Calculate degrees
  V(g)$degree <- degree(g, mode = "all")  # 'all' might be more appropriate depending on your graph's directedness
  
  # Normalize function
  normalise <- function (x, to = c(3, 15)) {
    if(max(x) == min(x)) return(rep(1, length(x)))  # avoid division by zero if all degrees are the same
    (x - min(x)) / (max(x) - min(x)) * (to[2] - to[1]) + to[1]
  }
  
  # Apply normalization to vertex degrees for plotting
  print("degree")
  print(V(g)$degree |> head())
  V(g)$size <- normalise(V(g)$degree)
  print("degree norm")
  print(V(g)$size |> head())
  
  # make the size inversely proportional to the graph density
  print(ecount(g))
  print(vcount(g))
  graph_density <- ecount(g) / ((vcount(g)*(vcount(g)-1))/2)
  
  # use graph_density to adjust text size, using some function of graph_density
  print(paste("density: ", graph_density))
  # text_size <- 1 / graph_density
  # text_size <- 7 * graph_density
  text_size <- (1/(4 - graph_density)) * 24
  
  
  print(paste("text size: ", text_size))
  
  # Compute layout
  # layout <- layout_nicely(string_igraph)
  # layout_nicely10 <- layout_nicely(g) * 10 # Fruchterman-Reingold layout with scaling factor of 10
  # Compute the Fruchterman-Reingold layout
  layout <- layout_with_fr(g)
  
  # Calculate degree for each node
  # print("degree string_igraph")
  # print(degree(g, mode = "all") |> head())
  # print("degree g")
  # print(degree(g, mode = "all")|> head())
  V(g)$degree <- degree(g, mode = "all")
  
  # Color mapping based on degree: Red to Yellow gradient
  max_degree <- max(V(g)$degree)
  min_degree <- min(V(g)$degree)
  V(g)$color_degree <- (V(g)$degree - min_degree) / (max_degree - min_degree)
  color_palette <- colorRampPalette(c("yellow", "red"))(100)
  V(g)$color <- color_palette[as.integer(round(V(g)$color_degree * 99)) + 1]
  
  # Use the computed layout in the ggraph function
  print(label_id)
  
  p_net <- ggraph(g, "stress") +
    geom_edge_link(aes(width = weights)) +
    scale_edge_width(range = c(.1, .4), guide = "none") +
    geom_node_point(aes(size = I(size)*1.1), alpha = .8, color = "black") +
    geom_node_point(aes(size = I(size)), alpha = .8 , color = V(g)$group)+
    # geom_node_text(aes(label = name), size=2, color="gray50", repel=T) + # ENSP...
    geom_node_text(aes(label = gene_symbol), size=text_size, color="gray20", repel=T) + 
    scale_color_manual(values = V(g)$group) +
    # theme_graph() +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(size = 20), 
      plot.background = element_rect(fill = "white", color = "white"),  # Sets plot background to white
      panel.background = element_rect(fill = "white", color = "white"),  # Sets panel background to white
      legend.background = element_rect(fill = "white", color = "white")  # Optional: sets legend background to white
    ) +
    labs(title = paste("Pathway ID:", label_id))
  
  print(p_net)
  filename <- paste0("../../images/untangleR/ppi_network_v2_", paste(current_subset_MCL_ID, collapse="_"), ".pdf")
  print(paste("saving: ", filename))
  ggsave(filename,
         plot = p_net
         , height = 5, width = 8
  )
}

plots <- list()

subset_MCL_ID <- c(22, 586, 836)  # IDs to subset on

# Call the function with all MCL_IDs collectively
joint_plot <- plot_network_for_mcl_id(subset_MCL_ID, mcl_enriched_long, string_db, string_id_df)
# print(joint_plot)

# Initialize a list to store plots for each individual MCL_ID
individual_plots <- list()

# Call the function for each MCL_ID individually
for (current_subset_MCL_ID in subset_MCL_ID) {
  individual_plot <- plot_network_for_mcl_id(list(current_subset_MCL_ID), mcl_enriched_long, string_db, string_id_df)
  individual_plots[[as.character(current_subset_MCL_ID)]] <- individual_plot
  print(individual_plot)
}
