## ggraph ----
library(ggplot2)
library(ggraph)
library(dplyr)
library(tidyr)
library(igraph)

# We developed \texttt{untanglR} to construct protein pathway network graphs from the outputs of \texttt{ProteoMCLustR} and \texttt{Archipelago} using the igraph package.
# The size of the nodes is determined by the degree of the node (the number of edges connected to it) divided by the graph density, defined as the ratio of the actual number of edges to the maximum possible number of edges in the graph. This method provides an intuitive size scaling that reflects the relative connectivity of each node within the context of the overall network structure.
# The color of the nodes (referred to as 'degree_colors' in the code) is assigned based on the degree of each node, ranging from red (for nodes with the fewest connections) to yellow (for nodes with the most connections).
# Edge weights and layout methods also influence the visual representation of the graph. The displayed plots employ the 'stress' layout method, which utilizes multidimensional scaling (MDS) to determine the positions of the nodes. This method seeks to place nodes such that the distance between them in the layout reflects the 'distance' between them in the network (often, the number of edges in the shortest path between them).
# In protein-protein interaction (PPI) networks, the degree of a protein node corresponds to the number of interactions that protein has with other proteins within the network. Proteins with higher degrees often play pivotal roles in cellular functions, indicating their essentiality. However, proteins with a small number of high-confidence interactions can also be vital, especially in the context of monogenic rare diseases where the implicated protein may interact selectively with only a few partners.


# NOT USED:
# The Fruchterman-Reingold layout is a force-directed algorithm for node placement in network visualization. It simulates a physical system to layout the network: nodes are treated as atomic particles that repel each other, and edges are treated as springs holding adjacent nodes together. The algorithm iteratively adjusts node positions to minimize the total system energy, producing aesthetically pleasing layouts for many types of graphs.
subset_MCL_ID <- 22
subset_MCL_ID <- 586
subset_MCL_ID <- c(22, 586)

# Read enriched pathway
mcl_enriched <-
  data.table::fread('../../data/post_ppi/skat_log_psig_csv.txt') 

names(mcl_enriched)
# Separate rows
mcl_enriched_long <- mcl_enriched  |>
  separate_rows(Genes_clean_csv, sep = ",") |>
  rename(Gene_ID = Genes_clean_csv)

# Check the result
head(mcl_enriched_long)


# # Define function to get a subnetwork from STRING database for a given set of protein IDs
# GetSubNetwork <- function(string_db, STRING_IDs) {
#   # Get the subnetwork and simplify the resulting graph
#   string_subgraph <- string_db$get_subnetwork(STRING_IDs)
#   string_subgraph <- igraph::simplify(string_subgraph)
# }

# Network plot ----
# <https://bioinformatics.stackexchange.com/questions/282/how-to-manipulate-protein-interaction-network-from-string-database-in-r>

# Load STRING data from file
string_db <- readRDS(file = "../ProteoMCLustR_github_clone/data/ppi/string_data_700.rds")

# Read protein info file into data frame and select relevant columns
string_id_df <-
  data.table::fread('../ProteoMCLustR_github_clone/data/ppi/9606.protein.info.v11.5.txt.gz') |>
  dplyr::select(#STRING_ID = protein_external_id,
    STRING_ID = "#string_protein_id", Gene_ID = preferred_name)

# get gene list as a vector
selected_genes <- mcl_enriched_long |> filter(MCL_ID == subset_MCL_ID) |> dplyr::pull(Gene_ID)


selected_genes <- mcl_enriched_long |> 
dplyr::filter(MCL_ID %in% subset_MCL_ID) |>
dplyr::pull(Gene_ID) 

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
# V(string_igraph)$name <- protein_names[match(V(string_igraph)$name, protein_names)] # ENSP....

# gets the symbol as used in plot, etc.
V(string_igraph)$gene_symbol <- string_id_df$Gene_ID[match(V(string_igraph)$name, string_id_df$STRING_ID)] # gene name

# Set edge weights using the combined score from STRINGdb
E(string_igraph)$weight <- selected_interactions$combined_score

# Compute the Fruchterman-Reingold layout with scaling factor of 10
# layout <- layout_nicely(string_igraph)

# scale ggraph ----
g <- string_igraph
V(g)$degree <- degree(g, mode = "in")

# this function is borrowed from the ambient package
normalise <- function (x, from = range(x), to = c(0, 1)) {
  x <- (x - from[1])/(from[2] - from[1])
  if (!identical(to, c(0, 1))) {
    x <- x * (to[2] - to[1]) + to[1]
  }
  x
}

V(g)$degree<- normalise(V(g)$degree, to = c(3, 100))

# weights <- normalise(E(string_igraph)$weight, to = c(.01, 1))
weights <- (E(string_igraph)$weight)

# set.seed(123)
# layout_nicely10 <- layout_nicely(g) * 10

# make the size inversely proportional to the graph density
graph_density <- ecount(g) / ((vcount(g)*(vcount(g)-1))/2)

# use graph_density to adjust text size, using some function of graph_density
text_size <- 1 / graph_density 

# Compute the Fruchterman-Reingold layout
layout <- layout_with_fr(string_igraph)

# Use the computed layout in the ggraph function
# ggraph(g, layout = layout) +
p_net <- ggraph(g, "stress") +
  geom_edge_link(aes(width = weights)) +
  scale_edge_width(range = c(.1, .4), guide = "none") +
  geom_node_point(aes(size = I((degree/graph_density)/20)*1.1, alpha = I(degree)/10), color = "black") +
  geom_node_point(aes(size = I((degree/graph_density)/20), alpha = I(degree)/10), color = V(g)$group)+ 
  # geom_node_text(aes(label = name), size=2, color="gray50", repel=T) + # ENSP...
  geom_node_text(aes(label = gene_symbol), size=text_size, color="gray20", repel=T) + 
  scale_color_manual(values = V(string_igraph)$group) +
  theme_graph()


filename <- paste0("../../data/untangleR/ppi_network_", paste(subset_MCL_ID, collapse="_"), ".pdf")
ggsave(filename, 
       plot = p_net
       #, height = 10, width = 10
       )



# # Plot the graph with the new group attribute
# plot(string_igraph, 
#      layout = layout.fruchterman.reingold(string_igraph), 
#      vertex.shape = "circle", 
#      edge.arrow.size = 0.1,
#      vertex.size = 5, 
#      vertex.label = NA, 
#      vertex.color = V(string_igraph)$group, 
#      edge.color = "blue", 
#      edge.width = (E(string_igraph)$weight)/1000)
# 
# library(networkD3)
# 
# # Create a color vector based on the degree of each vertex
# degree_colors <- colorRampPalette(c("blue", "red"))(vcount(string_igraph))
# 
# # Assign the group attribute to the vertices
# V(string_igraph)$group <- degree_colors
# 
# 
# # End ----
# 
# sessionInfo()
# 
# sink()