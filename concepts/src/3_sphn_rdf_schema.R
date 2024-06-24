
# Correct file path (ensure it's the correct one and the file exists)
# file_path <- "~/Downloads/NGBO.ttl"  # Adjust as necessary
file_path <- "../data/sphn_rdf_schema.ttl"

# if (!requireNamespace("devtools", quietly = TRUE))
  # install.packages("devtools")
# devtools::install_github("ropensci/rdflib")
library(rdflib)

# Parse the Turtle file
rdf_data <- rdf_parse(file_path, format = "turtle")

# hasGene query ----
# sparql_query <- "
# PREFIX sphn: <https://biomedit.ch/rdf/sphn-ontology/sphn#>
# PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
# 
# SELECT ?label
# WHERE {
#   sphn:hasGene rdfs:label ?label .
# }
# "

# GeneticVariation query ----
sparql_query <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-ontology/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT ?label
WHERE {
  sphn:GeneticVariation rdfs:label ?label .
}
"

# Execute the SPARQL query
query_results <- rdf_data %>%
  rdf_query(sparql_query)

# Print the results
print(query_results)

# test ----

# Extended SPARQL query to fetch relationships of sphn:GeneticVariation
sparql_query <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-ontology/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>

SELECT ?subject ?predicate ?object
WHERE {
  { ?subject ?predicate sphn:GeneticVariation }
  UNION
  { sphn:Gene ?predicate ?object }
}
"

# Execute the SPARQL query
relationships <- rdf_data %>%
  rdf_query(sparql_query)

# Print the results
print(relationships)
dt <- tibble::as_tibble(relationships)

# graph ----

# Install and load igraph if not already
# if (!requireNamespace("igraph", quietly = TRUE))
  # install.packages("igraph")
library(igraph)

# Create a graph from the query results
graph <- graph_from_data_frame(relationships, directed = TRUE)

# Customizing plot
# plot(graph,
#      vertex.label = V(graph)$name,
#      vertex.size = 10,
#      vertex.color = "skyblue",
#      edge.arrow.size = 0.5,
#      layout = layout_nicely(graph))

# Install and load visNetwork if not already
# if (!requireNamespace("visNetwork", quietly = TRUE))
  # install.packages("visNetwork")
# simple label ----
library(igraph)

simplify_labels <- function(labels) {
  sapply(labels, function(label) {
    # Check if it is an anonymous node
    if (startsWith(label, "_:")) {
      return("Anonymous")
    }
    
    # Extract last meaningful part after the last delimiter
    parts <- unlist(strsplit(label, split = "[/#]"))
    if (length(parts) > 0) {
      return(tail(parts, n = 1))
    } else {
      return(label)
    }
  })
}

library(visNetwork)

# Prepare data for visNetwork
data <- toVisNetworkData(graph)
data$nodes$label <- simplify_labels(data$nodes$id)

# Generate interactive graph
visNetwork(data$nodes, data$edges, width = "100%") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)




# V2 ----
# Load necessary libraries
if (!requireNamespace("rdflib", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
  devtools::install_github("ropensci/rdflib")
}
library(rdflib)

if (!requireNamespace("igraph", quietly = TRUE))
  install.packages("igraph")
library(igraph)

if (!requireNamespace("visNetwork", quietly = TRUE))
  install.packages("visNetwork")
library(visNetwork)

# Define the file path to your Turtle file
file_path <- "../data/sphn_rdf_schema.ttl"

# Parse the Turtle file
rdf_data <- rdf_parse(file_path, format = "turtle")

# Define the SPARQL query to retrieve relationships
sparql_query <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-ontology/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>

SELECT ?subject ?predicate ?object
WHERE {
  { sphn:GeneticVariation ?predicate ?object }
  UNION
  { ?subject ?predicate sphn:GeneticVariation }
  UNION
  { ?subject rdfs:subClassOf sphn:GeneticVariation }
  UNION
  { ?subject rdfs:subClassOf ?mid .
    ?mid owl:onProperty ?predicate .
    ?predicate rdfs:range ?object }
}
"

# Execute the SPARQL query
query_results <- rdf_data %>%
  rdf_query(sparql_query)

# Convert the results to a tibble and print
relationships <- tibble::as_tibble(query_results)
print(relationships)
dt <- tibble::as_tibble(relationships)

# Function to ensure unique and simplified labels
simplify_labels <- function(labels) {
  sapply(seq_along(labels), function(i) {
    label <- labels[i]
    if (startsWith(label, "_:")) {
      paste("Anonymous", i, sep = "_")  # Append index to make it unique
    } else {
      simple_label <- sub('.*[/#](.+)$', '\\1', label)
      paste(simple_label, i, sep = "_")  # Append index to make it unique
    }
  })
}

# Create a graph from the query results
graph <- graph_from_data_frame(relationships, directed = TRUE)
V(graph)$name <- simplify_labels(V(graph)$name)  # Apply label simplification

# Prepare data for visNetwork
data <- toVisNetworkData(graph)
data$nodes$label <- simplify_labels(data$nodes$id)

# Ensure node labels are present and correct
print(data$nodes)

# Generate interactive graph with visNetwork
visNetwork(data$nodes, data$edges, width = "100%", height = "800px") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visEdges(arrows = 'to', labelHighlightBold = TRUE) %>%
  visNodes(size = 20, font = list(size = 20)) %>%
  visInteraction(hover = TRUE, tooltipDelay = 200) %>%
  visLayout(randomSeed = 42)

# V3 ----

library(rdflib)

# Load your RDF data
file_path <- "../data/sphn_rdf_schema.ttl"

# Parse the Turtle file
rdf_data <- rdf_parse(file_path, format = "turtle")

# Serialize into n-triples for a simpler, line-based format
# nt_out <- tempfile(fileext = ".nt")
# nt_out <-  "../data/sphn_rdf_schema.nt"
# rdf_serialize(rdf_data, doc = nt_out, format = "ntriples")

# Serialize RDF data into JSON-LD
jsonld_out <- tempfile(fileext = ".json")
rdf_serialize(rdf_data, doc = jsonld_out, format = "jsonld")

# Re-load the JSON-LD to ensure it's correctly formatted
rdf_jsonld <- rdf_parse(jsonld_out, format = "jsonld")
class(rdf_jsonld)

# If you have the JSON-LD data, convert it to a dataframe
library(jsonlite)
json_data <- fromJSON(jsonld_out)
df <- as.data.frame(json_data)

# Now you can use dplyr or similar tools to manipulate your data
library(dplyr)
filtered_data <- df %>%
  filter(predicate == "http://purl.org/dc/elements/1.1/creator") %>%
  select(subject, object)

print(filtered_data)




library(dplyr)
library(jsonlite)

# Assuming your JSON-LD data is loaded into a DataFrame
df <- fromJSON(jsonld_out)  # parse JSON-LD data into a DataFrame
df <- as.data.frame(df)

# Filter the DataFrame to find entries related to GeneticVariation
filtered <- df |>
  # filter(grepl("GeneticVariation", X.graph..id))
filter(grepl("GeneticVariation", df))
names(df)
         # | # grepl("hasGenomicPosition", predicate) 
         # | grepl("hasChromosomalLocation", predicate)) 


# Print the filtered data to see the results
print(filtered_data)

