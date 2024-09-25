
# Correct file path (ensure it's the correct one and the file exists)
# file_path <- "~/Downloads/NGBO.ttl"  # Adjust as necessary
file_path <- "../data/sphn_schema.ttl"

# if (!requireNamespace("devtools", quietly = TRUE))
  # install.packages("devtools")
# devtools::install_github("ropensci/rdflib")
library(rdflib)

# Parse the Turtle file
rdf_data <- rdf_parse(file_path, format = "turtle")

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
query_results <- rdf_data |>
  rdf_query(sparql_query)

# Print the results
print(query_results)

# try 1 ----
sparql_query <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-ontology/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

SELECT ?property ?label ?datatype ?minCardinality ?maxCardinality
WHERE {
  ?property rdfs:domain sphn:SingleNucleotideVariation ;
            rdfs:range ?datatype ;
            rdfs:label ?label .
  
  OPTIONAL { ?property owl:minCardinality ?minCardinality . }
  OPTIONAL { ?property owl:maxCardinality ?maxCardinality . }
}
"



sparql_query <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-ontology/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

SELECT DISTINCT ?property ?label ?range ?domain
WHERE {
  ?property rdfs:label ?label ;
            rdfs:range ?range .
  OPTIONAL { ?property rdfs:domain ?domain . }
  FILTER (EXISTS { ?property rdfs:domain/rdfs:subClassOf* sphn:SingleNucleotideVariation } OR
          NOT EXISTS { ?property rdfs:domain })
}
"

query_results <- rdf_data |>
  rdf_query(sparql_query)

print(query_results)



query_results <- rdf_data |>
  rdf_query(sparql_query)

print(query_results)
