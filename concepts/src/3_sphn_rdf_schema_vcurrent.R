# Understanding RDF (Resource Description Framework) and OWL (Web Ontology Language) can indeed be complex. Here, I'll break down the terms and structures you've presented from snippet 3 regarding `sphn:SingleNucleotideVariation` to help clarify their meanings and how they function within RDF and OWL.
# 
# ### Overview of RDF and OWL Terms
# - **RDF (Resource Description Framework)**: A standard model for data interchange on the web, using triples (subject, predicate, object) to express data.
# - **OWL (Web Ontology Language)**: Built on top of RDF, OWL is used for defining structured, web-based ontologies that enable rich data integration and classification.
# 
# ### Breaking Down Snippet 3
# 
# #### Line by Line Explanation
# 1. **`owl:Class`**
# 	- This indicates that `sphn:SingleNucleotideVariation` is defined as a class in the ontology.
# 
# 2. **`rdfs:label "Single Nucleotide Variation"`**
# 	- A human-readable name for the class.
# 
# 3. **`rdfs:subClassOf`**
# 	- Indicates that `sphn:SingleNucleotideVariation` is a subclass of the complex class described in the following lines. This complex class is built using intersections of various restrictions.
# 
# 4. **`owl:intersectionOf`**
# 	- This is used to define a class by the intersection (logical AND) of several other classes. Here, it defines the necessary conditions for an individual to be considered an instance of `sphn:SingleNucleotideVariation`.
# 
# 5. **`owl:Restriction`**
# 	- A way to apply constraints on properties for individuals of a class. Each restriction is related to a specific property (e.g., `sphn:hasGenomicPosition`).
# 
# 6. **`owl:onProperty`**
# 	- Specifies the property that the restriction is applied to.
# 
# 7. **`owl:minCardinality "0"^^xsd:nonNegativeInteger`**
# 	- Specifies that the property (`sphn:hasGenomicPosition`, for instance) can be absent (minimum occurrence is 0).
# 
# 8. **`owl:maxCardinality "1"^^xsd:nonNegativeInteger`**
# 	- Indicates the property can occur at most once (maximum occurrence is 1).
# 
# 9. **`owl:someValuesFrom`**
# 	- Specifies that the property must have at least one value from the specified class (e.g., `sphn:GenomicPosition`).
# 
# #### Understanding the Restrictions
# - The `owl:Restriction` elements define what properties an instance of `sphn:SingleNucleotideVariation` must have, and what the characteristics of these properties are (like how many times they must appear and what types of values they should have).
# 
# #### Additional Terms
# - **`sphn:GeneticVariation`**
# 	- Indicates that `sphn:SingleNucleotideVariation` is related to the broader class `sphn:GeneticVariation`.
# 
# - **`owl:equivalentClass`**
# 	- Specifies that `sphn:SingleNucleotideVariation` is equivalent to another class identified by the term `so:0001483`. This creates an explicit semantic link between this class and a class defined possibly in another ontology.
# 
# - **`skos:definition`**
# 	- Provides a descriptive definition of the class, making the ontology more understandable and accessible.
# 
# ### Summary
# The RDF/OWL representation in snippet 3 provides a robust framework for defining the characteristics and constraints of the `sphn:SingleNucleotideVariation` class. It uses a combination of class definitions, property restrictions, and semantic linking to other ontology classes to establish a clear and enforceable schema for data that falls under this category.


# Load necessary library
library(rdflib)
library(dplyr)

# Parse the Turtle file
file_path <- "../data/sphn_rdf_schema.ttl"
rdf_data <- rdf_parse(file_path, format = "turtle")

# SPARQL query to retrieve related concepts and properties for GeneticVariation
sparql_query_1 <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-schema/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

SELECT ?property ?value
WHERE {
  sphn:SingleNucleotideVariation ?property ?value .
}
"

query_results <- rdf_data %>%
	rdf_query(sparql_query_1)

print(query_results)

sparql_query_2 <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-schema/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

SELECT ?property ?value ?superclass
WHERE {
  {
    sphn:SingleNucleotideVariation ?property ?value .
  } UNION {
    sphn:SingleNucleotideVariation rdfs:subClassOf ?superclass .
    ?superclass ?property ?value .
  }
}
"

query_results <- rdf_data %>%
	rdf_query(sparql_query_2)

print(query_results)

# split the URI property to make reading easier
df  <-query_results |> tidyr::separate(property, c("property_prefix", "property_fragment"), sep = "#", remove = F)

# 
# To proceed with analyzing the RDF ontology and building a requirement table for sphn:SingleNucleotideVariation, you've correctly identified the need to focus on properties ending in #intersectionOf. This step is crucial because owl:intersectionOf is often used in OWL ontologies to specify that a class is equivalent to the intersection (logical AND) of several other classes or conditions. These intersections typically define the necessary conditions (restrictions) that instances of a class must satisfy.
# 
# Clarifying the Terminology
# owl:intersectionOf: This OWL property is used to define a class as an intersection of several other classes or restrictions. In the context of RDF and OWL, it is used to express complex class definitions that require an instance to meet all specified conditions.
# Next Step: Querying and Analyzing owl:intersectionOf
# To refine your data extraction and focus on entries related to owl:intersectionOf, you can modify your SPARQL query to specifically look for this property and extract relevant information about the conditions that form part of the intersection. Here's how you could adjust the query:



# Run the modified query to focus on intersectionOf
sparql_query <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-schema/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

SELECT ?property ?value ?superclass
WHERE {
  {
    sphn:SingleNucleotideVariation ?property ?value .
    FILTER regex(str(?property), \"intersectionOf$\")
  } UNION {
    sphn:SingleNucleotideVariation rdfs:subClassOf ?superclass .
    ?superclass ?property ?value .
    FILTER regex(str(?property), \"intersectionOf$\")
  }
}
"

query_results <- rdf_data %>% rdf_query(sparql_query)

# Process the results
df <- query_results %>%
	mutate(property = as.character(property), value = as.character(value)) %>%
	tidyr::separate(property, c("property_prefix", "property_fragment"), sep = "#", remove = F) %>%
	filter(property_fragment == "intersectionOf")

# Print the processed data
print(df)

# Now that you've successfully extracted the intersectionOf parts that define the constraints for sphn:SingleNucleotideVariation, the next step is to delve deeper into these parts to understand what each intersection entails. This involves retrieving details about the conditions (restrictions) specified by these intersections, such as specific properties (owl:onProperty), cardinalities (owl:minCardinality, owl:maxCardinality), and the types of values allowed (owl:someValuesFrom).

library(rdflib)

# Assuming rdf_data is already loaded
sparql_query_trace <- "
PREFIX sphn: <https://biomedit.ch/rdf/sphn-schema/sphn#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

SELECT ?subject ?property ?object
WHERE {
  sphn:SingleNucleotideVariation ?property ?subject .
  ?subject ?prop ?object .
  FILTER(isBlank(?subject))
}
"

# Execute the SPARQL query
trace_query_results <- rdf_data %>%
	rdf_query(sparql_query_trace)

# Print the results
print(trace_query_results)


# Now that you have successfully identified the intersectionOf nodes related to sphn:SingleNucleotideVariation, the next step is to explore these nodes in more detail to extract all the necessary attributes and restrictions for synthetic data generation.
# 
# Goals for Detailed Exploration:
# 	Identify all properties and constraints linked to each intersectionOf node:
# 	What properties (owl:onProperty) are associated with these nodes?
# 	What are the cardinality constraints (owl:minCardinality, owl:maxCardinality)?
# 	Are there any datatype or class constraints (owl:someValuesFrom, owl:allValuesFrom)?


# new approach ----
library(rdflib)
library(dplyr)
library(tidyr)


# Load the RDF data from the Turtle file
file_path <- "../data/sphn_rdf_schema.ttl"
rdf_data <- rdf_parse(file_path, format = "turtle")

# SPARQL query to retrieve all triples
sparql_query <- "
SELECT ?subject ?predicate ?object
WHERE {
  ?subject ?predicate ?object .
}
"

# Execute the SPARQL query
query_results <- rdf_data %>%
	rdf_query(sparql_query)

# Convert results to a dataframe
rdf_df <- as.data.frame(query_results, stringsAsFactors = FALSE)

df_subclass <- rdf_df %>%
	filter(grepl("#GeneticVariation", object)) %>%
	filter(grepl("#subClassOf", predicate))

# we see that the first subject contins "#SingleNucleotideVariation".
# we now want to find ALL information which would relate to this concept since we need to make a datapoint that will fully adhere to the concept requirements. lets search the entire dataset for all mentions. 

#SingleNucleotideVariation
# Assuming rdf_df is already loaded and contains columns named subject, predicate, and object
df_snv <- rdf_df %>%
	filter(grepl("SingleNucleotideVariation", subject, ignore.case = TRUE) |
			 	grepl("SingleNucleotideVariation", predicate, ignore.case = TRUE) |
			 	grepl("SingleNucleotideVariation", object, ignore.case = TRUE))

# we now have a dataset which points to other entries. should we extend the search into suject or object from here?

# related_subjects <- unique(c(df_snv$subject, df_snv$object))
df_nodes <- df_snv |> filter(!grepl("SingleNucleotideVariation", subject, ignore.case = TRUE))

related_subjects <- unique(c(df_snv$subject))
related_subjects <- unique(c(df_nodes$subject))
related_subjects

# Extend the search to include these related subjects and objects
df_extended <- rdf_df %>% filter(subject %in% related_subjects | object %in% related_subjects)

# Print extended results
print(df_extended)
