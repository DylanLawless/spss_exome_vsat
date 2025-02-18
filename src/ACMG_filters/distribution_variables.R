dist_variables <- df
# Define the columns used for each criterion
columns_criteria <- list(
  PVS1 = c("IMPACT", "genotype", "Inheritance"),
  PS1 = c("CLIN_SIG", "ClinVar_CLNDN.y"),
  PS3 = c("genotype", "Inheritance"),
  PS5 = c("IMPACT", "sample", "SYMBOL"),
  PM2 = c("gnomAD_AF")
)

# Function to check if columns exist in the dataframe
check_columns_exist <- function(dist_variables, columns_list) {
  sapply(columns_list, function(columns) {
    all(columns %in% names(dist_variables))
  })
}

# Apply the function to check for column existence
column_existence <- check_columns_exist(dist_variables, columns_criteria)

# Print the results
print(column_existence)



# Function to check the class of each column
check_column_classes <- function(dist_variables, columns_list) {
  lapply(columns_list, function(columns) {
    sapply(columns, function(column) {
      if(column %in% names(dist_variables)) {
        class(dist_variables[[column]])
      } else {
        NA  # Return NA if column does not exist
      }
    })
  })
}

# Apply the function to get column classes
column_classes <- check_column_classes(dist_variables, columns_criteria)
print(column_classes)

# Load necessary library for plotting
library(ggplot2)

# Histogram for 'genotype' (used in both PVS1 and PS3)
if ("genotype" %in% names(dist_variables)) {
p1<-  ggplot(dist_variables, aes(x = genotype)) +
    geom_histogram(binwidth = 1, fill = "#fbd693", color = "black") +
    ggtitle("Genotype") +
    xlab("Genotype") +
    ylab("Frequency") +
    theme_bw()
}

# Histogram for 'gnomAD_AF' (used in PM2)
if ("gnomAD_AF" %in% names(dist_variables)) {
p2<-   ggplot(dist_variables, aes(x = gnomAD_AF)) +
    geom_histogram(fill = "#fbc393", color = "black", bins = 30) +
    ggtitle("gnomAD_AF") +
    xlab("gnomAD_AF") +
    ylab("Frequency") +
    theme_bw()
}

# Bar plot for 'IMPACT' (used in PVS1 and PS5)
if ("IMPACT" %in% names(dist_variables)) {
p3<-   ggplot(dist_variables, aes(x = IMPACT)) +
    geom_bar(fill = "#fbad85", color = "black") +
    ggtitle("IMPACT") +
    xlab("IMPACT") +
    ylab("Count") +
    theme_bw()
}

# Bar plot for 'Inheritance' (used in PVS1 and PS3)
if ("Inheritance" %in% names(dist_variables)) {
  p4<-   ggplot(dist_variables, aes(x = Inheritance)) +
    geom_bar(fill = "#ff9981", color = "black") +
    ggtitle("Inheritance") +
    xlab("Inheritance") +
    ylab("Count") +
    theme_bw()
}

# Bar plot for 'CLIN_SIG' (used in PS1)
if ("CLIN_SIG" %in% names(dist_variables)) {
  dist_variables$CLIN_SIG_trimmed <- substr(dist_variables$CLIN_SIG, 1, 20)
  p5<- ggplot(dist_variables, aes(x = CLIN_SIG_trimmed)) +
    geom_bar(fill = "#e08778", color = "black") +
    ggtitle("CLIN_SIG") +
    xlab("CLIN_SIG") +
    ylab("Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# Bar plot for 'ClinVar_CLNDN.y' (used in PS1)

if ("ClinVar_CLNDN.y" %in% names(dist_variables)) {
  dist_variables$ClinVar_CLNDN_trimmed <- substr(dist_variables$ClinVar_CLNDN.y, 1, 20)
  labels <- unique(dist_variables$ClinVar_CLNDN_trimmed)
  labels <- sort(labels)
  
  # Calculate indices to select labels evenly spaced
  if (length(labels) > 10) {
    indices <- round(seq(1, length(labels), length.out = 10))
    selected_labels <- labels[indices]
  } else {
    selected_labels <- labels
  }
  
  p6 <- ggplot(dist_variables, aes(x = ClinVar_CLNDN_trimmed)) +
    geom_bar(fill = "#b86572", color = "black") +
    ggtitle("ClinVar_CLNDN (n=10 samp)") +
    xlab("ClinVar_CLNDN") +
    ylab("Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_discrete(breaks = selected_labels, labels = selected_labels)  # Use the selected labels
}

# Check if 'SYMBOL' is present in the dataframe
if ("SYMBOL" %in% names(dist_variables)) {
  # Get unique labels and sort them alphabetically
  symbols <- unique(dist_variables$SYMBOL)
  symbols <- sort(symbols)
  
  # Calculate indices to select labels evenly spaced
  if (length(symbols) > 10) {
    # Generate indices for evenly spaced labels
    indices <- round(seq(1, length(symbols), length.out = 10))
    selected_symbols <- symbols[indices]
  } else {
    selected_symbols <- symbols
  }
  
  # Create the plot using the selected labels
  p7 <- ggplot(dist_variables, aes(x = SYMBOL)) +
    geom_bar(fill = "#575960", color = "black") +
    ggtitle("SYMBOL (n=10 samp)") +
    xlab("SYMBOL") +
    ylab("Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_discrete(breaks = selected_symbols, labels = selected_symbols)  # Use the selected labels
}


library(patchwork)

dist_patch <- (p1 + p2 + p3 + p4) / (p5 + p6 + p7) 

# "#fbd693"
# "#fbc393"
# "#fbad85"
# "#ff9981"
# "#e08778"
# "#b86572"
# "#7e676b"

ggsave(paste("../../images/", output_directory, file_suffix, "distribution_variables.pdf", sep = "") ,plot = dist_patch, width = 10, height = 10)
rm(dist_patch, dist_variables)
