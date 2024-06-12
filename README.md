# SPSS vsat

Rare variants in infection response protein pathway associated with sepsis in children.

<!-- ![](images/logo.webp) -->
<div style="display: flex; justify-content: space-between;">
  <img src="images/logo.webp" style="width: 50%;" alt="Logo of DNA wizard"/>
</div>

# Introduction
Welcome to the repository for our project on analyzing rare variants in infection response protein pathways associated with sepsis in children. 
This repository contains all the code, data, and resources used in "Rare variants in infection response protein pathway associated with sepsis in children".

Our project explores the genetic underpinnings of immune response in pediatric sepsis, employing a combination of novel computational tools and in-depth statistical analysis. 
These tools include ProteoMCLustR for protein pathway clustering, SkatRbrain and Archipelago for advanced statistical analyses, and ACMGuru, untangleR, and AutoDestructR for clinical genetics interpretation.

In this repository, you will find structured directories containing R scripts and source code that demonstrate our workflow from data preparation through to variant analysis and interpretation. 
Each tool and method used in our study is documented in detail to aid reproducibility and further research in this critical area of pediatric healthcare.

Contributions and insights from readers are highly encouraged, as we aim to foster a collaborative environment to enhance our understanding of sepsis at the genetic level.

# Script order
* **Interpret**
* 1. AMCGuru_singlecase
    - requires: stand_alone_vcf_to_table
* AMCGuru_post_ppi
* archipelag
* untangleR

* **Data prep**
* joint_pca
* ppi
* ProteoMCLustR_github_clone
 
* **VSAT**
* variant_level/
* gene_level/
 


# 1.Single variant analysis
## Code Structure and Documentation
This repository houses scripts crucial for interpreting variant data to determine pathogenicity in genetic studies. The directory structure and contents are designed to facilitate the analysis of chromosome-split gVCF files processed with GATK and annotated with Ensembl VEP alongside a set of annotation databases.

### Directory Structure
```
.
├── ACMGuru_singlecase_vcurrent.R      # Main script for ACMG-based variant interpretation
├── directory_structure.txt            # Documentation of the directory structure
├── stand_alone_vcf_to_table           # Directory containing scripts to process VCF files
│   ├── gather.R                       # Script to gather data from processed VCF
│   ├── genotype_clean.R               # Script to clean and format genotype data
│   ├── progress_bar.R                 # Displays progress during data processing
│   ├── stand_alone_vcf_to_table.R     # Main script to convert VCF to table format
│   └── vcf_to_tables.R                # Processes VCF files to extract variant data
└── sync.sh                            # Bash script to synchronize essential data files
```

### Script Descriptions
- **ACMGuru_singlecase_vcurrent.R**: Implements ACMG guidelines to interpret variants for single-case analysis. It utilizes multiple libraries such as `dplyr`, `ggplot2`, and `stringr` to process genetic data, annotate it with clinical significance, and visually represent the findings.
- **gather.R** and **vcf_to_tables.R**: These scripts are part of a pipeline within the `stand_alone_vcf_to_table` directory that processes VCF files, extracting relevant genetic information and structuring it into a usable format for further analysis.
- **genotype_clean.R**: Cleans and normalizes genotype information to ensure consistency across datasets, which is critical for accurate variant interpretation.
- **progress_bar.R**: Provides a visual indicator of progress when processing large genetic datasets, enhancing user experience during long-running operations.
- **stand_alone_vcf_to_table.R**: Coordinates the conversion of VCF files from initial read to final table output, ensuring all components function seamlessly.
- **sync.sh**: A utility script to synchronize or update genetic data files from remote or local sources to the current project directory, ensuring data consistency.

### Inputs
- **Chromosome-split gVCF Files**: These are the primary inputs, processed using GATK and annotated with Ensembl VEP and other annotation databases to ensure comprehensive genetic data analysis.

### Outputs
- **Results Tables**: Detailed tables containing interpreted genetic variants along with their clinical significance.
- **Visualisations**: Graphical representations of the analysis, providing insights into the pathogenicity of the variants.

### Usage
- To perform a full variant analysis, run the `ACMGuru_singlecase_vcurrent.R` script, which orchestrates the variant interpretation using ACMG guidelines, outputting potential pathogenic variants and detailed visualizations.
- Use `sync.sh` to ensure all necessary data files are present and up to date before starting your analysis.

