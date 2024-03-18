# GAP (Genotype-phenotype Association Predictor)

GAP is a R implementation of the research by Islam, UI and Assis, R : A Framework for learning associations between genes and phenotypes from multi-speices sequence alignments.

This package provides users the following functionalites: 1. Find Neural Network architectures for a given gene (GULO), that accounts for Cross validation error of zero where the network maps the multi-species exons of the gene to a specific binary trait. If a specific architecture can be mapp between the gene and trait perfectly, the positions within the exons of the gene that have greater significance is highlighted by dispalying a manhattan plot based on benjamini adjusted p-values. It runs Batch command with a list of genes to map them with the trait of interest based on the architecture that performs well from the previous function. It can also find all the genes that can be associated to the trait from the result of previous analysis.

# Required Packages

Functionality of GAP is dependant on these prebuilt R packages. R packages can be installed from the R console with the command: `install.packages("package_name")`.

1.  "ann2"
2.  "parallel" for second function only
3.  "dplyr"
4.  "gtools"
5.  "qqman"

# Tree Feature Extraction from Newick Format Tree
This Python script extracts tree features from a Newick format tree using the ETE Toolkit (ete3). It generates feature vectors representing species presence/absence within different branches of the tree.
Requirements

    Python 3.x
    ETE Toolkit (ete3)
    Pandas (pandas)

Usage

    Function: generate_features(tree_str)
        Input: tree_str (Newick format tree string)
        Output: Feature vectors for each branch in the tree

    Example:
        Replace newick_tree variable with your tree string.
        Run the script to display feature vectors for each branch.

Output

The code outputs feature vectors indicating species presence/absence within tree branches.

Notes

Ensure the Newick tree string format is accurate.
Species names should be represented as comma-separated alphanumeric strings in the Newick tree string.

# GAP Functions: NN Architecture Analysis

## Overview

The GAP package provides tools for exploring gene architectures using neural networks to map binary traits and identify associated genes.

### Setup

1. **Working Directory**:
   - Ensure the R environment's working directory contains:
     - `get_gene_archi.r` script.
     - `data-raw` subdirectory with relevant example datasets.

2. **Function Loading**:
   - Load GAP package functions using:
     ```R
     source("getGeneArchi.R")
     source("run_associated.R")
     ```

## Methods

### `get_gene_architecture()`

- **Description**:
  - Identifies architectures mapping genes to binary traits with zero cross-validation error.
- **Arguments**:
  - `path_source`: Source of the downloaded tool, i.e. '/users/username/downloads/_GAP/'.
  - `use_tree_features`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `hidden_layers`: Hidden layer count {0,1,2,3}.
- **Details**:
  - For `hidden_layers` > 0, explores permutations (1,2,3,..., number_of_species).
  - Displays a Manhattan plot for qualifying architectures.

### `find_associated_genes()`

- **Description**:
  - Identifies gene IDs trained with zero cross-validation error in a previous analysis.
- **Arguments**:
  - `path_to_result`: Path storing previous analysis results.
  - `start_range`, `end_range`: File name range in raw data for association.
  - `path_genes`: Repository of raw data for known species for neural network model training.
  - `path_corresponding_genes`: Repository of raw data for status-unknown species.
- **Flags**:
  - `use_tree_features_flag`: Flag for tree features (as previously explained).
  - `path_tree`: Directory for tree features.

## Usage

- **Execution**:
  - Replace placeholders with appropriate values.
  - Run functions in the R environment with specified arguments.

## Get Associated Genes

This command-line interface (CLI) command facilitates the execution of the `run_associated.R` script. It trains neural networks using a specified gene list against multi-species genetic data to map binary phenotypes. The results are then stored in a designated folder.

### Command Structure

```bash
<path_to_Rscript.exe> <path_to_run_associated.R> <path_to_tree_features> <path_to_gene_list> <path_to_results_folder> <path_to_corresponding_genes> number_of_species hidden_layers use_tree_features_flag starting_gene ending_gene
```
### Method Overview

This method requires execution from the shell due to resource constraints and follows a specific command structure:

- **Usage Details**:
  - `Tree Features`: Obtained from species tree approximation through genetic data available at `~/data-raw/feature_from_dendrogram_species`.
  - `Gene List`: Users can provide their genomic data directory or use the provided gene list at `~/data-raw/`.
  - `Results Path`: Path to store the obtained results.

- **Command Structure**:
  - `number_of_species`: Represents the count of species in the dataset (e.g., 34 in our datasets).
  - `starting_gene` and `ending_gene`: Define the range of genes to be trained from the list.





Sample Command:
```bash 
"C:/program files/R/R-4.2.1/bin/Rscript.exe" "C:/Users/username/Documents/GAP/R/run_associated.R" "C:/Users/username/Documents/GAP/data-raw/feature_from_dendrogram_species" "C:/Users/username/Documents/GAP/data-raw/known_genes/" "C:/Users/username/Documents/GAP/data-raw/unknown_genes/" "D:/results/" 34 0 FALSE 21030 21040.
```
where the package is located under the directory `C:/users/username/documents/GAP`
