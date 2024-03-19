# GAP (Genotype-phenotype Association Predictor)

GAP is a R implementation of the research by Islam, UI and Assis, R : A Framework for learning associations between genes and phenotypes from multi-speices sequence alignments.

This package provides users the following functionalites: 1. Find Neural Network architectures for a given gene (GULO), that accounts for Cross validation error of zero where the network maps the multi-species exons of the gene to a specific binary trait. If a specific architecture can be mapp between the gene and trait perfectly, the positions within the exons of the gene that have greater significance is highlighted by dispalying a manhattan plot based on benjamini adjusted p-values. It runs Batch command with a list of genes to map them with the trait of interest based on the architecture that performs well from the previous function. It can also find all the genes that can be associated to the trait from the result of previous analysis. `/data-raw` directory contain testing dataset where the range of gene number is 21001-21050.

# Required Packages

Functionality of GAP is dependant on these prebuilt R packages. R packages can be installed from the R console with the command: `install.packages("package_name")`.

1.  "ann2"
2.  "parallel" for second function only
3.  "dplyr"
4.  "gtools"
5.  "qqman"

# Tree Feature Extraction from Newick Format Tree
This Python script extracts tree features from a Newick format tree using the ETE Toolkit (ete3). It generates feature vectors representing species presence/absence within different branches of the tree.

**Requirements**

    Python 3
    ete3
    pandas

**Usage**

    Function: generate_features(tree_str)
        Input: tree_str (Newick format tree string)
        Output: Feature vectors for each branch in the tree

    Example Input:
        "((A,(B,C),D),E)"

**Output**

The code outputs feature vectors indicating species presence/absence within tree branches.

**Notes**

Ensure the Newick tree string format is accurate.
Species names should be represented as comma-separated alphanumeric strings in the Newick tree string.

# GAP methods for performing model training.

## Architecture with minimum CV error
This command-line interface (CLI) command executes the `runGeneArchi.R` script. It tries different architectures of neural networks based on the input parameteres to map multi-species genetic data to binary phenotypes. Architectures with minimum CV error are displayed alonside the manhattan plot containing nucleotide positional importance for exons within the gene(s).

**Command Structure**

```bash
Rscript <path_to_runGeneArchi.R> <path_source> boolean_tree_flag num_hidden_layers
```
**Arguments**:
  - `path_source`: Source of the downloaded tool, i.e. `/users/username/downloads/_GAP/`.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `num_hidden_layers`: Hidden layer count {0,1,2,3}.
    
**Details**:
  - For `num_hidden_layers` > 0, explores permutations (1,2,3,..., number_of_species).
  - Displays a Manhattan plot for qualifying architectures.
    
**Sample Command**:
```bash 
Rscript.exe "/users/username/downloads/_GAP/runGeneArchi.R" "/users/username/downloads/_GAP/" TRUE 2.
```
which runs the script to identify NN architectures with minimum CV error with tree-features included dataset for all 2-hidden layered architectures.


## Run selected architecture on a range of genes

This command-line interface (CLI) command executes the `run_associated.R` script. It trains neural network architectures with satsifactory CV scores found from previous function on a specified gene list against multi-species genetic data to map binary phenotypes. The results are then stored in the designated folder `results`, please create a directory named results under the GAP directory, i.e., `/users/username/downloads/_GAP/results` where the results will be stored after processing.

**Command Structure**

```bash
Rscript.exe <path_to_run_associated.R> <path_source> boolean_tree_flag num_hidden_layers starting_gene ending_gene
```
This method requires execution from the shell due to resource constraints and follows a specific command structure with parallel processing:

**Arguments**:
  - `path_source`: Source of the downloaded tool, i.e. `/users/username/downloads/_GAP/`.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `num_hidden_layers`: Hidden layer count {0,1,2,3}, use 0 <-no hidden layers, integer(n)<-single hidden layer with n nodes, string(a,b,..,e)<- 5 hidden layer architecture with each layer having nodes corresponding to its positions, 1st layer<-a nodes, 2nd layer<-b nodes, 5th layer<-e nodes etc.
  - `starting_gene`: specifies the gene number from where processing will start.
  - `ending_gene`: species the ending gene number of processing.

**Sample Command**:
```bash 
Rscript.exe "/users/username/downloads/_GAP/run_associated.R" "/users/username/downloads/_GAP/" TRUE 1 21030 21040.
```
where the tool is located under the directory `/users/username/downloads/_GAP/` and this command will execute all the single-hidden layered architectures for the specified gene list and store the results under the specified directory.

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


