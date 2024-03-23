# GAP (Genotype-phenotype Association Predictor)

GAP is a software implementation of Islam, Campelo dos Santos, Kanjilal, and Assis' (2024) machine learning framework for predicting genotype-phenotype associations from multi-species sequence alignments.

We introduce GAP, a novel machine learning framework for predicting a binary phenotype from gaps in a multi-species sequence alignment of a genomic region. In particular, GAP employs a neural network to predict the presence or absence of a phenotype solely from gaps in a multiple alignment, with optional consideration of phylogenetic relationships from a user-input species tree. GAP can be employed for three distinct problems: predicting a phenotype in one or more species from a known associated region, pinpointing which specific genomic positions within such a region are associated with a phenotype, and extracting a set of candidate genomic regions associated with a phenotype. Here, we demonstrate the utility of GAP by exploiting the well-known association between the gulonolactone (L-) oxidase *Gulo* gene and vitamin C synthesis.


For improvment suggestions and queries, email: uislam2022@fau.edu

# Cite GAP

Thank you for using GAP! 

Citation is appreciated: UI Islam, AL Campelo dos Santos, R Kanjilal, and R Assis' (2024) machine learning framework for predicting genotype-phenotype associations from multi-species sequence alignments 

# Getting Started

Running GAP scripts installation of Python3 (version > 3.9) and R (version > 4.0) is recommended into the system. GAP installs the required R and Python packages from within the script. If GAP fails to install a package itself and users encounter a package installation error, they can install that particular package manually. 
Commands to install R and Python packages from respective terminals:
```bash
# R package installation command
install.packages("package_name")

# Python package installation command
pip3 install package_name
```

# GAP Input

GAP accepts multi-species sequence data as input in a specified format. In a GAP input dataset, Each row represent a species alignment, where the first three columns are GENE/Transcript name, Species Name, and Binary Phenotype Status (1/ 0/ NA) respectively. From the fourth column till end, input data have alignments for a genomic region, consistiting of gaps(-) or nucleotide bases. Following is a representation of GAP input:

```bash

Gene_ID           species_name  phenotype_status  alignment_columns
ENSMUSG000059970  nile_tilapia  0                 - - - A - T - G C A A T C G C T A -
ENSMUSG000059970  mus_musculus  1                 A C G A C T A G C A A T C G C T A C
ENSMUSG000059970  sloth         NA                - - - A - T - G C A - - - - - C T A 
```




# Tree Feature Extraction from Newick Format Tree
This Python script extracts tree features from a Newick format tree using the ETE Toolkit (ete3). It generates feature vectors representing species presence/absence within different branches of the tree.


**Requirements**
```bash
  pip install ete3
  pip install pandas
```
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
This command executes the `runGeneArchi.R` script. It tries different NN architectures based on the input parameteres to map multi-species genetic data to binary phenotypes. Architectures with minimum CV error are displayed alonside the manhattan plot containing nucleotide positional importance for exons within the gene(s).

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
Rscript "/users/username/downloads/_GAP/runGeneArchi.R" "/users/username/downloads/_GAP/" FALSE 0
```
which runs the script to identify NN architectures with minimum CV error with alignment-only dataset for all 0-hidden layered architectures.


## Run selected architecture on a range of genes

This command-line interface (CLI) command executes the `run_associated.R` script. It applies NN architecture zero CV error found from previous function on a specified gene list against multi-species genetic data to map binary phenotypes. The training results are stored in the `\results` folder under parent directory.

**Command Structure**

```bash
Rscript <path_to_run_associated.R> <path_source> boolean_tree_flag num_hidden_layers starting_gene ending_gene
```
This method requires execution from the shell due to resource constraints and follows a specific command structure with parallel processing:

**Arguments**:
  - `path_source`: Source of the downloaded tool, i.e. `/users/username/downloads/_GAP/`.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `num_hidden_layers`: Hidden layer count and nodes for each layers; for hidden layers < 2, use any integer, where 0 means no hidden layer, and integer n means single hidden layer with n nodes, for hidden layers >= 2, use comma separated integers, i.e., 8,1,4; meaning 3 hidden-layered architecture with 1st, 2nd, and 3rd layer having 8, 1, and 4 nodes respectively. 
  - `starting_gene`: specifies the gene number from where processing will start.
  - `ending_gene`: species the ending gene number of processing.

**Sample Command**:
```bash 
Rscript "/users/username/downloads/_GAP/run_associated.R" "/users/username/downloads/_GAP/" TRUE 1 21030 21040.
```
where the tool is located under the directory `/users/username/downloads/_GAP/` and this command will execute all the single-hidden layered architectures for the specified gene list and store the results under the specified directory. To run list of genes with hidden layers >= 2, command of following format can be used where a 2 hidden layer architecture is implemented for all the genes in the list with 1st and 2nd layer having 15 and 4 nodes respectively.

**Sample Command**:
```bash 
Rscript "/users/username/downloads/_GAP/run_associated.R" "/users/username/downloads/_GAP/" TRUE 15,4 21030 21040.
```
## List Associated Genes

This command-line interface (CLI) command executes the `findAssociated.R` script. It analyses the training results from previous function on a provided list of genes to find genes with zero CV error, the list of such genes is stored in the file `associated_genes.csv`. under the designated folder `results` in parent directory.

**Command Structure**

```bash
Rscript <path_to_findAssociated.R> <path_source> boolean_tree_flag starting_gene ending_gene
```

**Arguments**:
  - `path_source`: Source of the downloaded tool, i.e. `/users/username/downloads/_GAP/`.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `starting_gene`: specifies the start of the range to be checked.
  - `ending_gene`: species the end of the range to be checked.

**Sample Command**:
```bash 
Rscript "/users/username/downloads/_GAP/findAssociated.R" "/users/username/downloads/_GAP/" FALSE 21030 21040.
```
where the tool is located under the directory `/users/username/downloads/_GAP/` and this command will find and return (if any) of the genes with the range [21030,21040] has minimum CV scores, where the tree_features are not used.

**Note**
1. For Windows OS, replace Rscript with the path to Rscript.exe while running the scripts.
2. To run it for any gene beyond GULO, replace the GULO number to any gene of choice in `runGeneArchi.R`.
