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

GAP accepts multi-species sequence data as input in a specified format. In a GAP input dataset, Each row represent a species alignment, where the first three columns are gene / transcript name, species name, and binary phenotype status (1/ 0/ NA) respectively. From the fourth column till the end, input data have alignments for a genomic region, consistiting of gaps(-) or nucleotide bases. Following is a representation of GAP input:

```bash

gene_id             species_name  phenotype_status  alignment_columns
ENSMUSG00000059970  nile_tilapia  0                 - - - A - T - G C A A T C G C T A -
ENSMUSG00000059970  mus_musculus  1                 A C G A C T A G C A A T C G C T A C
ENSMUSG00000059970  sloth         NA                - - - A - T - G C A - - - - - C T A 
```

GAP input files have numbered filename where each file contains aforementioned formatted multi-species sequence alignment. Sample data located in data-raw directory contains alignment files where data for phenotype status-known and status-unknown species are stored in known_species and unknown_species respectively.


# GAP functions

## Predict Binary Phenotype Status
This command executes the `runGeneArchi.R` script. It tries different NN architectures based on the input parameteres to map multi-species genetic data to binary phenotypes. Architectures with minimum CV error are displayed alonside the manhattan plot containing nucleotide positional importance for exons within the gene(s). Navigate to the diretory of GAP and initiate a terminal from that directory to run these scripts.

**Command Structure**

```bash
Rscript runGeneArchi.R <path_source> boolean_tree_flag num_hidden_layers
```
**Arguments**:
  - `path_source`: Source of the downloaded tool, in this case `./` as we are already in the GAP directory.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `num_hidden_layers`: Hidden layer count {0,1,2,3}.
    
**Details**:
  - For `num_hidden_layers` > 0, explores permutations (1,2,3,..., number_of_species).
  - Displays a Manhattan plot for qualifying architectures.
    
**Sample Command**:
```bash 
Rscript runGeneArchi.R ./ FALSE 0
```
which runs the script to identify NN architectures with minimum CV error with alignment-only dataset for all 0-hidden layered architectures to predict binary phenotypes.


## Run an architecture on genomic region(s)

This command-line interface (CLI) command executes the `run_associated.R` script. It applies NN architecture zero CV error found from previous function on a specified gene list against multi-species genetic data to map binary phenotypes. The training results are stored in the `\results` folder under parent directory.

**Command Structure**

```bash
Rscript run_associated.R <path_source> boolean_tree_flag num_hidden_layers starting_gene ending_gene
```
This method requires execution from the shell due to resource constraints and follows a specific command structure with parallel processing:

**Arguments**:
  - `boolean_tree_flag`: As introduced on the previous section
  - `num_hidden_layers`: Hidden layer count and nodes for each layers; for hidden layers < 2, use any integer, where 0 means no hidden layer, and integer n means single hidden layer with n nodes, for hidden layers >= 2, use comma separated integers, i.e., 8,1,4; meaning 3 hidden-layered architecture with 1st, 2nd, and 3rd layer having 8, 1, and 4 nodes respectively. 
  - `starting_gene`: specifies the gene number from where processing will start.
  - `ending_gene`: species the ending gene number of processing.

**Sample Command**:
```bash 
Rscript run_associated.R ./ TRUE 1 21030 21040
```
where the tool is located under the current terminal directory and this command will execute all the single-hidden layered architectures for the specified gene list and store the results under the specified directory. To run list of genes with hidden layers >= 2, command of following format can be used where a 2 hidden layer architecture is implemented for all the genes in the list with 1st and 2nd layer having 15 and 4 nodes respectively.

**Sample Command**:
```bash 
Rscript run_associated.R ./ TRUE 15,4 21030 21040
```
## List Associated Genes

This command-line interface (CLI) command executes the `findAssociated.R` script. It analyses the training results from previous function on a provided list of genes to find genes with zero CV error, the list of such genes is stored in the file `associated_genes.csv`. under the designated folder `results` in parent directory.

**Command Structure**

```bash
Rscript findAssociated.R <path_source> boolean_tree_flag starting_gene ending_gene
```

**Arguments**
 - `boolean_tree_flag`: As introduced on the previous section.
 - `starting_gene`: As introduced on the previous section.
 - `ending_gene`: As introduced on the previous section.

**Sample Command**:
```bash 
Rscript findAssociated.R "./" FALSE 21030 21040.
```
where the tool is located under the current terminal directory and this command will find and return (if any) of the genes with the range [21030,21040] has minimum CV scores, where the tree_features are not used.

**Note**
1. For Windows OS, replace Rscript with the path to Rscript.exe while running the scripts.
2. To run it for any gene beyond GULO, replace the GULO number to any gene of choice in `runGeneArchi.R`.

# User-defined Phylogeny

GAP includes the phylogeny for the 59 species examined. However, GAP can generate user-specific custom phylogeny features using the `extract_tree_feats.py` script. This script takes a species tree in newick format and generates output is a .csv file in the parent directory.

Input newick formatted species tree:
```  
"((gorilla,(chimp,human),baboon),orangutan);"
```
To execute this python script navigate to the parent directory and simply execute the following command in a terminal.

**Command**

```bash
python3 extract_tree_features.py
```

**Notes**

1. Ensure the Newick tree string format is accurate.
2. Species names should be represented as comma-separated alphanumeric strings in the Newick tree string.
3. Replace variable `path_tree` with the newly generated tree_features.csv path in all three R scripts.

