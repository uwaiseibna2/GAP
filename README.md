# GAP (Genotype-phenotype Association Predictor)

GAP is a software implementation of Islam, Campelo dos Santos, Kanjilal, and Assis' (2024) machine learning framework for predicting genotype-phenotype associations from multi-species sequence alignments.

In particular, GAP employs a neural network (NN) to predict the presence or absence of a phenotype solely from gaps in a multiple alignment, with optional consideration of phylogenetic relationships from a user-input species tree. GAP can be employed for three distinct problems: predicting a phenotype in one or more species from a known associated region, pinpointing which specific genomic positions within such a region are associated with a phenotype, and extracting a set of candidate genomic regions associated with a phenotype. 

For queries, email: uislam2022@fau.edu.

# Citing GAP

If you use GAP, then please cite: Islam UI, Campelo dos Santos AL, Kanjilal R, and Assis R. A machine learning framework for predicting genotype-phenotype associations from multi-species sequence alignments. Under review (2024). 

# Getting Started

Running GAP requires installation of R (version > 4.0) and Python (version > 3.9). GAP installs the required R and Python packages from within the script. If GAP encounters an error during installation of a package, they can be installed manually.
Commands to install R and Python packages from respective terminals:
```bash
# R package installation command
install.packages("package_name")

# Python package installation command
pip3 install package_name
```

# GAP Input
GAP takes in two required input files and an optional third input file. 

**Input 1**: tab-delimited file with n rows and 2 columns, where n is the number of species. The first column should contain the species name, and the second column refer to the phenotype status in that species, with a 0 indicating absence, a 1 indicating presence, and a NA indicating unknown. 


**Input 2**: FASTA file containing multiple sequence alignments for the n species at g genomic regions. Each header should contain the species name followed by a space and then an identifier for the genomic region (e.g., gene ID, gene name, genomic coordinates). Each region must contain sequences for all n species. If the region is completely absent in a species, then a sequence of gaps "-" can be used for that species.  

**Input 3 (optional)**: phylogenetic tree of the n species in Newick format. No distances should be included.




# Running GAP


GAP has three functions:
- PredictSpecies - predicts presence or absence of a phenotype in one or more species from a known association with a genomic region
- PredictPositions - predicts which positions within a genomic region are associated with a phenotype
- PredictGenes - predicts genomic regions that are associated with a phenotype

GAP is run in a terminal. Below are detailed instructions for using each of the GAP functions.

## PredictSpecies
The predictSpecies function takes in the paths to the input files and the identifier of a genomic region to be used for predictions. It trains a neural network on the specified genomic region, beginning with the simplest architecture and increasing in complexity, and stops if it identifies an architecture with a cross-validation error of zero. Otherwise, it trains on all architectures and selects the simplest architecture with the smallest cross-validation error. Then it uses the selected model to predict whether the phenotype of interest is absent (0) or present (1) in each species with unknown status (NA) from the first input file. The output of this function is a tab-delimited file named `Predictions.csv` in the `results` folder, which contains species names in the first column and predictions in the second column. 

**Command Structure**

```bash
Rscript PredictSpecies.R <path_source> boolean_tree_flag region_id <path_input_1> <path_input_2> <path_input3>
```
**Arguments**:
- `path_source`: Source of the GAP directory
- `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training
- `region_id`: Identifier for genomic region on which to train GAP
- `path_input_1`: Path to Input 1 file described in the `GAP Input` section above
- `path_input_2`: Path to Input 2 file described in the `GAP Input` section above
- `path_input_3`: Path to Input 3 file described in `GAP Input` section above

**Details**:
  - Output file `Predictions.csv` will be stored under the `results` folder in the parent directory.
    
## PredictPositions

The predictPositions function output a tab-delimited file containing predictive importance for each position in the genomic region, with position number in the first column and the Benjamini-Hochberg-adjusted p-value corresponding to predictive importance in the second column. Please note that, that the input number of genomic regions g=1 for this function with the transcript_id used for the previous function. This function requires results from the previous function to generate the positional importance, please make sure you have the results of the previous function before running this command.

**Command Structure**

```bash
Rscript PredictPostions.R <path_source>
```
This method requires execution from the shell due to resource constraints and follows a specific command structure with parallel processing:

**Arguments**:
  - `path_source`: Source of the downloaded tool/GAP directory.
**Details**:
  - Output file `PositionalPvals.csv` will be stored under the `results` folder in the parent directory.

## PredictGenes

The predictGenes function should takes in input file as a list of transcripts to examine, and then output a list of predicted genomic regions with minimum CV error. Please note that, the input number of genomic regions is g>1 for this function.  

**Command Structure**

```bash
Rscript PredictGenes.R <path_source> boolean_tree_flag <path_input_1> <path_input_2> <path_input_3> <path_transcripts_list>
```

**Arguments**
  - `path_source`: Source of the downloaded tool/GAP directory.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `path_input_1`: Path to the Input 1 file as described in the `GAP Input` section which list the phenotype status for a list of species.
  - `path_input_2`: Path to the Input 2 .fasta file containing cross-speices transcript alignments as described in the `GAP Input` section .
  - `path_input_3`: Path to the Input 3 file containing phylogeny fetures as described in the `GAP Input` section, the default set of phylogeny features are located in `data-raw/tree-features.csv`
  - `path_transcripts_list`: List of sample transcripts with a transcript_id in each row in a `.txt` file, if not provided, GAP will use the default sample list under `/data-raw` in parent directory.

**Details**:
  - Output file `associated.csv` stored under the `results` folder in the parent directory will list associated genes.
## User-defined Phylogeny

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

## Example Application of GAP



