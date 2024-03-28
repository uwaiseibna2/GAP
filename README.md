# GAP (Genotype-phenotype Association Predictor)

GAP is a software implementation of Islam, Campelo dos Santos, Kanjilal, and Assis' (2024) machine learning framework for predicting genotype-phenotype associations from multi-species sequence alignments.

We introduce GAP, a novel machine learning framework for predicting a binary phenotype from gaps in a multi-species sequence alignment of a genomic region. In particular, GAP employs a neural network (NN) to predict the presence or absence of a phenotype solely from gaps in a multiple alignment, with optional consideration of phylogenetic relationships from a user-input species tree. GAP can be employed for three distinct problems: predicting a phenotype in one or more species from a known associated region, pinpointing which specific genomic positions within such a region are associated with a phenotype, and extracting a set of candidate genomic regions associated with a phenotype. Here, we demonstrate the utility of GAP by exploiting the well-known association between the gulonolactone (L-) oxidase *Gulo* gene and vitamin C synthesis.


For improvment suggestions and queries, email: uislam2022@fau.edu

# Cite GAP

Thank you for using GAP! 

Citation is appreciated: UI Islam, AL Campelo dos Santos, R Kanjilal, and R Assis' (2024) machine learning framework for predicting genotype-phenotype associations from multi-species sequence alignments 

# Getting Started

Running GAP scripts require installation of Python3 (version > 3.9) and R (version > 4.0) into the system. GAP installs the required R and Python packages from within the script. If GAP fails to install a package itself and users encounter a package installation error, they can install that particular package manually. 
Commands to install R and Python packages from respective terminals:
```bash
# R package installation command
install.packages("package_name")

# Python package installation command
pip3 install package_name
```

# GAP Input

**Input 1**: tab-delimited file with n rows and 2 columns, where n is the number of species. The first column should contain the species name, and the second column refer to the phenotype status in that species, with a 0 indicating absence, a 1 indicating presence, and a NA indicating unknown. 


**Input 2**: FASTA file containing multiple sequence alignments for the n species at g genomic regions. Each header should contain the species name followed by a space and then an identifier for the genomic region (e.g., gene ID, gene name, genomic coordinates). Each region must contain sequences for all n species. If the region is completely absent in a species, then a sequence of gaps "-" can be used for that species.  

**Input 3 (optional)**: phylogenetic tree of the n species in Newick format. No distances should be included.




# GAP functions

All GAP functions are implemented from the terminal. To run GAP commands open a terminal on the parent directory of GAP. In UNIX based Operating system (linux and macOS), Rscript can be used to execute an R script from terminal, However in Windows, users need to replace the Rscript with the path to `Rscript.exe`.

## PredictPhenotype
The predictSpecies function take in the path of the input files, and then output a tab-delimited file containing predictions for each unknown species, with species name in the first column and predicted phenotype status (0 for absent or 1 for present) in the second column. Note that, input number of genomic regions, g=1 for this function. This functions trains NN on the specified genomic region starting with simple architecture and proceeding to complex ones and stops if it find one with zero cross-validation error. 

**Command Structure**

```bash
Rscript PredictPheno.R <path_source> boolean_tree_flag transcript_id <path_input_1> <path_input_2> <path_input3>
```
**Arguments**:
  - `path_source`: Source of the downloaded tool/GAP directory.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `transcript_id`: Ensemble transcript ID for specifying the gene on which to train GAP, here we have added the transcript_id for GULO.
  - `path_input_1`: Path to the Input 1 file as described in the `GAP Input` section which list the phenotype status for a list of species.
  - `path_input_2`: Path to the Input 2 .fasta file containing cross-speices transcript alignments as described in the `GAP Input` section .
  - `path_input_3: Path to the Input 3 file containing phylogeny fetures as described in the `GAP Input` section, the default set of phylogeny features are located in `data-raw/tree-features.csv`

**Details**:
  - Output will be stored under the `results` folder in the parent directory.
    
**Sample Command**:
```bash
#unix-based OS
Rscript PredictPheno.R ./ FALSE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa data-raw/tree-features.csv

#Windows OS
'C:/Program Files/.../Rscript.exe' PredictPheno.R ./ FALSE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa data-raw/tree-features.csv
```
which runs the script to identify NN architectures with minimum CV error by exploring different architecture, progressing from 0-hidden layer architecutre to 3-hidden layer architectures and stops the moment it finds an architecture with minimum CV error. Notice that this command excludes the tree features.


## PredictPositions

The predictPositions function output a tab-delimited file containing predictive importance for each position in the genomic region, with position number in the first column and the Benjamini-Hochberg-adjusted p-value corresponding to predictive importance in the second column. Please note that, that the input number of genomic regions g=1 for this function with the transcript_id used for the previous function. This function requires results from the previous function to generate the positional importance, please make sure you have the results of the previous function before running this command.

**Command Structure**

```bash
Rscript PredictPostions.R <path_source>
```
This method requires execution from the shell due to resource constraints and follows a specific command structure with parallel processing:

**Arguments**:
  - `path_source`: Source of the downloaded tool/GAP directory.

**Sample Command**:
```bash
#unix-based OS
Rscript PredictPositions.R ./

#Windows OS
'C:/Program Files/.../Rscript.exe' PredictPositions.R ./ 
```
where the tool is located under the current terminal (unix-based OS) directory and this command will provide the positions importance for each of the positions in the transcript for the optimal architecture found through PreditPhenotype function.

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
  - `path_input_3: Path to the Input 3 file containing phylogeny fetures as described in the `GAP Input` section, the default set of phylogeny features are located in `data-raw/tree-features.csv`
  - `path_transcripts_list`: List of sample transcripts with a transcript_id in each row in a `.txt` file, if not provided, GAP will use the default sample list under `/data-raw` in parent directory.

**Sample Command**:
```bash
#unix-based OS
Rscript PredictGenes.R "./" FALSE data-raw/species.txt data-raw/sample-dataset.fa data-raw/tree-features.csv data-raw/transcript_list

#Windows OS
'C:/Program Files/.../Rscript.exe' PredictGenes.R ./ FALSE data-raw/species.txt data-raw/sample-dataset.fa data-raw/tree-features.csv data-raw/transcript_list
```
where the tool is located under the current terminal directory and this command will find and return (if any) the ones within the listed genes having minimum CV scores, where the tree_features are not used.

**Note**
1. The default phylogeny features are located at `data-raw/tree-features.csv`, in case of user-defined phylogeny (see instuctions below to generate), change the path accordingly. 

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



