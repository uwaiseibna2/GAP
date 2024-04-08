# GAP (Genotype-phenotype Association Predictor)

GAP is a software implementation of Islam, Campelo dos Santos, Kanjilal, and Assis' (2024) learning genotype-phenotype associations from gaps in multi-species sequence alignments.

In particular, GAP employs a neural network to predict the presence or absence of a phenotype solely from gaps in a multiple alignment, with optional consideration of phylogenetic relationships from a user-input species tree. GAP can be employed for three distinct problems: predicting a phenotype in one or more species from a known associated region, pinpointing which specific genomic positions within such a region are associated with a phenotype, and extracting a set of candidate genomic regions associated with a phenotype. 

For queries, email: uislam2022@fau.edu.

# Citing GAP

If you use GAP, then please cite: Islam UI, Campelo dos Santos AL, Kanjilal R, and Assis R. Learning genotype-phenotype associations from gaps in multi-species sequence alignments. Under review (2024). 

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

**Input 3 (optional)**: phylogenetic tree of the n species in Newick format. No distances should be included, example format: "((gorilla,(chimp,human),baboon),orangutan);", where n=5.




# Running GAP


GAP has three functions:
- PredictSpecies - predicts presence or absence of a phenotype in one or more species from a known association with a genomic region
- PredictPositions - predicts which positions within a genomic region are associated with a phenotype
- PredictGenes - predicts genomic regions that are associated with a phenotype

GAP is run in a terminal. Below are detailed instructions for using each of the GAP functions.

## PredictSpecies

The PredictSpecies function takes in the paths to the input files and the identifier of a genomic region to be used for predictions. It trains a neural network on the specified genomic region, beginning with the simplest architecture and increasing in complexity, and stops if it identifies an architecture with a cross-validation error of zero. Otherwise, it trains on all architectures and selects the simplest architecture with the smallest cross-validation error. Then it uses the selected model to predict whether the phenotype of interest is absent (0) or present (1) in each species with unknown status (NA) from the first input file. The output of this function is a tab-delimited file named `Predictions.csv` in the `results` folder, which contains species names in the first column and predictions in the second column. 

**Command Structure**

```bash
Rscript PredictSpecies.R <path_source> boolean_tree_flag region_id <path_input_1> <path_input_2> <path_input_3>
```
**Arguments**:
- `path_source`: Source of the GAP directory
- `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training
- `region_id`: Identifier for genomic region on which to train GAP
- `path_input_1`: Path to Input 1 file described in the `GAP Input` section above
- `path_input_2`: Path to Input 2 file described in the `GAP Input` section above
- `path_input_3`: (optional) path to user-defined phylogenetic tree as `Input 3` described in `GAP Input`. If `boolean_tree_flag` is set to `TRUE` this input must be provided.

    
## PredictPositions

The predictPositions function output a tab-delimited file containing predictive importance for each position in the genomic region, with position number in the first column and the Benjamini-Hochberg-adjusted p-value corresponding to predictive importance in the second column. Please note that, that the input number of genomic regions g=1 for this function with the transcript_id used for the previous function. 

**Command Structure**

```bash
Rscript PredictPostions.R <path_source>
```

**Arguments**:
  - `path_source`: Source of the GAP directory.


## PredictGenes

The PredictGenes function idenfies takes the input files described in the GAP Input. GAP examines each transcript present in the Input 2 file, to identify potential association between the transcript and the said phenotype. The set of transcripts identified as associated is documented in the output file `associated.csv` under the results folder in the parent directory.

**Command Structure**

```bash
Rscript PredictGenes.R <path_source> boolean_tree_flag <path_input_1> <path_input_2> <path_input_3>
```

**Arguments**
  - `path_source`: Source of the GAP directory.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `path_input_1`: Path to Input 1 file described in the `GAP Input` section above.
  - `path_input_2`: Path to Input 2 file described in the `GAP Input` section above.
  - `path_input_3`: (optional) path to user-defined phylogenetic tree as `Input 3` described in `GAP Input`. If `boolean_tree_flag` is set to `TRUE` this input must be provided.

# Example Application of GAP

Commands to run the three GAP functions are discussed here. For each function, sample commands are provided alongside explanation for both alignment-only and tree-features based approach.

## PredictSpecies
**Sample Commands**:

- Alignment-only approach
  ```bash
  #unix-based OS
  Rscript PredictSpecies.R ./ FALSE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictSpecies.R ./ FALSE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa
  ```
- Tree-features included
  ```bash
  #unix-based OS
  Rscript PredictSpecies.R ./ TRUE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictSpecies.R ./ TRUE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
  ```
which runs the script to identify neural network architectures with minimum CV error by exploring different architecture, progressing from 0-hidden layer architecutre to 3-hidden layer architectures and stops the moment it finds an architecture with minimum CV error. Notice that the first set of commands exclude the tree features whereas the latter set include them. The predicted phenotypes for all species are stored in `Predictions.csv` under the `results` folder.

## PredictPositions
**Sample Commands**:

- Both alignment-only and tree-features included approach
  ```bash
  #unix-based OS
  Rscript PredictPositions.R ./ 
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictPositions.R ./ 
  ```
This function identifies positions within the sequence having p-values<= 0.05 within the alignment. PredictPositions needs the previous function implemented before in order to positions based on the identify optimal architecture idenfied in PredictSpecies. Results are stored in the `PositionalPvals.csv` file under the `results` directory.


## PredictGenes
**Sample Commands**:

- Alignment-only approach
    ```bash
    #unix-based OS
    Rscript PredictGenes.R "./" FALSE data-raw/species.txt data-raw/sample-dataset.fa data-raw/Transcript_list.txt 
    
    #Windows OS
    'C:/Program Files/.../Rscript.exe' PredictGenes.R ./ FALSE data-raw/species.txt data-raw/sample-dataset.fa data-raw/Transcript_list.txt
    ```

- Tree-features included
    ```bash
    #unix-based OS
    Rscript PredictGenes.R "./" TRUE data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
    
    #Windows OS
    'C:/Program Files/.../Rscript.exe' PredictGenes.R ./ TRUE data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
    ```
This function lists the associated genes based on the optimal architecture configuration and minimum CV error from the provided list of transcript ids. The associated transcript ids are stored in the `associated.csv` file under the `results` directory.
