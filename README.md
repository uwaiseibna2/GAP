# GAP (Genotype-phenotype Association Predictor)

GAP is a software implementation of Islam, Campelo dos Santos, Kanjilal, and Assis' (2024) learning genotype-phenotype associations from gaps in multi-species sequence alignments.

In particular, GAP employs a neural network to predict the presence or absence of a phenotype solely from gaps in a multiple alignment, with optional consideration of phylogenetic relationships from a user-defined species tree. GAP can be employed for three distinct problems: predicting a phenotype in one or more species from a known associated region, pinpointing which specific genomic positions within such a region are associated with a phenotype, and extracting a set of candidate genomic regions associated with a phenotype. 

For queries, email: uislam2022@fau.edu.

# Citing GAP

If you use GAP, please cite Islam UI, Campelo dos Santos AL, Kanjilal R, and Assis R. Learning genotype-phenotype associations from gaps in multi-species sequence alignments, under review (2024). 

# Getting Started

Running GAP requires the installation of R (version > 4.0) and Python (version > 3.9). GAP installs the required R and Python packages from within the script. If GAP encounters an error while installing a package, it can be installed manually.
Commands to install R and Python packages manually from respective terminals:
```bash
# R package installation command
install.packages("package_name")

# Python package installation command
pip3 install package_name
```

# GAP Input
GAP takes in two required input files and an optional third input file. 

**Input 1**: Tab-delimited text file with _n_ rows and two columns, where _n_ is the number of species. The first column should contain the species name, and the second column should refer to the phenotype status in that species, with 0 indicating absence, 1 indicating presence, and NA indicating unknown. 

**Input 2**: FASTA file containing multiple sequence alignments for the _n_ species at _g_ genomic regions. Each header should include the species name followed by a space and an identifier for the genomic region (e.g., gene ID, gene name, genomic coordinates). Each region must contain sequences for all _n_ species. If the region is entirely absent in a species, then a sequence of gaps "-" can be used for that species.

**Input 3 (optional)**: A text file containing the phylogenetic tree of the n species in Newick format. No distances should be included, and the species name should match the first input file.




# Running GAP

GAP has three functions:

- PredictSpecies - predicts the presence or absence of a phenotype in one or more species from a known association with a genomic region
- PredictPositions - predicts which positions within a genomic region are associated with a phenotype
- PredictGenes - predicts genomic regions that are associated with a phenotype

GAP is run in a terminal. Below are detailed instructions for using each of the GAP functions.

## PredictSpecies

The PredictSpecies function takes the paths to the GAP input files and the identifier of a genomic region. It trains a neural network on the specified genomic region, beginning with the simplest architecture and increasing in complexity, and stops if it identifies an architecture with a cross-validation error of zero. Otherwise, it trains on all architectures and selects the simplest architecture with the minimum cross-validation error. Then, it uses the selected model to predict whether the phenotype of interest is absent (0) or present (1) in each of the species with unknown status (NA) from the first input file. The output of this function is a tab-delimited file named `Predictions.csv` in the `results` folder, which contains species genomic region, names, and predictions for unknown species. 

**Command Structure**

```bash
Rscript PredictSpecies.R <path_source> boolean_tree_flag region_id <path_input_1> <path_input_2> <path_input_3>
```
**Arguments**:
- `path_source`: Source of the GAP directory
- `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training
- `region_id`: Identifier for the genomic region on which to train GAP
- `path_input_1`: Path to Input 1 file described above in the `GAP Input` section
- `path_input_2`: Path to Input 2 file described above in the `GAP Input` section
- `path_input_3`: (optional) path to the user-defined phylogenetic tree as `Input 3` described above in `GAP Input`. This input must be provided if `boolean_tree_flag` is set to `TRUE`.

    
## PredictPositions

The PredictPositions function takes the paths to the GAP input files and the identifier of a genomic region like the PredictSpecies function. Similarly, It trains a neural network on the specified genomic region, beginning with the simplest architecture and increasing in complexity, and stops if it identifies an architecture with a cross-validation error of zero. It trains on all architectures and selects the simplest architecture with the minimum cross-validation error if an architecture with a cross-validation error of zero is not found. GAP then uses the selected model to generate a tab-delimited file containing predictive importance for each position in the genomic region, with the position number in the first column and the Benjamini-Hochberg-adjusted p-value corresponding to predictive importance in the second column. This file with positional importance is stored in the file `PositionalPVals.csv` under the `results` directory.

**Command Structure**

```bash
Rscript PredictPositions.R <path_source> boolean_tree_flag region_id <path_input_1> <path_input_2> <path_input_3>
```
**Arguments**:
- `path_source`: Source of the GAP directory
- `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training
- `region_id`: Identifier for the genomic region on which to train GAP
- `path_input_1`: Path to Input 1 file described above in the `GAP Input` section
- `path_input_2`: Path to Input 2 file described above in the `GAP Input` section
- `path_input_3`: (optional) path to the user-defined phylogenetic tree as `Input 3` described in `GAP Input`. This input must be provided if `boolean_tree_flag` is set to `TRUE`.


## PredictGenes

The PredictGenes function takes as input the paths to the GAP input files described in the GAP Input. Unlike the previous functions, it does not require a genomic region, as it works on all the genomic regions present in the second input file. GAP examines these genomic regions to identify the potential association between each genomic region and the phenotype. The set of genomic regions identified as associated is documented in the output file `associated.csv` under the `results` folder in the parent directory.

**Command Structure**

```bash
Rscript PredictGenes.R <path_source> boolean_tree_flag <path_input_1> <path_input_2> <path_input_3>
```

**Arguments**
  - `path_source`: Source of the GAP directory.
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training.
  - `path_input_1`: Path to Input 1 file described above in the `GAP Input` section.
  - `path_input_2`: Path to Input 2 file described above in the `GAP Input` section.
  - `path_input_3`: (optional) path to the user-defined phylogenetic tree as `Input 3` described in `GAP Input`. This input must be provided if `boolean_tree_flag` is set to `TRUE`.

# Example Application of GAP

Commands to run the three GAP functions mentioned earlier are discussed below. Sample commands are provided for each function alongside explanations for alignment-only and tree-features-based approaches. The `phylogeny.txt` file under the `data-raw` directory contains a sample Newick formatted tree that can be followed while specifying the user-defined phylogenetic tree. To use the tree-features used in GAP analysis, replace the path to phylogeny features with the term `default`; in the following commands, `data-raw/phylogeny.txt` can be replaced by `default` to use these precomputed phylogeny features. 

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
This command runs the function PredictSpecies to identify optimal neural network architecture by exploring different parameters as described above. The simplest architecture with minimum CV error is used to predict the phenotype status of unknown species listed in the first input file. Note that the first set of commands excludes the tree features, whereas the latter set includes them, and the lengths of the arguments passed are different. 

## PredictPositions
**Sample Commands**:

- Alignment-only approach
  ```bash
  #unix-based OS
  Rscript PredictPositions.R ./ FALSE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictPositions.R ./ FALSE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa
  ```
- Tree-features included
  ```bash
  #unix-based OS
  Rscript PredictPositions.R ./ TRUE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictPositions.R ./ TRUE ENSMUST00000059970 data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
  ```
This command runs the function PredictPositions to identify optimal neural network architecture by exploring different parameters as described above. The simplest architecture with minimum CV error is used to measure the predictive importance of individual positions within the sequence of the given genomic region. All the positions within the sequence are listed alongside their p-values, and positions having p-values<= 0.05 are displayed separately for this architecture. 


## PredictGenes
**Sample Commands**:

- Alignment-only approach
    ```bash
    #unix-based OS
    Rscript PredictGenes.R "./" FALSE data-raw/species.txt data-raw/sample-dataset.fa 
    
    #Windows OS
    'C:/Program Files/.../Rscript.exe' PredictGenes.R ./ FALSE data-raw/species.txt data-raw/sample-dataset.fa 
    ```

- Tree-features included
    ```bash
    #unix-based OS
    Rscript PredictGenes.R "./" TRUE data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
    
    #Windows OS
    'C:/Program Files/.../Rscript.exe' PredictGenes.R ./ TRUE data-raw/species.txt data-raw/sample-dataset.fa data-raw/phylogeny.txt
    ```
This command runs the function PredictGenes to identify the genomic region(s) present in the first input file that can be associated with the given phenotype based on the simplest architecture with minimum cross-validation error. Regions identified as associated are listed as the output of this function.
