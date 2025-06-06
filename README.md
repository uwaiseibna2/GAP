# GAP (Genotype-phenotype Association Predictor)

GAP is a software implementation of Islam, Campelo dos Santos, Kanjilal, and Assis' (2025) learning genotype-phenotype associations from gaps in multi-species sequence alignments.

In particular, GAP employs a neural network to predict the presence or absence of a phenotype solely from gaps in a multiple alignment, with optional consideration of phylogenetic relationships from a user-defined species tree. GAP can be employed for three distinct problems: predicting a phenotype in one or more species from a known associated region, pinpointing which specific genomic positions within such a region are associated with a phenotype, and extracting a set of candidate genomic regions associated with a phenotype. 

For queries, email: uislam2022@fau.edu.

# Citing GAP

If you use GAP, please cite Islam UI, Campelo dos Santos AL, Kanjilal R, and Assis R. Learning genotype-phenotype associations from gaps in multi-species sequence alignments, Briefings in Bioinformatics, 2025. 

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

**Input 3 (optional)**: Text file containing the phylogenetic tree of the _n_ species in Newick format. No distances should be included, and the species names should match those in the first input file.

# Running GAP

GAP has three functions:

- PredictSpecies - predicts the presence or absence of a phenotype in one or more species from a known association with a genomic region
- PredictPositions - predicts which positions within a genomic region are associated with a phenotype
- PredictRegions - predicts genomic regions that are associated with a phenotype

GAP is run in a terminal. Below are detailed instructions for using each of the GAP functions.

## PredictSpecies

The PredictSpecies function takes as input the paths to the GAP input files and the identifier of a genomic region. It trains a neural network on the specified genomic region, beginning with the simplest architecture and increasing in complexity, and stops if it identifies an architecture with a cross-validation error of zero. Otherwise, it trains on all architectures and selects the simplest architecture with the minimum cross-validation error. Then, it uses the selected model to predict whether the phenotype of interest is absent (0) or present (1) in each of the species with unknown status (NA) from the first input file. The output of this function is a tab-delimited file named `predictions.csv` in the `results` folder, which contains the names of species with unknown status (NA) in the first column, and predictions for each of these species (0 for absence, 1 for presence) in the second column. 

**Command Structure**

```bash
Rscript PredictSpecies.R <path_source> boolean_tree_flag region_id <path_input_1> <path_input_2> <path_input_3> <cv_fold>
```
**Arguments**:
- `path_source`: Source of the GAP directory
- `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training
- `region_id`: Identifier for the genomic region on which to train GAP
- `path_input_1`: Path to Input 1 file described above in the `GAP Input` section
- `path_input_2`: Path to Input 2 file described above in the `GAP Input` section
- `path_input_3`: (optional) path to the user-defined phylogenetic tree as `Input 3` described above in `GAP Input`. This input must be provided if `boolean_tree_flag` is set to `TRUE`.
- `cv_fold`: (optional) determines n of n-fold cross-validation where the default used is leave-one-out, 5-fold CV is recommended for faster computation in case of large number of observations/species. 
  
## PredictPositions

The PredictPositions function takes as input the paths to the GAP input files and the identifier of a genomic region. It trains a neural network on the specified genomic region, beginning with the simplest architecture and increasing in complexity, and stops if it identifies an architecture with a cross-validation error of zero. Otherwise, it trains on all architectures and selects the simplest architecture with the minimum cross-validation error. Then, it obtains a _p_ value denoting the predictive importance of each position in the sequence of the genomic region. The output of this function is a tab-delimited file named `PositionalPVals.csv` in the `results` folder, which contains positions in the first column, and _p_ values for these positions in the second column. 

**Command Structure**

```bash
Rscript PredictPositions.R <path_source> boolean_tree_flag region_id <path_input_1> <path_input_2> <path_input_3> <cv_fold>
```
**Arguments**:
- `path_source`: Source of the GAP directory
- `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training
- `region_id`: Identifier for the genomic region on which to train GAP
- `path_input_1`: Path to Input 1 file described above in the `GAP Input` section
- `path_input_2`: Path to Input 2 file described above in the `GAP Input` section
- `path_input_3`: (optional) path to the user-defined phylogenetic tree as `Input 3` described in `GAP Input`. This input must be provided if `boolean_tree_flag` is set to `TRUE`.
- `cv_fold`: (optional) determines n of n-fold cross-validation where the default used is leave-one-out. 5-fold CV is recommended for faster computation in case of a large number of observations/species. 

## PredictRegions

The PredictRegions function takes as input the paths to the GAP input files. Unlike the previous functions, it does not require a genomic region, as it works on all the genomic regions present in the second input file. It trains a neural network on each genomic region, beginning with the simplest architecture and increasing in complexity, and stops if it identifies an architecture with a cross-validation error of zero. Otherwise, it trains on all architectures and selects the simplest architecture with the minimum cross-validation error. Then, the minimum cross-validation error across all genomic regions is identified, and those regions with this cross-validation error are predicted to be putatively associated with the phenotype. The output of this function is a file named `associated.csv` in the `results` folder, which contains the identifiers for all predicted genomic regions listed in a single column.

**Command Structure**

```bash
Rscript PredictRegions.R <path_source> boolean_tree_flag <path_input_1> <path_input_2> <path_input_3>
```

**Arguments**
  - `path_source`: Source of the GAP directory
  - `boolean_tree_flag`: Boolean (TRUE/FALSE) for tree feature inclusion in model training
  - `path_input_1`: Path to Input 1 file described above in the `GAP Input` section
  - `path_input_2`: Path to Input 2 file described above in the `GAP Input` section
  - `path_input_3`: (optional) path to the user-defined phylogenetic tree as `Input 3` described in `GAP Input`. This input must be provided if `boolean_tree_flag` is set to `TRUE`.

# Example Application of GAP

This example demonstrates the application of GAP to predict the association between the _Gulo_ gene and vitamin C synthesis with data from 59 vertebrates, following the strategy outlined in Islam _et al._ 2024. For enhanced computational efficiency, the second input file contains multiple alignments from  100 of the 22,476 genes used in the manuscript. Before beginning, ensure that the working directory is set to the directory containing GAP functions and input files.  

## PredictSpecies

To predict vitamin C synthesis in the 25 species with unknown status from the _Gulo_ multiple alignment, type:
  ```bash
  #unix-based OS
  Rscript PredictSpecies.R ./ FALSE ENSMUST00000059970 input/species.txt input/sample-dataset.fa
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictSpecies.R ./ FALSE ENSMUST00000059970 input/species.txt input/sample-dataset.fa
  ```
The output of this function is a tab-delimited file named `predictions.csv` in the `results` folder, which contains the names of species with unknown status (NA) in the first column, and predictions for each of these species (0 for absence, 1 for presence) in the second column. 

To predict vitamin C synthesis in the 25 species with unknown status from the _Gulo_ multiple alignment and phylogenetic tree relating the species with 5-fold cross-validation, type:
  ```bash
  #unix-based OS
  Rscript PredictSpecies.R ./ TRUE ENSMUST00000059970 input/species.txt input/sample-dataset.fa input/phylogeny.txt 5
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictSpecies.R ./ TRUE ENSMUST00000059970 input/species.txt input/sample-dataset.fa input/phylogeny.txt 5
  ```

## PredictPositions

To obtain _p_ values of positions in the _Gulo_ gene for predicting vitamin C synthesis from the multiple alignment, type: 
  ```bash
  #unix-based OS
  Rscript PredictPositions.R ./ FALSE ENSMUST00000059970 input/species.txt input/sample-dataset.fa
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictPositions.R ./ FALSE ENSMUST00000059970 input/species.txt input/sample-dataset.fa
  ```
To obtain _p_ values of positions in the _Gulo_ gene for predicting vitamin C synthesis from the multiple alignment and phylogenetic tree with 5-fold cross-validation, type: 
  ```bash
  #unix-based OS
  Rscript PredictPositions.R ./ TRUE ENSMUST00000059970 input/species.txt input/sample-dataset.fa input/phylogeny.txt 5
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictPositions.R ./ TRUE ENSMUST00000059970 input/species.txt input/sample-dataset.fa input/phylogeny.txt 5
  ```

The output of this function is a tab-delimited file named `PositionalPVals.csv` in the `results` folder, which contains positions in the first column, and _p_ values for these positions in the second column. 

## PredictRegions

To predict genomic regions that are associated with vitamin C synthesis from the multiple alignment, type:
  ```bash
  #unix-based OS
  Rscript PredictRegions.R "./" FALSE input/species.txt input/sample-dataset.fa 
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictRegions.R ./ FALSE input/species.txt input/sample-dataset.fa 
  ```

To predict genomic regions that are associated with vitamin C synthesis from the multiple alignment and phylogenetic tree, type:
  ```bash
  #unix-based OS
  Rscript PredictRegions.R "./" TRUE input/species.txt input/sample-dataset.fa input/phylogeny.txt
  
  #Windows OS
  'C:/Program Files/.../Rscript.exe' PredictRegions.R ./ TRUE input/species.txt input/sample-dataset.fa input/phylogeny.txt
  ```

The output of this function is a file named `associated.csv` in the `results` folder, which contains the identifiers for all predicted genomic regions in a single column.

## Shiny App

GAP is also available as shiny app with reduced functionality. To run the GAP shiny app, download all files in the repository and simply run the `Shiny_app.R` script from Rstudio and the interactive GUI can be used to explore all three functions with limited applications, but we recommend GAP in a terminal environment for extended use.

<img width="450" alt="GAP GUI" src="https://github.com/user-attachments/assets/e3e38658-fa25-4212-8250-7b9aa2cf949a" />

