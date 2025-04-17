# DRnew: Drug Response Data Analysis Tools

## Overview

DRnew is an R package that provides a comprehensive suite of tools for drug response data analysis, with a particular focus on working with the Connectivity Map (CMap) database. It enables researchers to efficiently extract, process, and match drug response signatures from their experimental data with reference profiles from CMap.

## Installation

### Development Version

```r
# Install directly from GitHub
# install.packages("devtools")
devtools::install_github("zeratulhx/DR_pipeline")
```

## Dependencies

DRnew depends on the following packages:

- **BiocManager**: For installing Bioconductor packages
- **cmapR**: For working with CMap GCTX files
- **data.table**: For efficient data manipulation
- **RCSM**: For signature matching methods (installed automatically from GitHub)

## Usage

### Basic Workflow

```r
library(DRnew)

# Step 1: Extract parameters from CMap siginfo file
params <- extract_cmap_parameters("path/to/siginfo_beta.txt")

# Step 2: Generate combinations of time points, doses, and cell lines
combinations <- generate_combinations_from_config("conf/cmap_options.conf")

# Step 3: Process combinations to extract data
process_combinations(combinations,
                     output_dir = "output",
                     geneinfo_file = "path/to/geneinfo_beta.txt",
                     siginfo_file = "path/to/siginfo_beta.txt",
                     gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx")

# Step 4: Match your drug response signature against reference profiles
run_analysis_with_combinations(combinations_file = "combinations.txt",
                               sig_file = "signature.txt",
                               out_dir = "results")
```

## Acquiring CMap Data

The Connectivity Map (CMap) dataset can be downloaded from the [Broad Institute's CMap website](https://clue.io/). You'll need to register for an account to access the data.

Required files include:
- `siginfo_beta.txt`: Contains information about all signatures in the database
- `geneinfo_beta.txt`: Contains information about all genes in the database
- `level5_beta_trt_cp_n720216x12328.gctx`: Contains the expression data for all signatures

You can use the provided script to download these files:

```r
source(system.file("scripts", "download_databases.R", package = "DRnew"))
download_cmap_databases(dest_dir = "databases")
```

## Documentation

For more detailed documentation and examples, refer to the package vignettes:

```r
browseVignettes("DRnew")
```

