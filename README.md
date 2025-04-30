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

# Step 1: Extract parameters from CMap siginfo file and generate configuration file
params <- extract_cmap_parameters("path/to/siginfo_beta.txt")

# Step 2: Process combinations to extract data
ref_df <- extract_cmap_data_from_config(
  config_file = "conf/cmap_options.conf",
  geneinfo_file = "/home/users/allstaff/pan.h/vast/DR_databases/geneinfo_beta.txt",
  siginfo_file = "/home/users/allstaff/pan.h/vast/DR_databases/siginfo_beta.txt",
  gctx_file = "/home/users/allstaff/pan.h/vast/DR_databases/level5_beta_all_n1201944x12328.gctx",
)

# Step 3: Match your drug response signature against reference profiles
process_signature_with_df(reference_df = ref_df,
                          signature_file = "signature.csv",
                          output_dir = "results",
                          methods="xsum",
                          topN=200)
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

