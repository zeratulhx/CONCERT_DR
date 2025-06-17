# CONCERTDR: CONtext-aware Cellular and Tissue-specific Expression for Drug Repurposing

## Overview

CONCERTDR is an R package that provides a comprehensive suite of tools for drug response data analysis, with a particular focus on working with the Connectivity Map (CMap) database. It enables researchers to efficiently extract, process, and match drug response signatures from their experimental data with reference profiles from CMap.

## Installation

### Dependencies

Before installing CONCERTDR, ensure you have these prerequisites:

1. R (version 4.1.0 or higher)
2. Bioconductor packages - will be installed automatically if missing
3. Required R packages (automatically installed with the package):
   - BiocManager
   - data.table
   - cmapR
   - methods
   - stats
   - utils

### Installation Steps

```r
# Install devtools if not already installed
install.packages("devtools")

# Install CONCERTDR from GitHub
devtools::install_github("zeratulhx/CONCERT_DR")
```

## Required Data Files

To use the full functionality of this package, you need access to the following CMap data files:

1. **siginfo_beta.txt**: Contains metadata about all signatures in the database
2. **geneinfo_beta.txt**: Contains information about all genes in the database
3. ***.gctx**: Contains the expression data for perturbation signatures
4. Optional***repurposing_drugs.txt**: Contains the clinical state of drugs for step 5.

### Obtaining CMap Data 

1. Navigate to [Connectivity Map Portal (clue.io)](https://clue.io/releases/data-dashboard)
2. Download the required files:
   - CMap L1000 Data: Level 5 (GCTx): `level5_beta_all_n1201944x12328.gctx`(recommended)
     or
     any of the GCTx file in your interest
   
   - Metadata files: `siginfo_beta.txt` and `geneinfo_beta.txt`
   
The description of the differences between these files can be found in the 'LINCS2020 Release Metadata Field Definitions.xlsx' provided in the website.

Alternatively, you can use our helper script to download the files:

```r
# Create a directory for the data
dir.create("cmap_data")

# Source the download script
source(system.file("scripts", "download_databases.R", package = "CONCERTDR"))

# This will guide you through the download process
download_cmap_data("cmap_data")
```

Also, the repurposing_drugs.txt is available at [https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt] on Broad Institute website.

## Complete Workflow Tutorial

### Step 1: Subsetting CMap for signatures of interest

First, we need to extract experimental parameters from the CMap database. This will create a file that contains all available time points, dosages, and cell lines for certain type of perturbation.

For this pipeline, you can simply filter out the parameters you are interested in by making a subset of the 'siginfo_beta.txt' file. Currently only four parameters are supported:
  - `pert_type`: Type of perturbation (e.g., "trt_cp" for chemical treatments)
  - `pert_itime`: Time of treatment (e.g., "6h", "24h")
  - `pert_idose`: Dosage of the treatment (e.g., "10uM", "100nM")
  - `cell`: Cell line used in the experiment (e.g., "A375", "MCF7")

For interactively selecting these parameters or see what's available inside the database, you can simply run:
```r
library(CONCERTDR)

filtered_siginfo <- subset_siginfo_beta("path/to/siginfo_beta.txt",
                                interactive = TRUE,
                                verbose = TRUE) 
```
or you can use the filters parameter in non-interactive mode:
```r
filtered_siginfo <- subset_siginfo_beta(
   "path/to/siginfo_beta.txt",
   interactive = FALSE,
   filters = list(
     pert_type = "trt_cp",
     pert_itime = c("6 h", "24 h"),
     pert_idose = "10 uM",
     cell_iname = c("A375", "MCF7")
     ))
```



You can also create a file by setting 'output=YOUR_PATH'. You can edit this file to select specific parameters or keep all of them.

The purt_type parameters is defined by LINCS2020 as follows:

![image](https://github.com/user-attachments/assets/02ef148d-736b-4c02-92b5-5fde0935db17)



### Step 2: Prepare Your Signature File

Create a tab-delimited signature file containing your gene signature. The file must have at least these two columns:
- `Gene`: Gene symbols matching those in CMap
- `log2FC`: Log2 fold change values (positive for up-regulated, negative for down-regulated)

Example `signature.txt`:
```
Gene    log2FC
STAT3    2.5
TP53    -1.8
MYC     3.2
BRCA1   -2.1
```

### Step 3: Extract Data from CMap

Now you can extract the expression data from CMap for your selected parameters:

```r


reference_df <- extract_cmap_data_from_siginfo(
  geneinfo_file = "path/to/geneinfo_beta.txt",
  siginfo_file = "path/to/siginfo_beta.txt",
  gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx",
  filter_quality = TRUE, # Use hiq=1 high quality signatures only
)
```



### Step 4: Match Your Signature Against Reference Profiles

```r
# Match your signature against the reference data
results <- process_signature_with_df(
  reference_df = reference_df,
  signature_file = "signature.txt",
  output_dir = "results",
  methods = c("ks", "xsum"),  # Choose methods 
  topN = 100,                 # For XCos/XSum methods
  permutations = 1000         # Statistical robustness
)

# Interactive exploration of results
print(results)           # Overview of results
summary(results)         # Top compounds across methods
head(results$results$ks) # Top hits for KS method

```

### Step 5(Optional): Annotate your drugs from cmap id

```r
annotate <- annotate_drug_results(
  results_df = results$results$ks,
  sig_info_file = filtered_siginfo,
  comp_info_file = "path/to/compoundinfo_beta.txt",
  drug_info_file = "/path/to/repurposing_drugs_20200324.txt"
)
```

## Available Signature Matching Methods

CONCERTDR provides access to various signature matching methods:

1. **KS Score** (`ks`): Kolmogorov-Smirnov statistic based method
2. **XCos Score** (`xcos`): Extreme cosine similarity method
3. **XSum Score** (`xsum`): Extreme sum score method
4. **GSEA Weight 0 Score** (`gsea0`): Gene Set Enrichment Analysis without weighting
5. **GSEA Weight 1 Score** (`gsea1`): Gene Set Enrichment Analysis with weight 1
6. **GSEA Weight 2 Score** (`gsea2`): Gene Set Enrichment Analysis with weight 2
7. **Zhang Score** (`zhang`): Method from Zhang et al.

## Function Reference

### Parameter Extraction and Management

| Function | Description |
|----------|-------------|
| `extract_cmap_parameters()` | Extracts experimental parameters from siginfo_beta.txt |
| `create_cmap_config_template()` | Creates a template configuration file |
| `read_cmap_config()` | Reads a configuration file |
| `generate_combinations_from_config()` | Generates parameter combinations |

### Data Extraction

| Function | Description |
|----------|-------------|
| `extract_cmap_data_from_config()` | Extracts data for all combinations into a data frame |
| `process_combinations()` | Processes combinations to extract data into individual files |
| `process_combinations_file()` | Processes combinations from a file (for SLURM jobs) |

### Signature Matching

| Function | Description |
|----------|-------------|
| `process_signature_with_df()` | Matches signature against reference data frame |
| `run_cmap_workflow()` | Runs the complete workflow from config to matching |

### Utilities

| Function | Description |
|----------|-------------|
| `select_items_interactive()` | Interactively select items from a list |
| `select_methods_interactive()` | Interactively select matching methods |
| `interactive_cmap_workflow()` | Run an interactive command-line workflow |
| `generate_slurm_script()` | Generate a SLURM submission script for HPC |


## License

This package is licensed under the MIT License - see the LICENSE file for details.
