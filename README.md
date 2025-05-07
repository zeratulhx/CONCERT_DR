# CONCERT-DR: CONtext-aware Cellular and Tissue-specific Expression for Drug Repurposing

## Overview

CONCERT-DR is an R package that provides a comprehensive suite of tools for drug response data analysis, with a particular focus on working with the Connectivity Map (CMap) database. It enables researchers to efficiently extract, process, and match drug response signatures from their experimental data with reference profiles from CMap.

## Installation

### Dependencies

Before installing DRnew, ensure you have these prerequisites:

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

# Install DRnew from GitHub
devtools::install_github("zeratulhx/CONCERT_DR")
```

## Required Data Files

To use the full functionality of this package, you need access to the following CMap data files:

1. **siginfo_beta.txt**: Contains metadata about all signatures in the database
2. **geneinfo_beta.txt**: Contains information about all genes in the database
3. **level5_beta_trt_cp_n720216x12328.gctx**: Contains the expression data for all signatures

### Obtaining CMap Data 

1. Navigate to [Connectivity Map Portal (clue.io)](https://clue.io/)
2. Download the required files:
   - CMap L1000 Data: Level 5 (GCTx): `level5_beta_trt_cp_n720216x12328.gctx`
   - Metadata files: `siginfo_beta.txt` and `geneinfo_beta.txt`

Alternatively, you can use our helper script to download the files:

```r
# Create a directory for the data
dir.create("cmap_data")

# Source the download script
source(system.file("scripts", "download_databases.R", package = "CONCERT_DR"))

# This will guide you through the download process
download_cmap_data("cmap_data")
```

## Complete Workflow Tutorial

### Step 1: Extract Parameters from CMap

First, we need to extract experimental parameters from the CMap database. This will create a configuration file that lists all available time points, dosages, and cell lines.

```r
library(DRnew)

# Extract parameters and create config file
params <- extract_cmap_parameters(
  siginfo_file = "path/to/siginfo_beta.txt",
  config_dir = "conf",
  config_filename = "cmap_options.conf"
)
```

This will create a file at `conf/cmap_options.conf` with all available parameters. You can edit this file to select specific parameters or keep all of them.

### Step 2: Generate Combinations

After editing the configuration file (if needed), generate combinations of parameters:

```r
# Generate combinations based on the config file
combinations <- generate_combinations_from_config("conf/cmap_options.conf")

# Examine the combinations
head(combinations)
```

### Step 3: Extract Data from CMap

Now you can extract the expression data for these combinations:

```r
# Option 1: Process combinations to create individual files
output_files <- process_combinations(
  combinations = combinations,
  output_dir = "output",
  geneinfo_file = "path/to/geneinfo_beta.txt",
  siginfo_file = "path/to/siginfo_beta.txt",
  gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx"
)

# Option 2: Extract all data into a single data frame
reference_df <- extract_cmap_data_from_config(
  config_file = "conf/cmap_options.conf",
  geneinfo_file = "path/to/geneinfo_beta.txt",
  siginfo_file = "path/to/siginfo_beta.txt",
  gctx_file = "path/to/level5_beta_trt_cp_n720216x12328.gctx"
)
```

### Step 4: Prepare Your Signature File

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

### Step 5: Match Your Signature Against Reference Profiles

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

## Available Signature Matching Methods

DRnew provides access to various signature matching methods:

1. **KS Score** (`ks`): Kolmogorov-Smirnov statistic based method
2. **XCos Score** (`xcos`): Extended cosine similarity method
3. **XSum Score** (`xsum`): Extended sum score method
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

## Example: Batch Processing with SLURM

For large datasets, you can use a high-performance computing cluster with SLURM:

```r
# Save combinations to a file
write.table(combinations, file = "combinations.txt", row.names = FALSE)

# Generate a SLURM submission script
generate_slurm_script(
  combinations_file = "combinations.txt",
  output_file = "run_cmap_extract.sh",
  job_name = "cmap_extract",
  wall_time = "04:00:00",
  memory = "24G"
)

# Submit the job from the command line:
# sbatch run_cmap_extract.sh
```


## License

This package is licensed under the MIT License - see the LICENSE file for details.
