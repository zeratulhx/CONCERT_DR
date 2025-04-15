#' Signature Matching Utility
#'
#' @description Functions to match drug response signatures against reference profiles
#' from the CMap database that were previously extracted using the CMap extract functions.
#'
#' @name signature_matching
NULL

#' Process drug response signatures against reference data
#'
#' @param signature_file Path to signature gene list with log2FC values
#' @param reference_file Path to reference data extracted from CMap
#' @param output_prefix Prefix for output files
#' @param permutations Number of permutations for statistical testing
#' @param methods Vector of method names to run (default: all)
#'        Options: "ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"
#'
#' @return List containing results from selected scoring methods
#' @export
process_signature <- function(signature_file, reference_file, output_prefix = "results",
                              permutations = 100, methods = c("ks", "xcos", "xsum", "gsea0",
                                                              "gsea1", "gsea2", "zhang")) {
  # Ensure RCSM is installed
  if (!requireNamespace("RCSM", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("This function requires 'devtools' and 'RCSM'. Please install them first:
           install.packages('devtools')
           devtools::install_github('Jasonlinchina/RCSM')")
    }
    message("Installing RCSM package from GitHub...")
    devtools::install_github("Jasonlinchina/RCSM")
  }
  
  # Load RCSM
  requireNamespace("RCSM", quietly = TRUE)
  
  # Read the reference data (output from CMap_split_gctx.R)
  message("Reading reference data from ", reference_file)
  ref_data <- utils::read.delim(reference_file, sep = "\t", row.names = 1, check.names = FALSE)

  # Convert to matrix format required by RCSM
  ref <- as.matrix(ref_data)

  # Read the signature file
  message("Reading signature data from ", signature_file)
  gene_data <- utils::read.delim(signature_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Check if required columns exist
  if (!all(c("Gene", "log2FC") %in% colnames(gene_data))) {
    stop("Signature file must contain 'Gene' and 'log2FC' columns")
  }

  # Separate up and down regulated genes based on log2FC values
  Up <- gene_data$Gene[gene_data$log2FC > 0]
  Down <- gene_data$Gene[gene_data$log2FC < 0]

  # Define query (all genes with their ranks)
  query_genes <- gene_data$Gene
  ranks <- rank(-gene_data$log2FC)  # Rank by decreasing log2FC
  query <- stats::setNames(ranks, query_genes)

  # Verify genes exist in reference
  common_up <- intersect(Up, rownames(ref))
  common_down <- intersect(Down, rownames(ref))

  message(sprintf("Found %d/%d up-regulated genes and %d/%d down-regulated genes in reference",
                  length(common_up), length(Up),
                  length(common_down), length(Down)))

  if (length(common_up) == 0 || length(common_down) == 0) {
    warning("No matching genes found in reference data. Check gene identifiers.")
    return(NULL)
  }

  # Available methods
  all_methods <- list(
    ks = function() {
      message("Running KS score...")
      RCSM::KSScore(refMatrix = ref, queryUp = common_up, queryDown = common_down, permuteNum = permutations)
    },

    xcos = function() {
      message("Running XCos score...")
      RCSM::XCosScore(refMatrix = ref, query = query[query_genes %in% rownames(ref)], topN = 4, permuteNum = permutations)
    },

    xsum = function() {
      message("Running XSum score...")
      RCSM::XSumScore(refMatrix = ref, queryUp = common_up, queryDown = common_down, topN = 4, permuteNum = permutations)
    },

    gsea0 = function() {
      message("Running GSEA weight 0 score...")
      RCSM::GSEAweight0Score(refMatrix = ref, queryUp = common_up, queryDown = common_down, permuteNum = permutations)
    },

    gsea1 = function() {
      message("Running GSEA weight 1 score...")
      RCSM::GSEAweight1Score(refMatrix = ref, queryUp = common_up, queryDown = common_down, permuteNum = permutations)
    },

    gsea2 = function() {
      message("Running GSEA weight 2 score...")
      RCSM::GSEAweight2Score(refMatrix = ref, queryUp = common_up, queryDown = common_down, permuteNum = permutations)
    },

    zhang = function() {
      message("Running Zhang score...")
      RCSM::ZhangScore(refMatrix = ref, queryUp = common_up, queryDown = common_down, permuteNum = permutations)
    }
  )

  # Validate selected methods
  invalid_methods <- setdiff(methods, names(all_methods))
  if (length(invalid_methods) > 0) {
    warning("Invalid methods specified: ", paste(invalid_methods, collapse = ", "),
            ". Will be ignored.")
    methods <- intersect(methods, names(all_methods))
  }

  # Run selected scoring methods
  results <- list()

  for (method in methods) {
    tryCatch({
      results[[method]] <- all_methods[[method]]()

      # Save individual results
      output_file <- paste0(output_prefix, "_", method, "_results.csv")
      utils::write.csv(results[[method]], file = output_file, row.names = FALSE)
      message("Saved ", method, " results to ", output_file)
    }, error = function(e) {
      warning("Error running ", method, " method: ", e$message)
    })
  }

  return(results)
}

#' Find reference files based on combinations.txt from CMap_generate_comb.R
#'
#' @param combinations_file Path to combinations.txt file
#' @return Data frame with reference file paths and condition info
#' @export
find_reference_files <- function(combinations_file) {
  # Check if combinations file exists
  if (!file.exists(combinations_file)) {
    stop("Combinations file not found: ", combinations_file)
  }

  # Read combinations file
  message("Reading combinations from ", combinations_file)
  combinations <- utils::read.table(combinations_file, header = TRUE, stringsAsFactors = FALSE)

  # Create a data frame to store file paths and condition info
  ref_files <- data.frame(
    time = combinations$itime,
    dose = combinations$idose,
    cell = combinations$cell,
    file_path = character(nrow(combinations)),
    exists = logical(nrow(combinations)),
    stringsAsFactors = FALSE
  )

  # For each combination, construct expected file path
  for (i in 1:nrow(ref_files)) {
    # Format the filename following the pattern from CMap_split_gctx.R
    filename <- paste0(
      "filtered_",
      gsub(" ", "_", ref_files$time[i]), "_",
      gsub(" ", "_", ref_files$dose[i]), "_",
      ref_files$cell[i], ".csv"
    )

    ref_files$file_path[i] <- filename
    ref_files$exists[i] <- file.exists(filename)
  }

  # Report on found files
  found_files <- sum(ref_files$exists)
  message(sprintf("Found %d/%d expected reference files", found_files, nrow(ref_files)))

  if (found_files == 0) {
    stop("No reference files found. Make sure to run process_combinations_file() first.")
  }

  # Return only rows with existing files
  return(ref_files[ref_files$exists, ])
}

#' Create a summary of top hits across all datasets
#'
#' @param all_results List of results from all reference datasets
#' @param output_file Path to output summary file
#' @return Invisibly returns the summary data frame
#' @keywords internal
create_summary <- function(all_results, output_file) {
  # Initialize summary data frame
  summary_data <- data.frame()

  # For each condition and method, get top compounds
  for (condition in names(all_results)) {
    results <- all_results[[condition]]

    if (is.null(results)) next

    for (method in names(results)) {
      # Get top 10 compounds or fewer if less available
      method_results <- results[[method]]

      if (nrow(method_results) == 0) next

      # Sort by score (assuming higher is better, adjust if needed)
      method_results <- method_results[order(method_results$Score, decreasing = TRUE), ]
      top_n <- min(10, nrow(method_results))

      # Add to summary
      for (i in 1:top_n) {
        if (i <= nrow(method_results)) {
          row <- method_results[i, ]
          summary_data <- rbind(summary_data, data.frame(
            Condition = condition,
            Method = method,
            Rank = i,
            Compound = as.character(row$Name),  # Adjust column name if different
            Score = row$Score,
            PValue = if ("P.Value" %in% colnames(row)) row$P.Value else NA
          ))
        }
      }
    }
  }

  # Save summary
  utils::write.csv(summary_data, file = output_file, row.names = FALSE)
  message("Summary of top hits saved to ", output_file)
  
  return(invisible(summary_data))
}

#' Run analysis with combinations file
#'
#' @param combinations_file Path to combinations.txt file (default: "combinations.txt")
#' @param sig_file Path to signature file (default: "signature.txt")
#' @param out_dir Output directory for results (default: "results")
#' @param methods Vector of method names to run (default: all)
#'        Options: "ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"
#' @param permutations Number of permutations for statistical testing
#' @return List of results from all analyses
#' @export
run_analysis_with_combinations <- function(combinations_file = "combinations.txt",
                                           sig_file = "signature.txt",
                                           out_dir = "results",
                                           methods = c("ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"),
                                           permutations = 100) {
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Find reference files based on combinations
  ref_files_df <- find_reference_files(combinations_file)

  # Process each reference file
  all_results <- list()

  message("Selected methods: ", paste(methods, collapse = ", "))

  for (i in 1:nrow(ref_files_df)) {
    ref_file <- ref_files_df$file_path[i]

    # Create a condition string from time, dose, and cell
    condition <- paste(
      gsub(" ", "_", ref_files_df$time[i]),
      gsub(" ", "_", ref_files_df$dose[i]),
      ref_files_df$cell[i],
      sep = "_"
    )

    message("\nProcessing reference file: ", ref_file)
    message("Condition: ", condition)

    # Create output prefix
    output_prefix <- file.path(out_dir, paste0("sig_match_", condition))

    # Process signature
    results <- process_signature(
      signature_file = sig_file,
      reference_file = ref_file,
      output_prefix = output_prefix,
      permutations = permutations,
      methods = methods
    )

    all_results[[condition]] <- results
  }

  # Create summary of top hits across all reference datasets
  create_summary(all_results, file.path(out_dir, "summary_results.csv"))

  return(all_results)
}

#' Run analysis with manually specified reference files
#'
#' @param ref_pattern Pattern to match reference files (default: "filtered_*.csv")
#' @param sig_file Path to signature file (default: "signature.txt")
#' @param out_dir Output directory for results (default: "results")
#' @param methods Vector of method names to run (default: all)
#'        Options: "ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"
#' @param permutations Number of permutations for statistical testing
#' @return List of results from all analyses
#' @export
run_analysis <- function(ref_pattern = "filtered_*.csv",
                         sig_file = "signature.txt",
                         out_dir = "results",
                         methods = c("ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"),
                         permutations = 100) {
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Get all reference files
  ref_files <- list.files(pattern = ref_pattern, full.names = TRUE)

  if (length(ref_files) == 0) {
    stop("No reference files found matching pattern: ", ref_pattern)
  }

  message("Found ", length(ref_files), " reference files")
  message("Selected methods: ", paste(methods, collapse = ", "))

  # Process each reference file
  all_results <- list()

  for (ref_file in ref_files) {
    message("\nProcessing reference file: ", ref_file)

    # Extract conditions from filename (e.g., filtered_6_h_10_uM_A375.csv)
    file_base <- tools::file_path_sans_ext(basename(ref_file))
    conditions <- gsub("filtered_", "", file_base)

    # Create output prefix
    output_prefix <- file.path(out_dir, paste0("sig_match_", conditions))

    # Process signature
    results <- process_signature(
      signature_file = sig_file,
      reference_file = ref_file,
      output_prefix = output_prefix,
      permutations = permutations,
      methods = methods
    )

    all_results[[conditions]] <- results
  }

  # Create summary of top hits across all reference datasets
  create_summary(all_results, file.path(out_dir, "summary_results.csv"))

  return(all_results)
}