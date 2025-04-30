#' Process drug response signatures against reference data in a dataframe
#'
#' @param signature_file Path to signature gene list with log2FC values
#' @param reference_df Dataframe containing reference data (combined across all conditions)
#' @param output_dir Directory for output files (default: "results")
#' @param permutations Number of permutations for statistical testing (default: 100)
#' @param methods Vector of method names to run (default: all)
#'        Options: "ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"
#' @param topN Integer; number of top-ranked genes to use for XCos and XSum methods (default: 4)
#' @param read_method Character; method to use for reading signature file ("auto", "fread", or "read.table") (default: "auto")
#'
#' @return List containing results for each method
#' @export
process_signature_with_df <- function(signature_file, reference_df, output_dir = "results",
                                      permutations = 100, methods = c("ks", "xcos", "xsum", "gsea0",
                                                                      "gsea1", "gsea2", "zhang"),
                                      topN = 4, read_method = "auto") {
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

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Read the signature file
  message("Reading signature data from ", signature_file)
  tryCatch({
    if (read_method == "auto") {
      # Auto-detect best method
      if (requireNamespace("data.table", quietly = TRUE)) {
        gene_data <- data.table::fread(signature_file, sep = "\t", header = TRUE,
                                       stringsAsFactors = FALSE, data.table = FALSE)
      } else {
        gene_data <- utils::read.delim(signature_file, header = TRUE, sep = "\t",
                                       stringsAsFactors = FALSE, comment.char = "")
      }
    } else if (read_method == "fread" && requireNamespace("data.table", quietly = TRUE)) {
      gene_data <- data.table::fread(signature_file, sep = "\t", header = TRUE,
                                     stringsAsFactors = FALSE, data.table = FALSE)
    } else {
      # Default to read.table
      gene_data <- utils::read.delim(signature_file, header = TRUE, sep = "\t",
                                     stringsAsFactors = FALSE, comment.char = "")
    }
  }, error = function(e) {
    stop("Error reading signature file: ", e$message)
  })

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

  # Prepare reference data for processing
  message("Preparing reference data for analysis...")

  # Check if reference_df has gene_symbol column
  if (!("gene_symbol" %in% colnames(reference_df))) {
    stop("Reference dataframe must contain a 'gene_symbol' column")
  }

  # Create reference matrix with genes as row names
  ref_data <- reference_df
  rownames(ref_data) <- ref_data$gene_symbol
  ref_data$gene_symbol <- NULL

  # Convert to matrix format required by RCSM
  ref <- as.matrix(ref_data)

  # Verify genes exist in reference
  common_up <- intersect(Up, rownames(ref))
  common_down <- intersect(Down, rownames(ref))

  message(sprintf("Found %d/%d up-regulated genes and %d/%d down-regulated genes in reference",
                  length(common_up), length(Up),
                  length(common_down), length(Down)))

  if (length(common_up) == 0 || length(common_down) == 0) {
    stop("No matching genes found in reference data. Please check your signature genes.")
  }

  # Initialize results list
  all_results <- list()

  # Available methods
  all_methods <- list(
    ks = function() {
      message("Running KS score...")
      RCSM::KSScore(refMatrix = ref, queryUp = common_up, queryDown = common_down, permuteNum = permutations)
    },

    xcos = function() {
      message(sprintf("Running XCos score with topN = %d...", topN))
      RCSM::XCosScore(refMatrix = ref, query = query[query_genes %in% rownames(ref)],
                      topN = topN, permuteNum = permutations)
    },

    xsum = function() {
      message(sprintf("Running XSum score with topN = %d...", topN))
      RCSM::XSumScore(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                      topN = topN, permuteNum = permutations)
    },

    gsea0 = function() {
      message("Running GSEA weight 0 score...")
      RCSM::GSEAweight0Score(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                             permuteNum = permutations)
    },

    gsea1 = function() {
      message("Running GSEA weight 1 score...")
      RCSM::GSEAweight1Score(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                             permuteNum = permutations)
    },

    gsea2 = function() {
      message("Running GSEA weight 2 score...")
      RCSM::GSEAweight2Score(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                             permuteNum = permutations)
    },

    zhang = function() {
      message("Running Zhang score...")
      RCSM::ZhangScore(refMatrix = ref, queryUp = common_up, queryDown = common_down,
                       permuteNum = permutations)
    }
  )

  # Validate selected methods
  invalid_methods <- setdiff(methods, names(all_methods))
  if (length(invalid_methods) > 0) {
    message("Invalid methods specified: ", paste(invalid_methods, collapse = ", "),
            ". Will be ignored.")
    methods <- intersect(methods, names(all_methods))
  }

  # Run selected scoring methods
  for (method in methods) {
    tryCatch({
      all_results[[method]] <- all_methods[[method]]()

      # Add compound names if not already included
      if (!("compound" %in% colnames(all_results[[method]]))) {
        all_results[[method]] <- cbind(compound = rownames(all_results[[method]]),
                                       all_results[[method]])
      }

      # Save individual results
      output_file <- file.path(output_dir, paste0("sig_match_", method, "_results.csv"))
      utils::write.csv(all_results[[method]], file = output_file, row.names = FALSE)
      message("Saved ", method, " results to ", output_file)
    }, error = function(e) {
      warning("Error running ", method, " method: ", e$message)
    })
  }

  # Create summary file with top hits across all methods
  summary_file <- file.path(output_dir, "summary_results.csv")
  create_summary_report(all_results, summary_file)
  message("Saved summary report to ", summary_file)

  return(all_results)
}
#' Create a summary report of top hits across all methods
#'
#' @param results_list List of results from different methods
#' @param output_file Path to output summary file
#' @param top_n Number of top hits to include (default: 20)
#'
#' @return Invisible NULL
#' @keywords internal
create_summary_report <- function(results_list, output_file, top_n = 20) {
  # Initialize summary dataframe
  summary_df <- data.frame(
    compound = character(),
    method = character(),
    Score = numeric(),
    pValue = numeric(),
    rank = integer(),
    stringsAsFactors = FALSE
  )

  # Extract top hits from each method
  for (method_name in names(results_list)) {
    result <- results_list[[method_name]]


    # Get top hits
    result$compound <- rownames(result)
    result$rank <- rank(-result$Score)
    top_hits <- result[result$rank <= top_n, ]

    if (nrow(top_hits) > 0) {
      top_hits$method <- method_name
      summary_df <- rbind(
        summary_df,
        top_hits[, c("compound", "method", "Score", "pValue", "rank")]
      )
    }
  }

  # Sort by method and rank
  summary_df <- summary_df[order(summary_df$method, summary_df$rank), ]

  # Write to file
  utils::write.csv(summary_df, file = output_file, row.names = FALSE)

  return(invisible(NULL))
}

#' Complete workflow from configuration to signature matching
#'
#' @param config_file Path to configuration file with selected parameters
#' @param signature_file Path to signature gene list with log2FC values
#' @param geneinfo_file Path to the gene info file
#' @param siginfo_file Path to the signature info file
#' @param gctx_file Path to the GCTX file
#' @param output_dir Directory for output files (default: "results")
#' @param methods Vector of method names to run (default: all)
#'        Options: "ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"
#' @param topN Integer; number of top-ranked genes to use for XCos and XSum methods (default: 4)
#' @param permutations Number of permutations for statistical testing (default: 100)
#' @param keep_all_genes Logical; whether to keep all genes when extracting CMap data (default: TRUE)
#' @param read_method Character; method to use for reading signature file ("auto", "fread", or "read.table") (default: "auto")
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return List containing results from all methods
#' @export
run_cmap_workflow <- function(config_file, signature_file,
                              geneinfo_file = "geneinfo_beta.txt",
                              siginfo_file = "siginfo_beta.txt",
                              gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                              output_dir = "results",
                              methods = c("ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang"),
                              topN = 4,
                              permutations = 100,
                              keep_all_genes = TRUE,
                              read_method = "auto",
                              verbose = TRUE) {

  # Step 1: Extract CMap data based on config file
  if (verbose) message("Step 1: Extracting CMap data based on configuration...")
  reference_df <- extract_cmap_data_from_config(
    config_file = config_file,
    geneinfo_file = geneinfo_file,
    siginfo_file = siginfo_file,
    gctx_file = gctx_file,
    keep_all_genes = keep_all_genes,
    verbose = verbose
  )

  if (nrow(reference_df) == 0) {
    stop("No data extracted from CMap. Please check your configuration and input files.")
  }

  # Step 2: Process signature against the extracted data
  if (verbose) message("\nStep 2: Processing signature against reference data...")
  results <- process_signature_with_df(
    signature_file = signature_file,
    reference_df = reference_df,
    output_dir = output_dir,
    permutations = permutations,
    methods = methods,
    topN = topN,
    read_method = read_method
  )

  if (verbose) message("\nWorkflow completed successfully!")
  return(results)
}
