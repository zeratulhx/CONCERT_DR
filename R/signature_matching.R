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
#' @param save_files Logical; whether to save results to files (default: TRUE)
#'
#' @return A structured list containing:
#'   \item{results}{List containing results for each method}
#'   \item{summary}{Data frame of top hits across all methods}
#'   \item{gene_data}{Original signature gene data}
#'   \item{settings}{List of parameters used for the analysis}
#'   \item{common_genes}{Counts of genes found in reference data}
#'   Use print() or summary() methods to get a quick overview of results
#'
#' @export
process_signature_with_df <- function(signature_file, reference_df, output_dir = "results",
                                      permutations = 100, methods = c("ks", "xcos", "xsum", "gsea0",
                                                                      "gsea1", "gsea2", "zhang"),
                                      topN = 4, read_method = "auto", save_files = TRUE) {
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

  # Create output directory if it doesn't exist and files will be saved
  if (save_files && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Record start time for performance tracking
  start_time <- Sys.time()

  # Read the signature file
  message("Reading signature data from ", signature_file)
  tryCatch({
    if (read_method == "auto") {
      # Auto-detect best method
      if (requireNamespace("data.table", quietly = TRUE)) {
        gene_data <- data.table::fread(signature_file, header = TRUE,
                                       stringsAsFactors = FALSE, data.table = FALSE)
      } else {
        gene_data <- utils::read.delim(signature_file, header = TRUE, sep = "\t",
                                       stringsAsFactors = FALSE, comment.char = "")
      }
    } else if (read_method == "fread" && requireNamespace("data.table", quietly = TRUE)) {
      gene_data <- data.table::fread(signature_file, header = TRUE,
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

  # Calculate percent overlap to assess gene coverage
  pct_up <- round(length(common_up) / length(Up) * 100, 1)
  pct_down <- round(length(common_down) / length(Down) * 100, 1)

  message(sprintf("Found %d/%d up-regulated genes (%g%%) and %d/%d down-regulated genes (%g%%) in reference",
                  length(common_up), length(Up), pct_up,
                  length(common_down), length(Down), pct_down))

  if (length(common_up) == 0 || length(common_down) == 0) {
    stop("No matching genes found in reference data. Please check your signature genes.")
  }

  # Initialize results list
  all_results <- list()

  # Store metadata about the analysis
  settings <- list(
    signature_file = signature_file,
    permutations = permutations,
    methods = methods,
    topN = topN,
    time_started = start_time,
    reference_dimensions = dim(ref),
    reference_genes = length(rownames(ref))
  )

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
      result_df <- all_methods[[method]]()

      # Add compound names if not already included
      if (!("compound" %in% colnames(result_df))) {
        result_df <- cbind(compound = rownames(result_df), result_df)
      }

      # Add rank column for easier interpretation
      result_df$rank <- rank(-result_df$Score)

      # Store in results list
      all_results[[method]] <- result_df

      # Save individual results if requested
      if (save_files) {
        output_file <- file.path(output_dir, paste0("sig_match_", method, "_results.csv"))
        utils::write.csv(result_df, file = output_file, row.names = FALSE)
        message("Saved ", method, " results to ", output_file)
      }
    }, error = function(e) {
      warning("Error running ", method, " method: ", e$message)
      all_results[[method]] <- data.frame(
        compound = character(0),
        Score = numeric(0),
        pValue = numeric(0),
        error = character(0)
      )
      all_results[[method]]$error <- e$message
    })
  }

  # Create summary of top hits across all methods
  summary_df <- create_summary_from_results(all_results)

  # Save summary file if requested
  if (save_files) {
    summary_file <- file.path(output_dir, "summary_results.csv")
    utils::write.csv(summary_df, file = summary_file, row.names = FALSE)
    message("Saved summary report to ", summary_file)
  }

  # Record end time
  end_time <- Sys.time()
  time_taken <- difftime(end_time, start_time, units = "mins")

  # Update settings with completion info
  settings$time_completed <- end_time
  settings$time_taken_mins <- as.numeric(time_taken)

  # Create a structured result object
  result_object <- structure(
    list(
      results = all_results,
      summary = summary_df,
      gene_data = gene_data,
      settings = settings,
      common_genes = list(
        up = list(found = common_up, count = length(common_up), percent = pct_up),
        down = list(found = common_down, count = length(common_down), percent = pct_down)
      )
    ),
    class = "cmap_signature_result"
  )

  # Return the complete result object
  return(result_object)
}

#' Create a summary from signature matching results
#'
#' @param results_list List of results from different methods
#' @param top_n Number of top hits to include (default: 20)
#'
#' @return Data frame with summary of top hits across methods
#'
#' @keywords internal
create_summary_from_results <- function(results_list, top_n = 20) {
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

    # Skip if error or empty
    if (nrow(result) == 0 || "error" %in% colnames(result)) {
      next
    }

    # Make sure rank column exists
    if (!"rank" %in% colnames(result)) {
      result$rank <- rank(-result$Score)
    }

    # Get top hits
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

  # Add global rank across methods
  if (nrow(summary_df) > 0) {
    # Calculate a weighted score considering both Score and p-value
    summary_df$weighted_score <- summary_df$Score * (1 - summary_df$pValue)

    # Rank across all methods
    summary_df$global_rank <- rank(-summary_df$weighted_score)

    # Remove the temporary weighted score column
    summary_df$weighted_score <- NULL
  }

  return(summary_df)
}

#' Print method for cmap_signature_result objects
#'
#' @param x A cmap_signature_result object
#' @param ... Additional arguments (not used)
#'
#' @return x invisibly
#'
#' @export
print.cmap_signature_result <- function(x, ...) {
  cat("CMap Signature Matching Results\n")
  cat("===============================\n\n")

  # Print basic information
  cat(sprintf("Signature file: %s\n", x$settings$signature_file))
  cat(sprintf("Analysis completed: %s\n", format(x$settings$time_completed)))
  cat(sprintf("Time taken: %.2f minutes\n\n", x$settings$time_taken_mins))

  # Print methods used
  cat("Methods used:\n")
  for (method in names(x$results)) {
    result <- x$results[[method]]
    if ("error" %in% colnames(result)) {
      cat(sprintf("  - %s: ERROR - %s\n", method, result$error[1]))
    } else {
      cat(sprintf("  - %s: %d results\n", method, nrow(result)))
    }
  }

  # Print gene coverage
  cat("\nGene coverage:\n")
  cat(sprintf("  Up-regulated: %d/%d genes (%.1f%%)\n",
              x$common_genes$up$count, length(x$gene_data$Gene[x$gene_data$log2FC > 0]),
              x$common_genes$up$percent))
  cat(sprintf("  Down-regulated: %d/%d genes (%.1f%%)\n",
              x$common_genes$down$count, length(x$gene_data$Gene[x$gene_data$log2FC < 0]),
              x$common_genes$down$percent))

  # Print top compounds from summary (if available)
  if (nrow(x$summary) > 0) {
    top_n <- min(5, nrow(x$summary))
    top_compounds <- x$summary[order(x$summary$global_rank)[1:top_n], ]

    cat("\nTop compounds across all methods:\n")
    for (i in 1:nrow(top_compounds)) {
      cat(sprintf("  %d. %s (method: %s, score: %.4f, p-value: %.4f)\n",
                  i, top_compounds$compound[i], top_compounds$method[i],
                  top_compounds$Score[i], top_compounds$pValue[i]))
    }
  }

  cat("\nUse the following to explore the results:\n")
  cat("  - $results: List of result data frames for each method\n")
  cat("  - $summary: Summary of top hits across all methods\n")
  cat("  - $gene_data: Original signature gene data\n")
  cat("  - $settings: Analysis settings and metadata\n")
  cat("  - $common_genes: Genes found in the reference data\n")
  cat("\nExample: result_obj$results$ks to view KS score results\n")

  invisible(x)
}

#' Summary method for cmap_signature_result objects
#'
#' @param object A cmap_signature_result object
#' @param top_n Number of top compounds to show (default: 10)
#' @param ... Additional arguments (not used)
#'
#' @return A data frame with the top compounds
#'
#' @export
summary.cmap_signature_result <- function(object, top_n = 10, ...) {
  if (nrow(object$summary) == 0) {
    message("No summary data available.")
    return(NULL)
  }

  # Get top compounds across all methods
  top_compounds <- object$summary[order(object$summary$global_rank), ]

  # Limit to requested number
  top_compounds <- head(top_compounds, top_n)

  # Print a summary
  cat("Top", top_n, "compounds across all methods:\n\n")
  print(top_compounds[, c("compound", "method", "Score", "pValue", "global_rank")])

  # Return the data invisibly
  invisible(top_compounds)
}

#' Plot method for cmap_signature_result objects
#'
#' @param x A cmap_signature_result object
#' @param method Which method to plot (default: first method in results)
#' @param plot_type Type of plot: "scores", "volcano", or "heatmap" (default: "scores")
#' @param top_n Number of compounds to include (default: 20)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot object
#'
#' @export
plot.cmap_signature_result <- function(x, method = NULL, plot_type = "scores", top_n = 20, ...) {
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it with: install.packages('ggplot2')")
  }

  # If method is not specified, use the first available method
  if (is.null(method)) {
    if (length(x$results) == 0) {
      stop("No results available to plot.")
    }
    method <- names(x$results)[1]
  }

  # Check if the specified method exists
  if (!method %in% names(x$results)) {
    stop("Method '", method, "' not found in results. Available methods: ",
         paste(names(x$results), collapse = ", "))
  }

  # Get the data for the specified method
  result_df <- x$results[[method]]

  # Check if there was an error with this method
  if ("error" %in% colnames(result_df)) {
    stop("Cannot plot results for method '", method, "' due to error: ", result_df$error[1])
  }

  # Different plot types
  if (plot_type == "scores") {
    # Bar plot of top scores
    top_df <- result_df[order(-result_df$Score), ][1:min(top_n, nrow(result_df)), ]

    p <- ggplot2::ggplot(top_df, ggplot2::aes(x = reorder(compound, Score), y = Score)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste("Top", top_n, "Compounds by", method, "Score"),
        x = "Compound",
        y = paste(method, "Score")
      ) +
      ggplot2::theme_minimal()

  } else if (plot_type == "volcano") {
    # Volcano plot of scores vs p-values
    p <- ggplot2::ggplot(result_df, ggplot2::aes(x = Score, y = -log10(pValue))) +
      ggplot2::geom_point(ggplot2::aes(color = abs(Score) > 0.5 & pValue < 0.05)) +
      ggplot2::scale_color_manual(values = c("grey", "red"),
                                  name = "Significant",
                                  labels = c("No", "Yes")) +
      ggplot2::labs(
        title = paste("Volcano Plot of", method, "Results"),
        x = paste(method, "Score"),
        y = "-log10(p-value)"
      ) +
      ggplot2::theme_minimal()

  } else if (plot_type == "heatmap") {
    # Check if pheatmap is available
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      stop("Package 'pheatmap' is required for heatmap plotting. Please install it with: install.packages('pheatmap')")
    }

    # Get top compounds
    top_df <- result_df[order(-result_df$Score), ][1:min(top_n, nrow(result_df)), ]

    # Create a heatmap of gene expression for top compounds
    # This would need the actual expression data, not just the scores
    # We'll show a message and return NULL for now
    message("Heatmap plotting requires the actual gene expression data for the top compounds.")
    message("This feature is not fully implemented in this version.")
    return(NULL)

  } else {
    stop("Invalid plot_type: '", plot_type, "'. Must be one of: 'scores', 'volcano', 'heatmap'")
  }

  return(p)
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
