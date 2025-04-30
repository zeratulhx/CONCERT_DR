#' Extract Combined Data from CMap GCTX Based on Config File
#'
#' @description Extract expression data from CMap GCTX files based on parameters
#' specified in a configuration file, returning a combined data frame.
#'
#' @param config_file Path to configuration file with selected parameters
#' @param geneinfo_file Path to the gene info file
#' @param siginfo_file Path to the signature info file
#' @param gctx_file Path to the GCTX file
#' @param keep_all_genes Logical; whether to keep all genes across combinations (TRUE) or only common genes (FALSE) (default: TRUE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A data frame with combined results from all combinations, with annotation
#'         columns indicating time, dose, and cell line for each column
#'
#' @examples
#' \dontrun{
#' # First create and edit a configuration file
#' create_cmap_config_template()
#' # Edit the conf/cmap_options_template.conf file
#'
#' # Then extract data directly based on the config
#' results <- extract_cmap_data_from_config("conf/cmap_options_template.conf")
#'
#' # Write to CSV if needed
#' write.csv(results, file = "cmap_results.csv")
#' }
#'
#' @export
extract_cmap_data_from_config <- function(config_file,
                                          geneinfo_file = "geneinfo_beta.txt",
                                          siginfo_file = "siginfo_beta.txt",
                                          gctx_file = "level5_beta_trt_cp_n720216x12328.gctx",
                                          keep_all_genes = TRUE,
                                          verbose = TRUE) {

  # Check if config file exists
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }

  # Generate combinations from config file
  if (verbose) message("Reading configuration from: ", config_file)
  combinations <- generate_combinations_from_config(config_file, verbose = verbose)

  if (nrow(combinations) == 0) {
    stop("No valid combinations found in config file.")
  }

  # Read gene info file
  if (verbose) message("Reading gene info file...")
  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      geneinfo_df <- data.table::fread(geneinfo_file, sep = "\t", header = TRUE)
    } else {
      # If data.table is not available, use read.table with more robust parameters
      geneinfo_df <- utils::read.table(geneinfo_file, sep = "\t", header = TRUE,
                                       stringsAsFactors = FALSE, quote = "",
                                       comment.char = "", fill = TRUE)
    }
    result <- get_rid(geneinfo_df)
    rid <- result$rid
    genenames <- result$genenames

    if (verbose) message("Found ", length(rid), " landmark genes")
  }, error = function(e) {
    stop("Error reading gene info file: ", e$message)
  })

  # Read siginfo file
  if (verbose) message("Reading signature info file...")
  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      sig_info <- data.table::fread(siginfo_file, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE)
    } else {
      # If data.table is not available, use read.table with more robust parameters
      sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE, quote = "",
                                    comment.char = "", fill = TRUE)
    }
    sig_info <- sig_info[sig_info$pert_type == "trt_cp", ]
    sig_info <- sig_info[sig_info$is_hiq == 1, ]

    if (verbose) message("Found ", nrow(sig_info), " high-quality treatment signatures")
  }, error = function(e) {
    stop("Error reading signature info file: ", e$message)
  })

  # Initialize results list to store all dataframes
  result_list <- list()

  # Process each combination
  for (i in 1:nrow(combinations)) {
    combination <- combinations[i, ]
    itime <- combination$itime
    idose <- combination$idose
    cell <- combination$cell

    if (verbose) {
      message(sprintf("\n--- Combination %d of %d ---", i, nrow(combinations)))
      message(sprintf("Processing: time=%s, dose=%s, cell=%s", itime, idose, cell))
    }

    # Filter data
    filtered_df <- sig_info[
      sig_info$pert_itime == itime &
        sig_info$pert_idose == idose &
        sig_info$cell_iname == cell,
    ]

    cid <- filtered_df$sig_id

    if (length(cid) > 0) {
      # Parse the gctx file using CmapR
      pert_data <- cmapR::parse_gctx(
        fname = gctx_file,
        cid = cid,
        rid = rid
      )

      # Get the data matrix
      mat <- cmapR::mat(pert_data)

      if (verbose) message(sprintf("Data dimensions: %d x %d", nrow(mat), ncol(mat)))

      if (nrow(mat) > 0) {
        # Convert to data frame and set row names as gene names
        df <- as.data.frame(mat)
        rownames(df) <- genenames

        # Add metadata columns for this combination
        metadata <- data.frame(
          sample_id = colnames(df),
          time = itime,
          dose = idose,
          cell = cell,
          stringsAsFactors = FALSE
        )

        # Add to result list
        result_list[[i]] <- list(data = df, metadata = metadata)
      } else {
        if (verbose) message("No data found for this combination")
      }
    } else {
      if (verbose) message("No matching signatures found for this combination")
    }
  }

  # Check if any results were found
  if (length(result_list) == 0) {
    warning("No data found for any of the specified combinations.")
    return(data.frame())
  }

  # Combine all results
  if (verbose) message("\nCombining results...")

  # Initialize variables to hold the combined data
  combined_data <- NULL
  combined_metadata <- NULL

  # Process each result based on the keep_all_genes parameter
  if (keep_all_genes) {
    # Strategy for keeping ALL genes across all combinations

    # First, collect all unique gene names
    all_genes <- character()
    for (result in result_list) {
      if (!is.null(result$data) && ncol(result$data) > 0) {
        all_genes <- union(all_genes, rownames(result$data))
      }
    }

    if (verbose) message("Total unique genes across all combinations: ", length(all_genes))

    # Create the combined matrix with all genes and all samples
    total_samples <- sum(sapply(result_list, function(r) {
      if (!is.null(r$data)) ncol(r$data) else 0
    }))

    # Initialize combined data with NAs
    combined_matrix <- matrix(NA, nrow = length(all_genes), ncol = total_samples)
    rownames(combined_matrix) <- all_genes

    # Collect all sample IDs
    all_sample_ids <- character()

    # Fill in the data for each combination
    col_start <- 1
    for (result in result_list) {
      if (!is.null(result$data) && ncol(result$data) > 0) {
        data_df <- result$data
        num_cols <- ncol(data_df)

        # Add sample IDs
        all_sample_ids <- c(all_sample_ids, colnames(data_df))

        # Fill in data for current combination
        for (gene in rownames(data_df)) {
          idx <- match(gene, all_genes)
          if (!is.na(idx)) {
            combined_matrix[idx, col_start:(col_start + num_cols - 1)] <-
              as.numeric(data_df[gene, ])
          }
        }

        # Combine metadata
        if (is.null(combined_metadata)) {
          combined_metadata <- result$metadata
        } else {
          combined_metadata <- rbind(combined_metadata, result$metadata)
        }

        # Update column start for next combination
        col_start <- col_start + num_cols
      }
    }

    # Set column names
    colnames(combined_matrix) <- all_sample_ids

    # Convert to data frame
    combined_data <- as.data.frame(combined_matrix)

  } else {
    # Original behavior - only keep genes common to all combinations

    for (result in result_list) {
      if (!is.null(result$data) && ncol(result$data) > 0) {
        if (is.null(combined_data)) {
          combined_data <- result$data
          combined_metadata <- result$metadata
        } else {
          # Find common genes
          common_genes <- intersect(rownames(combined_data), rownames(result$data))
          if (verbose) message("Common genes with previous combinations: ", length(common_genes))

          # Append columns to the existing dataframe, keeping only common genes
          combined_data <- cbind(combined_data[common_genes, ], result$data[common_genes, ])
          combined_metadata <- rbind(combined_metadata, result$metadata)
        }
      }
    }
  }

  # Final check
  if (is.null(combined_data)) {
    warning("No valid data found to combine.")
    return(data.frame())
  }

  # Create result metadata
  if (verbose) {
    message(sprintf("Final combined data dimensions: %d genes × %d samples",
                    nrow(combined_data), ncol(combined_data)))
  }

  # Add gene symbols as a column
  final_result <- data.frame(gene_symbol = rownames(combined_data), combined_data,
                             check.names = FALSE, stringsAsFactors = FALSE)

  # Add metadata as attributes
  attr(final_result, "metadata") <- combined_metadata

  return(final_result)
}
