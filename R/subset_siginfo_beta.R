#' Subset siginfo_beta file interactively
#'
#' This function allows interactive exploration and filtering of the siginfo_beta file
#' based on available values in key columns (pert_type, pert_itime, pert_idose, cell_iname).
#'
#' @param siginfo_file Path to siginfo_beta.txt file
#' @param output_file Path to save the filtered siginfo file (optional)
#' @param interactive Logical; whether to run in interactive mode (default: TRUE)
#' @param filters List of pre-defined filters to apply in non-interactive mode
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A filtered data frame of the siginfo file
#'
#' @examples
#' \dontrun{
#' # Interactive mode - will prompt for selections
#' filtered_siginfo <- subset_siginfo_beta("path/to/siginfo_beta.txt")
#'
#' # Non-interactive mode with pre-defined filters
#' filtered_siginfo <- subset_siginfo_beta(
#'   "path/to/siginfo_beta.txt",
#'   interactive = FALSE,
#'   filters = list(
#'     pert_type = "trt_cp",
#'     pert_itime = c("6 h", "24 h"),
#'     pert_idose = "10 uM",
#'     cell_iname = c("A375", "MCF7")
#'   )
#' )
#' }
#'
#' @export
subset_siginfo_beta <- function(siginfo_file,
                                output_file = NULL,
                                interactive = TRUE,
                                filters = NULL,
                                verbose = TRUE) {

  # Check if file exists
  if (!file.exists(siginfo_file)) {
    stop("Siginfo file not found: ", siginfo_file)
  }

  # Read the siginfo file
  if (verbose) message("Reading siginfo file: ", siginfo_file)

  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      sig_info <- data.table::fread(siginfo_file, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE, data.table = FALSE)
    } else {
      sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE, quote = "",
                                    comment.char = "", fill = TRUE)
    }
  }, error = function(e) {
    stop("Error reading siginfo file: ", e$message)
  })

  if (verbose) {
    message(sprintf("Loaded siginfo data with %d rows and %d columns",
                    nrow(sig_info), ncol(sig_info)))
  }

  # Define the key columns to filter on
  key_columns <- c("pert_type", "pert_itime", "pert_idose", "cell_iname")

  # Check if all key columns exist
  missing_cols <- setdiff(key_columns, names(sig_info))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
    key_columns <- intersect(key_columns, names(sig_info))
  }

  # Apply filters
  filtered_data <- sig_info

  if (interactive) {
    # Interactive mode - let user select filters
    for (col in key_columns) {
      # Get unique values
      unique_vals <- sort(unique(filtered_data[[col]]))

      if (length(unique_vals) == 0) {
        next
      }

      # Display available values
      cat("\n", paste(rep("=", 50), collapse = ""), "\n")
      cat("Column:", col, "\n")
      cat("Available values:\n")

      for (i in seq_along(unique_vals)) {
        cat(sprintf("  %d: %s\n", i, unique_vals[i]))
      }

      # Get user selection
      cat("\nSelect values to keep (comma-separated numbers, 'all' for all, 'skip' to skip):\n")
      selection <- readline(prompt = "> ")

      if (tolower(selection) == "skip" || selection == "") {
        next
      } else if (tolower(selection) == "all") {
        # Keep all values
        next
      } else {
        # Parse selection
        indices <- as.integer(unlist(strsplit(selection, ",")))
        indices <- indices[!is.na(indices)]

        # Validate indices
        valid_indices <- indices[indices >= 1 & indices <= length(unique_vals)]

        if (length(valid_indices) > 0) {
          selected_vals <- unique_vals[valid_indices]
          filtered_data <- filtered_data[filtered_data[[col]] %in% selected_vals, ]

          if (verbose) {
            message(sprintf("Filtered %s to %d value(s), %d rows remaining",
                            col, length(selected_vals), nrow(filtered_data)))
          }
        }
      }
    }

  } else {
    # Non-interactive mode - use provided filters
    if (!is.null(filters)) {
      for (col in names(filters)) {
        if (col %in% names(filtered_data)) {
          filtered_data <- filtered_data[filtered_data[[col]] %in% filters[[col]], ]

          if (verbose) {
            message(sprintf("Filtered %s to %d value(s), %d rows remaining",
                            col, length(filters[[col]]), nrow(filtered_data)))
          }
        }
      }
    }
  }

  # Show summary of filtered data
  if (verbose) {
    cat("\n", paste(rep("=", 50), collapse = ""), "\n")
    cat("Filtering complete!\n")
    cat(sprintf("Original rows: %d\n", nrow(sig_info)))
    cat(sprintf("Filtered rows: %d (%.1f%%)\n",
                nrow(filtered_data),
                nrow(filtered_data) / nrow(sig_info) * 100))

    # Show distribution of remaining values
    cat("\nDistribution of remaining values:\n")
    for (col in key_columns) {
      if (col %in% names(filtered_data)) {
        val_counts <- table(filtered_data[[col]])
        cat(sprintf("\n%s (%d unique values):\n", col, length(val_counts)))

        # Show up to 10 values
        show_n <- min(length(val_counts), 10)
        for (i in 1:show_n) {
          cat(sprintf("  %s: %d\n", names(val_counts)[i], val_counts[i]))
        }
        if (length(val_counts) > 10) {
          cat("  ... and", length(val_counts) - 10, "more\n")
        }
      }
    }
  }

  # Save filtered data if output file is specified
  if (!is.null(output_file)) {
    if (verbose) message("\nSaving filtered data to: ", output_file)

    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir) && output_dir != ".") {
      dir.create(output_dir, recursive = TRUE)
    }

    # Write the file
    utils::write.table(filtered_data, file = output_file,
                       sep = "\t", row.names = FALSE, quote = FALSE)

    if (verbose) {
      message(sprintf("Saved %d rows to %s", nrow(filtered_data), output_file))
    }
  }

  return(filtered_data)
}


#' Explore siginfo_beta file structure
#'
#' This function provides a quick overview of the siginfo_beta file structure,
#' showing available columns and sample values.
#'
#' @param siginfo_file Path to siginfo_beta.txt file
#' @param n_examples Number of example values to show for each column (default: 5)
#'
#' @return Invisibly returns the column information
#'
#' @examples
#' \dontrun{
#' explore_siginfo_beta("path/to/siginfo_beta.txt")
#' }
#'
#' @export
explore_siginfo_beta <- function(siginfo_file, n_examples = 5) {

  # Check if file exists
  if (!file.exists(siginfo_file)) {
    stop("Siginfo file not found: ", siginfo_file)
  }

  # Read the first 1000 rows to get a sense of the data
  cat("Reading siginfo file to explore structure...\n")

  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      sig_info <- data.table::fread(siginfo_file, sep = "\t", header = TRUE,
                                    nrows = 1000, stringsAsFactors = FALSE,
                                    data.table = FALSE)
    } else {
      sig_info <- utils::read.table(siginfo_file, sep = "\t", header = TRUE,
                                    nrows = 1000, stringsAsFactors = FALSE,
                                    quote = "", comment.char = "", fill = TRUE)
    }
  }, error = function(e) {
    stop("Error reading siginfo file: ", e$message)
  })

  # Get column names and types
  col_info <- data.frame(
    column = names(sig_info),
    type = sapply(sig_info, class),
    n_unique = sapply(sig_info, function(x) length(unique(x))),
    stringsAsFactors = FALSE
  )

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SIGINFO FILE STRUCTURE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("File: %s\n", siginfo_file))
  cat(sprintf("Columns: %d\n", ncol(sig_info)))
  cat(sprintf("Sample rows: %d\n", nrow(sig_info)))

  # Key columns
  key_cols <- c("pert_type", "pert_itime", "pert_idose", "cell_iname",
                "pert_iname", "is_hiq", "distil_cc_q75")

  cat("\nKEY COLUMNS:\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (col in key_cols) {
    if (col %in% names(sig_info)) {
      unique_vals <- unique(sig_info[[col]])
      n_unique <- length(unique_vals)

      cat(sprintf("\n%s (type: %s, unique values: %d)\n",
                  col, class(sig_info[[col]])[1], n_unique))

      # Show examples
      if (n_unique <= n_examples * 2) {
        # Show all values if not too many
        cat("  Values: ", paste(sort(unique_vals), collapse = ", "), "\n")
      } else {
        # Show first few examples
        examples <- head(sort(unique_vals), n_examples)
        cat("  Examples: ", paste(examples, collapse = ", "), "...\n")
      }
    }
  }

  cat("\nALL COLUMNS:\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  print(col_info)

  invisible(col_info)
}


#' Create a custom configuration file from filtered siginfo
#'
#' After filtering siginfo_beta, this function extracts the remaining
#' time points, doses, and cell lines to create a new configuration file.
#'
#' @param filtered_siginfo Filtered siginfo data frame
#' @param output_file Path to save the configuration file
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return Path to the created configuration file (invisibly)
#'
#' @examples
#' \dontrun{
#' # First filter the siginfo
#' filtered <- subset_siginfo_beta("siginfo_beta.txt",
#'                                filters = list(pert_type = "trt_cp"))
#'
#' # Then create a config from the filtered data
#' create_config_from_filtered_siginfo(filtered, "conf/filtered_config.conf")
#' }
#'
#' @export
create_config_from_filtered_siginfo <- function(filtered_siginfo,
                                                output_file = "conf/filtered_options.conf",
                                                verbose = TRUE) {

  # Extract unique values from filtered data
  times <- sort(unique(filtered_siginfo$pert_itime))
  doses <- sort(unique(filtered_siginfo$pert_idose))
  cells <- sort(unique(filtered_siginfo$cell_iname))

  # Create options list
  options_list <- list(
    times = times,
    doses = doses,
    cells = cells
  )

  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }

  # Use the existing write_config_file function from the package
  write_config_file(options_list, output_file)

  if (verbose) {
    message("\nCreated configuration file from filtered siginfo:")
    message("  File: ", output_file)
    message("  Time points: ", length(times))
    message("  Doses: ", length(doses))
    message("  Cell lines: ", length(cells))
  }

  return(invisible(output_file))
}
