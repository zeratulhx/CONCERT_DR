#' CMap Utility Functions
#'
#' @description Utility functions for interactive selection of parameters and other
#' CMap workflow helpers.
#'
#' @name cmap_utils
NULL

#' Select items from a list interactively
#'
#' @param all_items Character vector of all available items
#' @param item_type String describing the type of items
#' @return Character vector of selected items
#' @export
select_items_interactive <- function(all_items, item_type) {
  cat(sprintf("\nAvailable %s:\n", item_type))
  for (i in seq_along(all_items)) {
    cat(sprintf("%d: %s\n", i, all_items[i]))
  }

  cat(sprintf("\nEnter numbers for %s to include (comma-separated, or 'a' for all): ", item_type))
  selection <- readline(prompt = "")

  if (tolower(selection) == "a") {
    return(all_items)
  } else {
    # Handle empty input
    if (selection == "") {
      message("No selection made. Using first item.")
      return(all_items[1])
    }
    
    # Parse selection
    indices <- as.integer(unlist(strsplit(selection, ",")))
    
    # Validate indices
    valid_indices <- indices[indices >= 1 & indices <= length(all_items)]
    
    if (length(valid_indices) == 0) {
      message("No valid selection made. Using first item.")
      return(all_items[1])
    }
    
    return(all_items[valid_indices])
  }
}

#' Get input with a default value
#'
#' @param prompt Text prompt to display
#' @param default Default value to use if no input is provided
#' @return User input or default value
#' @export
get_input_with_default <- function(prompt, default) {
  cat(sprintf("%s [%s]: ", prompt, default))
  input <- readline(prompt = "")
  if (input == "") {
    return(default)
  } else {
    return(input)
  }
}

#' Generate a SLURM submission script
#'
#' @param combinations_file Path to the combinations file
#' @param output_file Path to save the SLURM script
#' @param job_name Job name for SLURM (default: "cmap_analysis")
#' @param wall_time Wall time limit (default: "06:00:00")
#' @param memory Memory per node (default: "24G")
#' @param cpus_per_task CPUs per task (default: 2)
#' @param make_executable Whether to make the script executable (default: TRUE)
#' @return Path to the generated script (invisibly)
#' @export
generate_slurm_script <- function(combinations_file, output_file = "run_analysis.sh",
                                 job_name = "cmap_analysis", wall_time = "06:00:00",
                                 memory = "24G", cpus_per_task = 2,
                                 make_executable = TRUE) {
  
  # Check if combinations file exists
  if (!file.exists(combinations_file)) {
    stop("Combinations file not found: ", combinations_file)
  }
  
  # Count the number of combinations (lines minus header)
  total_combinations <- length(readLines(combinations_file)) - 1
  
  # Create SLURM script
  slurm_script <- c(
    "#!/bin/bash",
    paste0("#SBATCH --job-name=", job_name),
    paste0("#SBATCH --output=", job_name, "_%A_%a.out"),
    paste0("#SBATCH --error=", job_name, "_%A_%a.err"),
    paste0("#SBATCH --time=", wall_time),
    "#SBATCH --nodes=1",
    "#SBATCH --ntasks=1",
    paste0("#SBATCH --cpus-per-task=", cpus_per_task),
    paste0("#SBATCH --mem=", memory),
    paste0("#SBATCH --array=1-", total_combinations),
    "# Run the R script with the combinations file",
    paste0("Rscript -e \"CONCERTDR::process_combinations_file('", combinations_file, "')\"")
  )
  
  # Write SLURM script to file
  writeLines(slurm_script, output_file)
  
  # Make the script executable if requested
  if (make_executable) {
    Sys.chmod(output_file, mode = "755")
    message(sprintf("Made script executable: %s", output_file))
  }
  
  message(sprintf("SLURM submission script generated: %s", output_file))
  message(sprintf("To submit the job: sbatch %s", output_file))
  
  return(invisible(output_file))
}

#' Interactive selection of signature matching methods
#'
#' @return Vector of selected method names
#' @export
select_methods_interactive <- function() {
  methods <- c("ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang")
  descriptions <- c(
    "KS Score - Kolmogorov-Smirnov statistic based method",
    "XCos Score - Extended cosine similarity method",
    "XSum Score - Extended sum score method",
    "GSEA Weight 0 Score - Gene Set Enrichment Analysis without weighting",
    "GSEA Weight 1 Score - Gene Set Enrichment Analysis with weight 1",
    "GSEA Weight 2 Score - Gene Set Enrichment Analysis with weight 2",
    "Zhang Score - Method from Zhang et al."
  )

  cat("Available signature matching methods:\n")
  for (i in 1:length(methods)) {
    cat(i, ": ", descriptions[i], "\n", sep = "")
  }
  cat("8: Select all methods\n")

  # Get user selection
  selection_input <- readline(prompt = "Enter numbers of methods to use (comma-separated, e.g. 1,3,5): ")

  # Handle empty input
  if (selection_input == "") {
    message("No selection made. Using all methods.")
    return(methods)
  }

  # Parse selection
  selection <- as.numeric(unlist(strsplit(selection_input, ",")))

  # Process selection
  if (length(selection) == 1 && selection == 8) {
    return(methods)
  } else {
    # Remove out-of-range values
    selection <- selection[selection >= 1 & selection <= length(methods)]
    if (length(selection) == 0) {
      message("No valid methods selected. Using KS Score as default.")
      return("ks")
    }
    return(methods[selection])
  }
}

#' Interactive CMap analysis workflow
#'
#' Provides an interactive command-line interface to the CMap analysis workflow.
#'
#' @return Invisibly returns NULL
#' @export
interactive_cmap_workflow <- function() {
  # Ask the user which mode they want to run in
  cat("=== CMap Signature Matching Tool ===\n")
  cat("Please select run mode:\n")
  cat("1: Use combinations.txt from generate_combinations_from_config\n")
  cat("2: Find reference files manually (filtered_*.csv)\n")

  # Get user selection for mode
  mode_selection <- readline(prompt = "Enter selection (1-2, default: 1): ")

  if (mode_selection == "") mode_selection <- "1"
  mode_selection <- as.numeric(mode_selection)

  # Get signature file
  sig_file <- readline(prompt = "Enter signature file path (default: signature.txt): ")
  if (sig_file == "") sig_file <- "signature.txt"

  # Get output directory
  out_dir <- readline(prompt = "Enter output directory (default: results): ")
  if (out_dir == "") out_dir <- "results"

  # Ask about methods
  cat("Do you want to select specific scoring methods? (y/n, default: n): ")
  select_methods <- readline(prompt = "")

  if (tolower(select_methods) == "y") {
    methods <- select_methods_interactive()
  } else {
    methods <- c("ks", "xcos", "xsum", "gsea0", "gsea1", "gsea2", "zhang")
  }

  # Get permutations
  perm_input <- readline(prompt = "Enter number of permutations (default: 100): ")
  permutations <- if (perm_input == "") 100 else as.numeric(perm_input)

  # Run the appropriate mode
  if (mode_selection == 1) {
    # Use combinations.txt mode
    comb_file <- readline(prompt = "Enter combinations file path (default: combinations.txt): ")
    if (comb_file == "") comb_file <- "combinations.txt"

    run_analysis_with_combinations(
      combinations_file = comb_file,
      sig_file = sig_file,
      out_dir = out_dir,
      methods = methods,
      permutations = permutations
    )
  } else {
    # Manual mode
    ref_pattern <- readline(prompt = "Enter reference file pattern (default: filtered_*.csv): ")
    if (ref_pattern == "") ref_pattern <- "filtered_*.csv"

    run_analysis(
      ref_pattern = ref_pattern,
      sig_file = sig_file,
      out_dir = out_dir,
      methods = methods,
      permutations = permutations
    )
  }

  message("Analysis complete!")
  return(invisible(NULL))
}
