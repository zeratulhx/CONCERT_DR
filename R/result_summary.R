#' Annotate CONCERTDR Results with Drug Information
#'
#' @description
#' This function annotates signature matching results from CONCERTDR with drug
#' information including compound names, aliases, mechanism of action (MOA),
#' and clinical phase data. It performs fuzzy string matching to find the best
#' matches when exact matches are not available.
#'
#' @param results_df Data frame containing signature matching results with a 'compound' column
#' @param sig_info_file Path to siginfo_beta.txt file or data frame with signature information
#' @param comp_info_file Path to compound information file or data frame
#' @param drug_info_file Path to drug information file or data frame with clinical phase data
#' @param output_file Optional path to save annotated results (default: NULL, no file saved)
#' @param fuzzy_threshold Minimum similarity score for fuzzy matching (default: 90)
#' @param perfect_match_only Logical; if TRUE, only accept perfect matches (score=100) for fuzzy matching (default: FALSE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return Data frame with annotated results including drug names, aliases, MOA, and clinical phases
#'
#' @examples
#' \dontrun{
#' # Read signature matching results
#' results <- read.csv("sig_match_xsum_results.csv")
#'
#' # Annotate with drug information
#' annotated_results <- annotate_drug_results(
#'   results_df = results,
#'   sig_info_file = "siginfo_beta.txt",
#'   comp_info_file = "compound_info.txt",
#'   drug_info_file = "drug_info.txt",
#'   output_file = "annotated_results.csv"
#' )
#'
#' # View results with drug information
#' head(annotated_results[, c("compound", "Score", "name", "moa", "phase")])
#' }
#'
#' @export
annotate_drug_results <- function(results_df,
                                  sig_info_file,
                                  comp_info_file,
                                  drug_info_file,
                                  output_file = NULL,
                                  fuzzy_threshold = 90,
                                  perfect_match_only = FALSE,
                                  verbose = TRUE) {

  # Check required packages
  required_packages <- c("RecordLinkage")
  missing_packages <- setdiff(required_packages, installed.packages()[,"Package"])

  if (length(missing_packages) > 0) {
    if (verbose) {
      message("Installing required packages for fuzzy matching: ", paste(missing_packages, collapse = ", "))
    }
    install.packages(missing_packages, quiet = !verbose)
  }

  # Load required library for fuzzy matching
  if (!requireNamespace("RecordLinkage", quietly = TRUE)) {
    stop("Package 'RecordLinkage' is required for fuzzy string matching. Please install it with: install.packages('RecordLinkage')")
  }

  if (verbose) message("Starting drug annotation process...")

  # Input validation
  if (!"compound" %in% colnames(results_df)) {
    stop("Results data frame must contain a 'compound' column")
  }

  # Helper function to read files, use data frames, or check environment
  read_or_use_df <- function(file_or_df, file_description, skip_rows = 0, env_var_name = NULL) {
    # If it's a data frame, return it directly
    if (is.data.frame(file_or_df)) {
      if (verbose) message("Using provided data frame for ", file_description)
      return(file_or_df)
    }

    # If it's a character (file path)
    if (is.character(file_or_df)) {
      # Check if there's a corresponding variable in the environment
      if (!is.null(env_var_name) && exists(env_var_name, envir = .GlobalEnv)) {
        env_data <- get(env_var_name, envir = .GlobalEnv)
        if (is.data.frame(env_data)) {
          if (verbose) message("Found ", file_description, " already loaded in environment as '", env_var_name, "'")
          return(env_data)
        }
      }

      # Check if file exists
      if (!file.exists(file_or_df)) {
        stop(file_description, " file not found: ", file_or_df)
      }

      # Read the file
      if (verbose) message("Reading ", file_description, " from file: ", file_or_df)

      if (skip_rows > 0) {
        data <- utils::read.table(file_or_df, sep = "\t", header = TRUE,
                                  stringsAsFactors = FALSE, quote = "",
                                  comment.char = "", fill = TRUE, skip = skip_rows)
      } else {
        data <- utils::read.table(file_or_df, sep = "\t", header = TRUE,
                                  stringsAsFactors = FALSE, quote = "",
                                  comment.char = "", fill = TRUE)
      }

      # Optionally cache in environment for future use
      if (!is.null(env_var_name)) {
        assign(env_var_name, data, envir = .GlobalEnv)
        if (verbose) message("Cached ", file_description, " in environment as '", env_var_name, "'")
      }

      return(data)
    }

    stop(file_description, " must be either a file path (character) or a data.frame")
  }

  # Read input files (check environment first)
  if (verbose) message("Reading signature information...")
  sig_info <- read_or_use_df(sig_info_file, "Signature info", env_var_name = "sig_info_cached")

  if (verbose) message("Reading compound information...")
  df_comp <- read_or_use_df(comp_info_file, "Compound info", env_var_name = "comp_info_cached")

  if (verbose) message("Reading drug information...")
  df_drug <- read_or_use_df(drug_info_file, "Drug info", skip_rows = 9, env_var_name = "drug_info_cached")

  # Make a copy of the results to avoid modifying the original
  df_output <- results_df

  # Step 1: Extract compound info data (names, MOA) for annotation
  if (verbose) message("Mapping compound identifiers to perturbation IDs...")

  # Create mapping from sig_id to pert_id
  if (!"sig_id" %in% colnames(sig_info) || !"pert_id" %in% colnames(sig_info)) {
    stop("sig_info must contain 'sig_id' and 'pert_id' columns")
  }

  # Extract sig_id from compound column
  # Extract core BRD identifier (e.g., BRD-K22134346 from BRD-K22134346-001-11-6)
  df_output$sig_id <- sapply(strsplit(df_output$compound, ":"), function(x) {
    if (length(x) >= 2) {
      # Extract the BRD code (second part)
      brd_code <- x[2]
      # Extract only the core BRD identifier (BRD-[letter][numbers])
      # This removes suffixes like -001-11-6
      core_brd <- sub("^(BRD-[A-Z][0-9]+).*", "\\1", brd_code)
      return(core_brd)
    } else {
      return(x[1])  # Fallback to the whole string
    }
  })

  # Map sig_id to pert_id
  sig_to_pert <- stats::setNames(sig_info$cmap_name, sig_info$pert_id)
  df_output$pert_name <- sig_to_pert[df_output$sig_id]

  # Map compound information
  if (verbose) message("Mapping compound names and MOA...")

  # Map aliases
  if ("compound_aliases" %in% colnames(df_comp) && "pert_id" %in% colnames(df_comp)) {
    pert_to_aliases <- stats::setNames(df_comp$compound_aliases, df_comp$pert_id)
    df_output$aliases <- pert_to_aliases[df_output$sig_id]
  }

  # Step 2: Annotate compounds with clinical phase data
  if (verbose) message("Annotating clinical phase information...")

  # Convert drug names to lowercase for matching
  if ("pert_iname" %in% colnames(df_drug)) {
    df_drug$pert_iname <- tolower(df_drug$pert_iname)
  } else {
    warning("'pert_iname' column not found in drug info file")
    return(df_output)
  }

  # Initialize new columns
  df_output$phase <- NA
  df_output$is_poss_pair <- 0
  df_output$pair <- NA
  df_output$pair_score <- NA

  # Process each row
  for (i in 1:nrow(df_output)) {
    if (verbose && i %% 100 == 0) {
      message("Processing row ", i, " of ", nrow(df_output))
    }

    # Get names for current row
    name1 <- if (!is.na(df_output$pert_name[i])) tolower(as.character(df_output$pert_name[i])) else ""
    name2 <- if (!is.na(df_output$aliases[i])) tolower(as.character(df_output$aliases[i])) else ""

    # Try exact matches first
    exact_match_found <- FALSE

    if (name1 != "" && name1 %in% df_drug$pert_iname) {
      matched_drug <- df_drug[df_drug$pert_iname == name1, ]
      if (nrow(matched_drug) > 0 && "clinical_phase" %in% colnames(matched_drug)) {
        df_output$phase[i] <- matched_drug$clinical_phase[1]
        exact_match_found <- TRUE
      }
    } else if (name2 != "" && name2 %in% df_drug$pert_iname) {
      matched_drug <- df_drug[df_drug$pert_iname == name2, ]
      if (nrow(matched_drug) > 0 && "clinical_phase" %in% colnames(matched_drug)) {
        df_output$phase[i] <- matched_drug$clinical_phase[1]
        exact_match_found <- TRUE
      }
    }

    # If no exact match, try fuzzy matching
    if (!exact_match_found && (name1 != "" || name2 != "")) {
      score1 <- 0
      score2 <- 0
      pair1 <- ""
      pair2 <- ""

      # Fuzzy match for name1
      if (name1 != "") {
        similarities1 <- RecordLinkage::levenshteinSim(name1, df_drug$pert_iname)
        max_idx1 <- which.max(similarities1)
        if (length(max_idx1) > 0) {
          score1 <- similarities1[max_idx1] * 100  # Convert to percentage
          pair1 <- df_drug$pert_iname[max_idx1]
        }
      }

      # Fuzzy match for name2
      if (name2 != "") {
        similarities2 <- RecordLinkage::levenshteinSim(name2, df_drug$pert_iname)
        max_idx2 <- which.max(similarities2)
        if (length(max_idx2) > 0) {
          score2 <- similarities2[max_idx2] * 100  # Convert to percentage
          pair2 <- df_drug$pert_iname[max_idx2]
        }
      }

      # Choose the best match
      if (score1 > fuzzy_threshold || score2 > fuzzy_threshold) {
        if (score1 > score2) {
          matched_name <- pair1
          score <- score1
        } else {
          matched_name <- pair2
          score <- score2
        }

        # Apply matching criteria
        apply_match <- if (perfect_match_only) {
          score == 100
        } else {
          score > fuzzy_threshold
        }

        if (apply_match) {
          df_output$is_poss_pair[i] <- 1
          matched_drug <- df_drug[df_drug$pert_iname == matched_name, ]
          if (nrow(matched_drug) > 0 && "clinical_phase" %in% colnames(matched_drug)) {
            df_output$phase[i] <- matched_drug$clinical_phase[1]
          }
          df_output$pair[i] <- matched_name
          df_output$pair_score[i] <- score
        }
      }
    }
  }

  # Save to file if requested
  if (!is.null(output_file)) {
    if (verbose) message("Saving annotated results to: ", output_file)
    utils::write.csv(df_output, file = output_file, row.names = FALSE)
  }

  if (verbose) {
    # Print summary statistics
    total_rows <- nrow(df_output)
    with_phase <- sum(!is.na(df_output$phase))
    fuzzy_matches <- sum(df_output$is_poss_pair == 1, na.rm = TRUE)

    message("\nAnnotation Summary:")
    message("Total compounds: ", total_rows)
    message("With clinical phase: ", with_phase, " (", round(with_phase/total_rows*100, 1), "%)")
    message("Fuzzy matches used: ", fuzzy_matches, " (", round(fuzzy_matches/total_rows*100, 1), "%)")
  }

  return(df_output)
}

#' Extract Compound ID from Signature Matching Results
#'
#' @description
#' Helper function to extract compound identifiers from complex compound strings
#' in signature matching results. This function can be customized based on the
#' specific format of compound identifiers in your data.
#'
#' @param compound_strings Vector of compound identifier strings
#' @param method Method for extraction: "split_colon" (default), "split_underscore", or "regex"
#' @param regex_pattern Regular expression pattern for extraction (used when method="regex")
#' @param part_index Which part to extract when splitting (default: 2)
#'
#' @return Vector of extracted compound identifiers
#'
#' @examples
#' \dontrun{
#' # Example compound strings
#' compounds <- c("CVD001_HEPG2_6H:BRD-K03652504-001-01-9:10.0497",
#'                "CVD001_HEPG2_6H:BRD-A37828317-001-03-0:10")
#'
#' # Extract using colon splitting (default)
#' ids <- extract_compound_id(compounds)
#'
#' # Extract using custom regex
#' ids <- extract_compound_id(compounds, method = "regex",
#'                           regex_pattern = "BRD-[A-Z0-9-]+")
#' }
#'
#' @export
extract_compound_id <- function(compound_strings,
                                method = "split_colon",
                                regex_pattern = NULL,
                                part_index = 2) {

  if (method == "split_colon") {
    return(sapply(strsplit(compound_strings, ":"), function(x) {
      if (length(x) >= part_index) {
        return(x[part_index])
      } else {
        return(x[1])
      }
    }))
  } else if (method == "split_underscore") {
    return(sapply(strsplit(compound_strings, "_"), function(x) {
      if (length(x) >= part_index) {
        return(x[part_index])
      } else {
        return(x[1])
      }
    }))
  } else if (method == "regex") {
    if (is.null(regex_pattern)) {
      stop("regex_pattern must be provided when method='regex'")
    }
    matches <- regmatches(compound_strings, regexpr(regex_pattern, compound_strings))
    return(ifelse(matches == "", compound_strings, matches))
  } else {
    stop("Invalid method. Must be 'split_colon', 'split_underscore', or 'regex'")
  }
}

#' Fuzzy String Matching for Drug Names
#'
#' @description
#' Performs fuzzy string matching to find the best matches between query drug names
#' and a reference list of drug names. Uses Levenshtein distance for similarity calculation.
#'
#' @param query_names Vector of query drug names
#' @param reference_names Vector of reference drug names to match against
#' @param method Similarity method: "levenshtein" (default), "jaro", or "jarowinkler"
#' @param threshold Minimum similarity threshold (0-100)
#' @param top_n Number of top matches to return for each query (default: 1)
#'
#' @return Data frame with query names, matched names, and similarity scores
#'
#' @examples
#' \dontrun{
#' query <- c("aspirin", "ibuprofen", "acetaminophen")
#' reference <- c("aspirin", "ibuprofen", "acetaminophen", "naproxen", "diclofenac")
#'
#' matches <- fuzzy_drug_match(query, reference, threshold = 80)
#' print(matches)
#' }
#'
#' @export
fuzzy_drug_match <- function(query_names,
                             reference_names,
                             method = "levenshtein",
                             threshold = 80,
                             top_n = 1) {

  if (!requireNamespace("RecordLinkage", quietly = TRUE)) {
    stop("Package 'RecordLinkage' is required for fuzzy matching. Please install it with: install.packages('RecordLinkage')")
  }

  results <- data.frame()

  for (i in seq_along(query_names)) {
    query <- query_names[i]

    if (is.na(query) || query == "") {
      next
    }

    # Calculate similarities
    if (method == "levenshtein") {
      similarities <- RecordLinkage::levenshteinSim(query, reference_names) * 100
    } else if (method == "jaro") {
      similarities <- RecordLinkage::jarowinkler(query, reference_names) * 100
    } else if (method == "jarowinkler") {
      similarities <- RecordLinkage::jarowinkler(query, reference_names, W = 0.1) * 100
    } else {
      stop("Invalid method. Must be 'levenshtein', 'jaro', or 'jarowinkler'")
    }

    # Find top matches above threshold
    above_threshold <- which(similarities >= threshold)
    if (length(above_threshold) > 0) {
      sorted_indices <- above_threshold[order(similarities[above_threshold], decreasing = TRUE)]
      top_indices <- head(sorted_indices, top_n)

      for (idx in top_indices) {
        results <- rbind(results, data.frame(
          query_name = query,
          matched_name = reference_names[idx],
          similarity_score = similarities[idx],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(results)
}
