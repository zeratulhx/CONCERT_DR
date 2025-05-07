#' @title Demonstrate CONCERTDR Workflow
#' @description
#' This function demonstrates the full workflow of the CONCERTDR package
#' with example data. Users can run this to see how the package works
#' in practice. It extracts a small subset of data from a GCTX file
#' to create a realistic reference for signature matching.
#'
#' @param demo_dir Directory to store demonstration files
#' @param use_minimal Logical; use minimal settings for faster demo
#' @param gctx_file Path to GCTX file (if available)
#' @param geneinfo_file Path to geneinfo file (if available)
#' @param siginfo_file Path to siginfo file (if available)
#' @export
demonstrate_workflow <- function(demo_dir = "CONCERTDR_demo",
                                 use_minimal = TRUE,
                                 gctx_file="inst/extdata/level5_beta_trt_cp_n720216x12328.gctx",
                                 geneinfo_file="inst/extdata/geneinfo_beta.txt",
                                 siginfo_file="inst/extdata/siginfo_beta.txt") {
  # Create demo directory if it doesn't exist
  if (!dir.exists(demo_dir)) {
    dir.create(demo_dir, recursive = TRUE)
    message("Created demonstration directory: ", demo_dir)
  }

  # Create subdirectories
  conf_dir <- file.path(demo_dir, "conf")
  output_dir <- file.path(demo_dir, "output")
  results_dir <- file.path(demo_dir, "results")

  for (dir in c(conf_dir, output_dir, results_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
  }


  message("\n--- STEP 1: Creating Configuration Template ---")
  # Create a template config file
  config_file <- create_cmap_config_template(
    dest_dir = conf_dir,
    template_name = "cmap_options_template.conf",
    overwrite = TRUE
  )

  message("\n--- STEP 3: Creating Example Signature ---")
  # Load example signature file with some genes
  example_sig_file <- file.path(demo_dir, "signature.txt")

  message("\n--- STEP 4: Creating Reference Data ---")
  ref_file <- file.path(demo_dir, "reference_data.csv")

  # Check if we have real CMap data files
  have_real_data <- !is.null(gctx_file) && !is.null(geneinfo_file) && !is.null(siginfo_file) &&
    file.exists(gctx_file) && file.exists(geneinfo_file) && file.exists(siginfo_file)

  if (have_real_data) {
    message("Using real CMap data files to extract reference data")

    # Use the package function to process the combinations
    extracted_files <- try({
      # Process the combinations using the package function
      ref_df <- extract_cmap_data_from_config(
        config_file = config_file,
        geneinfo_file = geneinfo_file,
        siginfo_file = siginfo_file,
        gctx_file = gctx_file,
      )
    })
  }
  message("\n--- STEP 5: Signature Matching ---")
  # Run signature matching with the reference data
  methods_to_use <- c("ks", "xsum")  # Use just two methods for demo

  result <- try({
    process_signature_with_df(
      signature_file = example_sig_file,
      reference_df = ref_df,
      output_dir = results_dir,
      permutations = 10,  # Small number for demo
      methods = methods_to_use
    )
  })

  if (inherits(result, "try-error")) {
    message("Note: Signature matching requires the RCSM package. To install:")
    message("    devtools::install_github('Jasonlinchina/RCSM')")
  } else {
    message("Successfully ran signature matching with methods: ", paste(methods_to_use, collapse = ", "))
    message("Results saved to: ", file.path(results_dir))
  }

  message("\n--- WORKFLOW DEMONSTRATION COMPLETE ---")
  message("Demonstration files are in: ", demo_dir)

  invisible(list(
    config_file = config_file,
    signature_file = example_sig_file,
    reference_file = ref_file,
    results_dir = results_dir
  ))
}
