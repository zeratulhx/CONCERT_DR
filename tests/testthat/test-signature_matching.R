test_that("signature matching works with test data", {
  # Create test directory
  test_temp_dir <- tempdir()

  # Get paths to test files - make sure these exist
  test_signature <- system.file("extdata", "test_data", "signature.txt", package = "DRnew")
  test_reference <- system.file("extdata", "test_data", "filtered_6_h_10_uM_HEPG2.csv", package = "DRnew")

  # Skip if files don't exist or required package not available
  skip_if_not(file.exists(test_signature), "Test signature file not found")
  skip_if_not(file.exists(test_reference), "Test reference file not found")
  skip_if_not(requireNamespace("RCSM", quietly = TRUE), "RCSM package not available")

  # Create output directory for results
  results_dir <- file.path(test_temp_dir, "results")
  dir.create(results_dir, recursive = TRUE)

  # Run the function with simplified methods to reduce dependencies
  result <- process_signature(
    signature_file = test_signature,
    reference_file = test_reference,
    output_prefix = file.path(results_dir, "test_output"),
    permutations = 10,  # Use small number for faster tests
    methods = c("ks")   # Just use one method for faster tests
  )

  # Check that results have expected structure
  expect_type(result, "list")
  expect_true("ks" %in% names(result))

  # Check output file exists
  expect_true(file.exists(file.path(results_dir, "test_output_ks_results.csv")))
})
