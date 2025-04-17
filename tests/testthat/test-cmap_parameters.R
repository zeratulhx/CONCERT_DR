<<<<<<< HEAD
context("CMap Parameter Extraction")

test_that("write_config_file works", {
  # Create a temporary file
  temp_file <- tempfile(fileext = ".conf")
  on.exit(unlink(temp_file), add = TRUE)

  # Create a simple options list
  options_list <- list(
    times = c("6 h", "24 h"),
    doses = c("10 uM", "1 uM"),
    cells = c("MCF7", "A375")
  )

  # Write the config file
  result <- write_config_file(options_list, temp_file)

  # Check that the file exists
  expect_true(file.exists(temp_file))

  # Read the file contents
  lines <- readLines(temp_file)

  # Check that key sections exist
  expect_true(any(grepl("\\[times\\]", lines)))
  expect_true(any(grepl("\\[doses\\]", lines)))
  expect_true(any(grepl("\\[cells\\]", lines)))
  expect_true(any(grepl("\\[selected\\]", lines)))
})

test_that("read_cmap_config works", {
  # Create a temporary file
  temp_file <- tempfile(fileext = ".conf")
  on.exit(unlink(temp_file), add = TRUE)

  # Create a simple config file
  writeLines(c(
    "# CMap test config",
    "[times]",
    "1 = 6 h",
    "2 = 24 h",
    "",
    "[doses]",
    "1 = 10 uM",
    "2 = 1 uM",
    "",
    "[cells]",
    "1 = MCF7",
    "2 = A375",
    "",
    "[selected]",
    "times = 1,2",
    "doses = 1",
    "cells = 2"
  ), temp_file)

  # Read the config
  config <- read_cmap_config(temp_file)

  # Check the structure
  expect_type(config, "list")
  expect_true(all(c("all_options", "selected") %in% names(config)))

  # Check selected values
  expect_equal(config$selected$times, c("6 h", "24 h"))
  expect_equal(config$selected$doses, c("10 uM"))
  expect_equal(config$selected$cells, c("A375"))
=======
test_that("CMap parameter extraction and combination generation work with real data", {
  # Create a temporary directory for our test artifacts
  test_temp_dir <- tempdir()
  config_file <- file.path(test_temp_dir, "real_data_test.conf")

  # 1. Test extract_cmap_parameters
  test_siginfo <- system.file("extdata", "siginfo_beta.txt", package = "DRnew")

  # Skip if data file doesn't exist
  skip_if_not(file.exists(test_siginfo), "Test siginfo file not found")

  # Run the parameter extraction
  params <- extract_cmap_parameters(
    siginfo_file = test_siginfo,
    config_dir = test_temp_dir,
    config_filename = "real_data_test.conf",
    verbose = FALSE
  )

  # Verify expected values are in the parameters
  expect_true("6 h" %in% params$times)
  expect_true("10 uM" %in% params$doses)
  expect_true("MCF7" %in% params$cells)

  # Check that config file was created
  expect_true(file.exists(config_file))

  # 2. Test generate_combinations_from_config with the file we just created

  # Edit the config file to limit selections (for faster testing)
  # Read the config file
  config_lines <- readLines(config_file)

  # Find the [selected] section
  selected_idx <- which(grepl("\\[selected\\]", config_lines))
  if (length(selected_idx) > 0) {
    # Update to select just the first option for each parameter
    # Find where each parameter is defined in the selected section
    times_idx <- grep("^times\\s*=", config_lines[selected_idx:length(config_lines)])
    doses_idx <- grep("^doses\\s*=", config_lines[selected_idx:length(config_lines)])
    cells_idx <- grep("^cells\\s*=", config_lines[selected_idx:length(config_lines)])

    # Update the selections (adding the offset from selected_idx)
    if (length(times_idx) > 0)
      config_lines[selected_idx + times_idx[1] - 1] <- "times = 1"
    if (length(doses_idx) > 0)
      config_lines[selected_idx + doses_idx[1] - 1] <- "doses = 1"
    if (length(cells_idx) > 0)
      config_lines[selected_idx + cells_idx[1] - 1] <- "cells = 1"

    # Write the updated config back
    writeLines(config_lines, config_file)
  }

  # Generate combinations
  combinations <- generate_combinations_from_config(config_file, verbose = FALSE)

  # Check results based on known config contents
  expect_s3_class(combinations, "data.frame")
  expect_true(all(c("itime", "idose", "cell") %in% names(combinations)))
  expect_equal(nrow(combinations), 1)  # Should have 1 combination (1 time x 1 dose x 1 cell)

  # 3. Save combinations to a file for potential further testing
  combinations_file <- file.path(test_temp_dir, "test_combinations.txt")
  write.table(combinations, file = combinations_file, row.names = FALSE)
  expect_true(file.exists(combinations_file))
>>>>>>> bundle_remote/test
})
