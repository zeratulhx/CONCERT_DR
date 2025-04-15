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
})
