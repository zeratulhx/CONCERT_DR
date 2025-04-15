context("CMap Combinations Generation")

test_that("generate_combinations_from_config works", {
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

  # Generate combinations
  suppressMessages(
    combinations <- generate_combinations_from_config(temp_file, verbose = FALSE)
  )

  # Check the structure
  expect_s3_class(combinations, "data.frame")
  expect_equal(nrow(combinations), 2) # 2 times × 1 dose × 1 cell
  expect_true(all(c("itime", "idose", "cell") %in% names(combinations)))

  # Check the values
  expect_equal(as.character(combinations$itime), c("6 h", "24 h"))
  expect_equal(as.character(combinations$idose), c("10 uM", "10 uM"))
  expect_equal(as.character(combinations$cell), c("A375", "A375"))
})
