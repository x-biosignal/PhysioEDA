library(testthat)
library(PhysioEDA)

test_that("edaArtifact detects threshold violations", {
  set.seed(42)
  x <- make_eda()
  # Inject out-of-range values
  raw <- SummarizedExperiment::assay(x, "raw")
  raw[10, 1] <- -0.5    # below lower bound
  raw[20, 1] <- 100.0   # above upper bound
  SummarizedExperiment::assay(x, "raw") <- raw

  x <- edaArtifact(x, methods = "threshold", correct = "none")

  mask <- S4Vectors::metadata(x)$eda_artifacts$mask
  expect_true(mask[10, 1])
  expect_true(mask[20, 1])
})

test_that("edaArtifact detects gradient artifacts", {
  set.seed(42)
  x <- make_eda()
  sr <- samplingRate(x)  # 10 Hz
  raw <- SummarizedExperiment::assay(x, "raw")
  # Inject a spike: jump of 5 uS in one sample (50 uS/sec >> 10 uS/sec limit)
  raw[100, 1] <- raw[99, 1] + 5.0
  SummarizedExperiment::assay(x, "raw") <- raw

  x <- edaArtifact(x, methods = "gradient", correct = "none")

  mask <- S4Vectors::metadata(x)$eda_artifacts$mask
  # Both samples around the spike should be flagged
  expect_true(mask[100, 1])
})

test_that("edaArtifact detects flatline segments", {
  set.seed(42)
  x <- make_eda()
  sr <- samplingRate(x)  # 10 Hz
  raw <- SummarizedExperiment::assay(x, "raw")
  # Create a flatline of 6 seconds (60 samples at 10 Hz, > 5 sec threshold)
  raw[200:260, 1] <- 5.0
  SummarizedExperiment::assay(x, "raw") <- raw

  x <- edaArtifact(x, methods = "flatline", flatline_sec = 5,
                   correct = "none")

  mask <- S4Vectors::metadata(x)$eda_artifacts$mask
  # Flatline region should be flagged
  expect_true(all(mask[200:260, 1]))
})

test_that("edaArtifact interpolate correction produces no artifacts", {
  set.seed(42)
  x <- make_eda()
  raw <- SummarizedExperiment::assay(x, "raw")
  # Inject out-of-range values
  raw[50, 1] <- -1.0
  raw[51, 1] <- -1.0
  SummarizedExperiment::assay(x, "raw") <- raw

  x <- edaArtifact(x, methods = "threshold", correct = "interpolate")

  cleaned <- SummarizedExperiment::assay(x, "cleaned")
  # Corrected values should be interpolated (within range)
  expect_true(cleaned[50, 1] > 0)
  expect_true(cleaned[51, 1] > 0)
  # Original raw should still have the bad values
  raw_after <- SummarizedExperiment::assay(x, "raw")
  expect_equal(raw_after[50, 1], -1.0)
})

test_that("edaArtifact NA correction replaces artifacts with NA", {
  set.seed(42)
  x <- make_eda()
  raw <- SummarizedExperiment::assay(x, "raw")
  raw[50, 1] <- -1.0
  SummarizedExperiment::assay(x, "raw") <- raw

  x <- edaArtifact(x, methods = "threshold", correct = "na")

  cleaned <- SummarizedExperiment::assay(x, "cleaned")
  expect_true(is.na(cleaned[50, 1]))
})

test_that("edaArtifact correct='none' does not add output assay", {
  set.seed(42)
  x <- make_eda()

  x <- edaArtifact(x, correct = "none")

  nms <- SummarizedExperiment::assayNames(x)
  expect_false("cleaned" %in% nms)
  # But metadata should still be present
  expect_true(!is.null(S4Vectors::metadata(x)$eda_artifacts))
})

test_that("edaArtifact summary contains counts per method", {
  set.seed(42)
  x <- make_eda()
  raw <- SummarizedExperiment::assay(x, "raw")
  raw[10, 1] <- -1.0  # threshold violation
  SummarizedExperiment::assay(x, "raw") <- raw

  x <- edaArtifact(x, methods = c("threshold", "gradient"),
                   correct = "none")

  summary_df <- S4Vectors::metadata(x)$eda_artifacts$summary
  expect_true(is.data.frame(summary_df))
  expect_true(all(c("channel", "method", "n_artifacts", "pct") %in%
                    colnames(summary_df)))
  # Should have rows for each method x channel combination
  expect_equal(nrow(summary_df), 2)  # 1 channel x 2 methods
  # Threshold should have detected at least the injected artifact
  thresh_row <- summary_df[summary_df$method == "threshold", ]
  expect_gte(thresh_row$n_artifacts, 1)
})

test_that("edaArtifact works with multiple channels", {
  set.seed(42)
  x <- make_eda(n_channels = 3)
  raw <- SummarizedExperiment::assay(x, "raw")
  raw[10, 1] <- -1.0
  raw[20, 2] <- 100.0
  SummarizedExperiment::assay(x, "raw") <- raw

  x <- edaArtifact(x, methods = "threshold", correct = "interpolate")

  mask <- S4Vectors::metadata(x)$eda_artifacts$mask
  expect_equal(ncol(mask), 3)
  expect_true(mask[10, 1])
  expect_true(mask[20, 2])
  # Channel 3 should have no threshold artifacts from injected values
  expect_false(mask[10, 3])
  expect_false(mask[20, 3])

  # Cleaned assay should exist
  expect_true("cleaned" %in% SummarizedExperiment::assayNames(x))
})

test_that("edaArtifact mask dimensions match input", {
  set.seed(42)
  x <- make_eda(n_channels = 2)

  x <- edaArtifact(x, correct = "none")

  mask <- S4Vectors::metadata(x)$eda_artifacts$mask
  raw_dim <- dim(SummarizedExperiment::assay(x, "raw"))
  expect_equal(dim(mask), raw_dim)
})

test_that("edaArtifact custom output_assay name works", {
  set.seed(42)
  x <- make_eda()

  x <- edaArtifact(x, correct = "interpolate", output_assay = "my_clean")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("my_clean" %in% nms)
  expect_false("cleaned" %in% nms)
})

test_that("edaArtifact gradient_max parameter is respected", {
  set.seed(42)
  x <- make_eda()
  raw <- SummarizedExperiment::assay(x, "raw")
  # Create a moderate gradient that exceeds a strict threshold
  raw[100, 1] <- raw[99, 1] + 0.5
  SummarizedExperiment::assay(x, "raw") <- raw

  # With a very strict gradient_max, this should be detected
  x_strict <- edaArtifact(x, methods = "gradient", gradient_max = 0.01,
                           correct = "none")
  mask_strict <- S4Vectors::metadata(x_strict)$eda_artifacts$mask
  expect_true(mask_strict[100, 1])

  # With a lenient gradient_max, it should not be detected
  x_lenient <- edaArtifact(x, methods = "gradient", gradient_max = 10.0,
                            correct = "none")
  mask_lenient <- S4Vectors::metadata(x_lenient)$eda_artifacts$mask
  expect_false(mask_lenient[100, 1])
})
