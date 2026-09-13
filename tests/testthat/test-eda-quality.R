library(testthat)
library(PhysioEDA)

test_that("edaQuality returns 'good' for clean simulated data", {
  set.seed(42)
  x <- make_eda()
  qual <- edaQuality(x)

  expect_equal(qual$quality_label, "good")
})

test_that("edaQuality returns correct columns", {
  set.seed(42)
  x <- make_eda()
  qual <- edaQuality(x)

  expected_cols <- c("channel", "mean_sc", "sd_sc", "min_sc", "max_sc",
                     "pct_negative", "pct_flatline", "pct_artifact",
                     "snr_db", "quality_score", "quality_label")
  expect_equal(colnames(qual), expected_cols)
})

test_that("edaQuality returns one row per channel", {
  set.seed(42)
  x <- make_eda()
  qual <- edaQuality(x)

  expect_equal(nrow(qual), 1)
})

test_that("edaQuality handles multiple channels", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  qual <- edaQuality(x)

  expect_equal(nrow(qual), 2)
  expect_s3_class(qual, "data.frame")
})

test_that("edaQuality quality_score is between 0 and 100", {
  set.seed(42)
  x <- make_eda()
  qual <- edaQuality(x)

  expect_true(all(qual$quality_score >= 0))
  expect_true(all(qual$quality_score <= 100))
})

test_that("edaQuality numeric columns are numeric", {
  set.seed(42)
  x <- make_eda()
  qual <- edaQuality(x)

  numeric_cols <- c("mean_sc", "sd_sc", "min_sc", "max_sc",
                    "pct_negative", "pct_flatline", "pct_artifact",
                    "snr_db", "quality_score")
  for (col in numeric_cols) {
    expect_true(is.numeric(qual[[col]]),
                info = paste("Column", col, "should be numeric"))
  }
})

test_that("edaQuality detects poor quality signal", {
  set.seed(42)
  x <- make_eda()
  # Replace with mostly negative/zero values to trigger penalties
  bad_signal <- matrix(rep(-1, 6000), ncol = 1)
  SummarizedExperiment::assay(x, "raw") <- bad_signal
  qual <- edaQuality(x)

  expect_true(qual$quality_score < 70)
  expect_true(qual$pct_negative == 100)
})
