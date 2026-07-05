library(testthat)
library(PhysioEDA)

test_that("edaFeatures returns data.frame with correct columns", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x)
  feat <- edaFeatures(x)

  expect_s3_class(feat, "data.frame")
  expected_cols <- c("channel", "scr_count", "scr_rate_per_min",
                     "mean_amplitude", "mean_scl", "scl_sd",
                     "auc_phasic", "ns_scr_freq")
  expect_equal(colnames(feat), expected_cols)
})

test_that("edaFeatures has one row per channel", {
  set.seed(42)
  x <- make_eda(n_channels = 3)
  x <- edaDecompose(x)
  feat <- edaFeatures(x)

  expect_equal(nrow(feat), 3)
})

test_that("edaFeatures mean_scl is positive", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x)
  feat <- edaFeatures(x)

  expect_true(all(feat$mean_scl > 0))
})

test_that("edaFeatures accepts pre-computed peaks", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x)
  peaks <- edaPeaks(x)
  feat <- edaFeatures(x, peaks = peaks)

  expect_equal(feat$scr_count, nrow(peaks))
})

test_that("edaFeatures supports time window", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x)
  feat <- edaFeatures(x, window = c(100, 300))

  expect_s3_class(feat, "data.frame")
})
