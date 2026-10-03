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
                     "auc_phasic", "ns_scr_count", "ns_scr_freq",
                     "cda_nscr", "cda_ampsum", "cda_iscr")
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

test_that("edaFeatures exposes CDA/DDA driver metrics when a driver assay exists", {
  fx <- make_eda_scr_clean(seed = 1)
  dec <- edaDecompose(fx$pe, method = "cda")
  feat <- edaFeatures(dec)

  expect_true(all(c("cda_nscr", "cda_ampsum", "cda_iscr") %in% names(feat)))
  expect_false(is.na(feat$cda_nscr[1]))
  expect_gte(feat$cda_nscr[1], 0)
  expect_gte(feat$cda_iscr[1], 0)

  # the exposed metrics equal a direct .cda_metrics call on the driver assay
  drv <- SummarizedExperiment::assay(dec, "driver")[, 1]
  dm <- PhysioEDA:::.cda_metrics(drv, samplingRate(dec))
  expect_equal(feat$cda_nscr[1], dm$nSCR)
  expect_equal(feat$cda_iscr[1], dm$ISCR)
})

test_that("edaFeatures returns NA CDA metrics when no driver assay exists", {
  x <- make_eda()
  x <- edaDecompose(x, method = "median")   # no driver assay
  feat <- edaFeatures(x)
  expect_true(all(c("cda_nscr", "cda_ampsum", "cda_iscr") %in% names(feat)))
  expect_true(is.na(feat$cda_nscr[1]))
  expect_true(is.na(feat$cda_iscr[1]))
})
