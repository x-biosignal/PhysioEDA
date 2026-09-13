library(testthat)
library(PhysioEDA)

test_that("edaPeaks returns data.frame with correct columns", {
  set.seed(42)
  x <- make_eda()
  peaks <- edaPeaks(x)

  expect_s3_class(peaks, "data.frame")
  expected_cols <- c("channel", "onset_sample", "onset_sec",
                     "peak_sample", "peak_sec", "amplitude",
                     "rise_time", "recovery_time")
  expect_equal(colnames(peaks), expected_cols)
})

test_that("edaPeaks detects SCRs in simulated signal", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 1, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.001)
  x <- edaDecompose(x)
  peaks <- edaPeaks(x)

  expect_true(nrow(peaks) > 0)
})

test_that("edaPeaks respects amplitude_min threshold", {
  set.seed(42)
  x <- make_eda()
  peaks <- edaPeaks(x, amplitude_min = 0.01)

  if (nrow(peaks) > 0) {
    expect_true(all(peaks$amplitude >= 0.01))
  }
})

test_that("edaPeaks threshold method works", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 1, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.001)
  x <- edaDecompose(x)
  peaks <- edaPeaks(x, method = "threshold")

  expect_s3_class(peaks, "data.frame")
})

test_that("edaPeaks auto-selects phasic assay", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 1, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.001)
  x <- edaDecompose(x)

  # Call without specifying assay_name; should auto-select "phasic"
  peaks <- edaPeaks(x)

  expect_s3_class(peaks, "data.frame")
  expected_cols <- c("channel", "onset_sample", "onset_sec",
                     "peak_sample", "peak_sec", "amplitude",
                     "rise_time", "recovery_time")
  expect_equal(colnames(peaks), expected_cols)
})

test_that("edaPeaks handles multi-channel data", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 2, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.001)
  x <- edaDecompose(x)
  peaks <- edaPeaks(x)

  expect_s3_class(peaks, "data.frame")
  if (nrow(peaks) > 0) {
    expect_true(all(peaks$channel %in% c("EDA1", "EDA2")))
  }
})

test_that("edaPeaks returns empty data.frame for flat signal", {
  x <- PhysioExperiment(
    assays = list(raw = matrix(rep(5.0, 1000), ncol = 1)),
    colData = S4Vectors::DataFrame(label = "EDA1", type = "EDA"),
    samplingRate = 10
  )
  peaks <- edaPeaks(x)

  expect_s3_class(peaks, "data.frame")
  expect_equal(nrow(peaks), 0)
  expected_cols <- c("channel", "onset_sample", "onset_sec",
                     "peak_sample", "peak_sec", "amplitude",
                     "rise_time", "recovery_time")
  expect_equal(colnames(peaks), expected_cols)
})
