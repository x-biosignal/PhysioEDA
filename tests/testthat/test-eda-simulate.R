library(testthat)
library(PhysioEDA)

test_that("edaSimulate returns PhysioExperiment with correct dimensions", {
  set.seed(42)
  x <- edaSimulate(n_time = 1000, n_channels = 1, sr = 10)

  expect_s4_class(x, "PhysioExperiment")
  expect_equal(nrow(SummarizedExperiment::assay(x)), 1000)
  expect_equal(ncol(SummarizedExperiment::assay(x)), 1)
  expect_equal(samplingRate(x), 10)
})

test_that("edaSimulate creates multi-channel data", {
  set.seed(42)
  x <- edaSimulate(n_time = 500, n_channels = 3, sr = 10)

  expect_equal(ncol(SummarizedExperiment::assay(x)), 3)
})

test_that("edaSimulate embeds SCRs in signal", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 1, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.001)

  data <- SummarizedExperiment::assay(x)[, 1]
  # Signal should have variation above noise level (SCRs present)
  expect_true(max(data) - min(data) > 0.5)
})

test_that("edaSimulate signal is positive (physiological constraint)", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 1, sr = 10,
                   scl_level = 5.0, noise_sd = 0.01)

  data <- SummarizedExperiment::assay(x)[, 1]
  expect_true(all(data > 0))
})

test_that("edaSimulate respects seed for reproducibility", {
  x1 <- edaSimulate(n_time = 100, sr = 10, seed = 123)
  x2 <- edaSimulate(n_time = 100, sr = 10, seed = 123)

  expect_equal(SummarizedExperiment::assay(x1), SummarizedExperiment::assay(x2))
})

test_that("edaSimulate stores ground-truth tonic/phasic/onsets (WS5-25)", {
  sr <- 10
  x <- edaSimulate(n_time = 600, n_channels = 1, sr = sr, scr_count = 5,
                   scr_amplitude = 0.8, noise_sd = 0.005, seed = 7)
  truth <- S4Vectors::metadata(x)$eda_truth
  expect_type(truth, "list")
  expect_named(truth, c("tonic", "phasic", "onsets"))

  raw <- SummarizedExperiment::assay(x, "raw")
  expect_equal(dim(truth$tonic), dim(raw))
  expect_equal(dim(truth$phasic), dim(raw))
  # raw = tonic + phasic + noise (+ positivity clip); the noise sd is small
  resid <- raw - truth$tonic - truth$phasic
  expect_lt(stats::sd(resid), 0.02)

  on <- truth$onsets
  expect_s3_class(on, "data.frame")
  expect_equal(nrow(on), 5L)
  expect_true(all(c("channel", "sample", "time_sec", "amplitude") %in% names(on)))
  expect_equal(on$time_sec, (on$sample - 1) / sr)
  expect_true(all(on$sample >= 1 & on$sample <= nrow(raw)))
})

test_that("edaSimulate ground-truth onsets scale with channels and scr_count", {
  x <- edaSimulate(n_time = 500, n_channels = 3, sr = 10, scr_count = 4, seed = 1)
  on <- S4Vectors::metadata(x)$eda_truth$onsets
  expect_equal(nrow(on), 3L * 4L)
  expect_setequal(unique(on$channel), 1:3)
})
