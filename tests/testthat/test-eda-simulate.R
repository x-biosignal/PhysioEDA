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
