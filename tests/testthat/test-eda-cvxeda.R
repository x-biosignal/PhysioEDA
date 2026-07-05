library(testthat)
library(PhysioEDA)

test_that(".cvxeda_decompose returns correct structure", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  res <- PhysioEDA:::.cvxeda_decompose(signal, sr)

  expect_type(res, "list")
  expect_named(res, c("tonic", "phasic", "driver", "alpha", "gamma",
                       "iterations", "converged"))
})

test_that(".cvxeda_decompose produces correct dimensions", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)
  n <- length(signal)

  res <- PhysioEDA:::.cvxeda_decompose(signal, sr)

  expect_length(res$tonic, n)
  expect_length(res$phasic, n)
  expect_length(res$driver, n)
})

test_that("tonic + phasic approximates original signal (high correlation)", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  res <- PhysioEDA:::.cvxeda_decompose(signal, sr)
  reconstructed <- res$tonic + res$phasic

  # Reconstruction should be very close (tonic is set as signal - phasic)
  expect_equal(reconstructed, signal, tolerance = 1e-10)

  # Correlation should be essentially perfect
  r <- cor(reconstructed, signal)
  expect_gt(r, 0.999)
})

test_that("driver is non-negative", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  res <- PhysioEDA:::.cvxeda_decompose(signal, sr)

  expect_true(all(res$driver >= 0))
})

test_that("higher alpha produces sparser driver", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  res_low <- PhysioEDA:::.cvxeda_decompose(signal, sr, alpha = 0.001)
  res_high <- PhysioEDA:::.cvxeda_decompose(signal, sr, alpha = 0.1)

  nnz_low <- sum(res_low$driver > 1e-10)
  nnz_high <- sum(res_high$driver > 1e-10)

  # Higher alpha should produce fewer non-zero elements (or equal)
  expect_lte(nnz_high, nnz_low)
})

test_that("algorithm converges", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  res <- PhysioEDA:::.cvxeda_decompose(signal, sr, max_iter = 50)

  # Either converged flag is TRUE or finished before max_iter

  expect_true(res$converged || res$iterations < 50)
})

test_that("tonic is smooth (lower variance of differences than raw)", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  # Use low alpha to ensure non-trivial phasic extraction
  res <- PhysioEDA:::.cvxeda_decompose(signal, sr, alpha = 0.001)

  # Tonic should have smaller first-differences (smoother)
  raw_diff_var <- var(diff(signal))
  tonic_diff_var <- var(diff(res$tonic))

  expect_lt(tonic_diff_var, raw_diff_var)
})

test_that(".cvxeda_decompose returns parameters used", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  res <- PhysioEDA:::.cvxeda_decompose(signal, sr, alpha = 0.05, gamma = 0.2)

  expect_equal(res$alpha, 0.05)
  expect_equal(res$gamma, 0.2)
})

test_that(".running_mean produces smoothed output of correct length", {
  x <- sin(seq(0, 4 * pi, length.out = 200)) + rnorm(200, sd = 0.5)

  smoothed <- PhysioEDA:::.running_mean(x, 11)

  expect_length(smoothed, length(x))
  # Smoothed should have less variance
  expect_lt(var(smoothed), var(x))
})
