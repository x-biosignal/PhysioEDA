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


# ---- Exact cvxEDA quadratic program (WS5-25) --------------------------------

# Ground-truth SCR fixture from edaSimulate(): known onsets, tonic, and phasic.
qp_truth_fixture <- function(seed) {
  pe <- edaSimulate(n_time = 600, sr = 10, scr_count = 5, scr_amplitude = 0.9,
                    noise_sd = 0.003, seed = seed)
  list(pe = pe, truth = S4Vectors::metadata(pe)$eda_truth)
}

# Cluster driver samples into distinct SCR onset times (peaks > threshold).
qp_driver_onsets <- function(driver, sr) {
  n <- length(driver)
  t <- (seq_len(n) - 1) / sr
  thr <- 0.2 * max(driver)
  is_pk <- driver > thr & driver >= c(driver[-1], 0) & driver >= c(0, driver[-n])
  pk_t <- t[is_pk]
  det <- numeric(0)
  for (pt in pk_t) if (length(det) == 0 || min(abs(det - pt)) > 1) det <- c(det, pt)
  det
}

test_that(".bateman_arma and .spline_tonic_basis have the expected structure", {
  ba <- PhysioEDA:::.bateman_arma(n = 50, dt = 0.2, tau_rise = 0.75, tau_decay = 2)
  expect_equal(dim(ba$A), c(50L, 50L))
  expect_equal(dim(ba$M), c(50L, 50L))
  # ARMA is a 3-tap banded system: first two rows are the startup (all zero).
  expect_true(all(ba$A[1:2, ] == 0))
  expect_true(all(ba$M[1:2, ] == 0))
  expect_equal(ba$M[3, 1:3], c(1, 2, 1))          # MA taps on cols i-2..i

  B <- PhysioEDA:::.spline_tonic_basis(n = 200, dt = 0.2, delta_knot = 8)
  expect_equal(nrow(B), 200L)
  expect_gt(ncol(B), 1L)
  expect_true(all(B >= 0))                          # triangular tent functions
})

test_that(".cvxeda_qp_core solves the QP with a non-negative driver", {
  skip_if_not_installed("quadprog")
  for (seed in c(7, 11)) {
    fx <- qp_truth_fixture(seed)
    raw <- SummarizedExperiment::assay(fx$pe, "raw")[, 1]
    # decimate to the QP rate the wrapper would use (sr 10 -> 5 Hz)
    grp <- rep(seq_len(300), each = 2)
    y <- as.numeric(tapply(raw, grp, mean))
    core <- PhysioEDA:::.cvxeda_qp_core(y, sr = 5, alpha = 0.003,
                                        gamma = 0.02, delta_knot = 8)

    expect_length(core$tonic, length(y))
    expect_true(all(core$driver >= 0))              # non-negative driver
    # QP objective is non-increasing relative to the driver-free start
    expect_lte(core$objective, core$obj_init + 1e-6)
    expect_true(core$converged)
  }
})

test_that(".cvxeda_qp reports the QP engine and native-length components", {
  skip_if_not_installed("quadprog")
  fx <- qp_truth_fixture(11)
  raw <- SummarizedExperiment::assay(fx$pe, "raw")[, 1]
  res <- PhysioEDA:::.cvxeda_qp(raw, sr = 10, alpha = 0.003, gamma = 0.02,
                                delta_knot = 8)
  expect_equal(res$engine, "qp")
  expect_equal(res$solver, "quadprog")
  expect_length(res$tonic, length(raw))
  expect_length(res$phasic, length(raw))
  expect_length(res$driver, length(raw))
  expect_true(all(res$driver >= 0))
  expect_true(is.finite(res$objective))
})

test_that(".cvxeda_qp falls back cleanly when no QP solver is available", {
  fx <- qp_truth_fixture(11)
  raw <- SummarizedExperiment::assay(fx$pe, "raw")[, 1]
  # Simulate quadprog being absent.
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "quadprog")) FALSE
      else base::requireNamespace(package, ...)
    }, .package = "base")
  res <- expect_warning(
    PhysioEDA:::.cvxeda_qp(raw, sr = 10, alpha = 0.003, gamma = 0.02),
    "falling back")
  expect_equal(res$engine, "wiener_fallback")
  expect_equal(res$solver, "none")
  expect_length(res$tonic, length(raw))
  expect_length(res$phasic, length(raw))
  expect_true(all(res$driver >= 0))
})
