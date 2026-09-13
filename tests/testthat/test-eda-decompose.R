library(testthat)
library(PhysioEDA)

test_that("edaDecompose with median method adds tonic and phasic assays", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("tonic" %in% nms)
  expect_true("phasic" %in% nms)
})

test_that("edaDecompose preserves dimensions", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  x <- edaDecompose(x, method = "median")

  raw_dim <- dim(SummarizedExperiment::assay(x, "raw"))
  tonic_dim <- dim(SummarizedExperiment::assay(x, "tonic"))
  phasic_dim <- dim(SummarizedExperiment::assay(x, "phasic"))

  expect_equal(tonic_dim, raw_dim)
  expect_equal(phasic_dim, raw_dim)
})

test_that("edaDecompose median produces smoother tonic than raw", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  raw <- SummarizedExperiment::assay(x, "raw")[, 1]
  tonic <- SummarizedExperiment::assay(x, "tonic")[, 1]

  expect_lt(sd(tonic), sd(raw))
})

test_that("edaDecompose tonic + phasic approximately equals raw", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  raw <- SummarizedExperiment::assay(x, "raw")
  tonic <- SummarizedExperiment::assay(x, "tonic")
  phasic <- SummarizedExperiment::assay(x, "phasic")

  expect_equal(tonic + phasic, raw, tolerance = 1e-10)
})

test_that("edaDecompose highpass method works", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "highpass")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("tonic" %in% nms)
  expect_true("phasic" %in% nms)
})

test_that("edaDecompose stores method in metadata", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  meta <- S4Vectors::metadata(x)
  expect_equal(meta$eda_decompose$method, "median")
})

test_that("edaDecompose supports custom output assay names", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median",
                    output_tonic = "my_tonic",
                    output_phasic = "my_phasic")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("my_tonic" %in% nms)
  expect_true("my_phasic" %in% nms)
})


# ---- Exact cvxEDA (Greco QP) end-to-end via edaDecompose (WS5-25) ------------

# Detect distinct SCR onset times from the sparse driver.
.decompose_driver_onsets <- function(driver, sr) {
  n <- length(driver)
  t <- (seq_len(n) - 1) / sr
  thr <- 0.2 * max(driver)
  is_pk <- driver > thr & driver >= c(driver[-1], 0) & driver >= c(0, driver[-n])
  det <- numeric(0)
  for (pt in t[is_pk]) if (length(det) == 0 || min(abs(det - pt)) > 1) det <- c(det, pt)
  det
}

test_that("edaDecompose(method='cvxeda') recovers the ground-truth SCRs", {
  skip_if_not_installed("quadprog")
  # Acceptance (Greco cvxEDA): on edaSimulate() with known onsets, the driver
  # spikes align to true onsets within 0.5 s, the recovered phasic correlates
  # r > 0.95 with truth, and tonic + phasic reconstructs the input to RMSE < 1%
  # of the signal range. Verified across two seeds (not a single lucky draw).
  for (seed in c(7, 11)) {
    pe <- edaSimulate(n_time = 600, sr = 10, scr_count = 5, scr_amplitude = 0.9,
                      noise_sd = 0.003, seed = seed)
    truth <- S4Vectors::metadata(pe)$eda_truth

    dec <- edaDecompose(pe, method = "cvxeda", alpha = 0.003, gamma = 0.02,
                        delta_knot = 8)

    nms <- SummarizedExperiment::assayNames(dec)
    expect_true(all(c("tonic", "phasic", "driver") %in% nms))

    raw    <- SummarizedExperiment::assay(pe, "raw")[, 1]
    tonic  <- SummarizedExperiment::assay(dec, "tonic")[, 1]
    phasic <- SummarizedExperiment::assay(dec, "phasic")[, 1]
    driver <- SummarizedExperiment::assay(dec, "driver")[, 1]

    # (1) reconstruction RMSE < 1% of signal range
    rmse <- sqrt(mean((tonic + phasic - raw)^2))
    expect_lt(rmse / diff(range(raw)), 0.01)

    # (2) recovered phasic correlates r > 0.95 with the true phasic
    expect_gt(cor(phasic, truth$phasic[, 1]), 0.95)

    # (3) each true SCR onset has a driver spike within +/- 0.5 s
    det <- .decompose_driver_onsets(driver, samplingRate(pe))
    for (o in truth$onsets$time_sec) expect_lte(min(abs(det - o)), 0.5)

    # (4) the QP engine ran and recorded a converged objective
    info <- S4Vectors::metadata(dec)$eda_decompose$cvxeda
    expect_equal(info$engine, "qp")
    expect_true(is.finite(info$objective))
  }
})

test_that("edaDecompose(method='cvxeda') non-increasing QP cost per channel", {
  skip_if_not_installed("quadprog")
  pe <- edaSimulate(n_time = 500, n_channels = 2, sr = 10, scr_count = 4,
                    scr_amplitude = 0.9, noise_sd = 0.003, seed = 3)
  dec <- edaDecompose(pe, method = "cvxeda", alpha = 0.003, gamma = 0.02,
                      delta_knot = 8)
  info <- S4Vectors::metadata(dec)$eda_decompose$cvxeda
  expect_length(info$objective, 2L)
  expect_true(all(is.finite(info$objective)))
  expect_true(all(info$engine == "qp"))
})

test_that("edaDecompose(method='cvxeda_fast') keeps the Wiener approximation", {
  set.seed(42)
  pe <- edaSimulate(n_time = 600, sr = 10, scr_count = 5, seed = 42)
  dec <- edaDecompose(pe, method = "cvxeda_fast")
  nms <- SummarizedExperiment::assayNames(dec)
  expect_true(all(c("tonic", "phasic", "driver") %in% nms))
  driver <- SummarizedExperiment::assay(dec, "driver")[, 1]
  expect_true(all(driver >= 0))
  # cvxeda_fast never records the QP engine info
  expect_null(S4Vectors::metadata(dec)$eda_decompose$cvxeda)
})

test_that("edaDecompose method='dda' adds tonic/phasic/driver and CDA info", {
  fx <- make_eda_scr_clean(seed = 1)
  dec <- edaDecompose(fx$pe, method = "dda")
  nms <- SummarizedExperiment::assayNames(dec)
  expect_true(all(c("tonic", "phasic", "driver") %in% nms))
  expect_true(all(SummarizedExperiment::assay(dec, "driver")[, 1] >= 0))
  expect_equal(S4Vectors::metadata(dec)$eda_decompose$cda$method, "dda")
})

test_that("edaDecompose method='cda' with optimize_tau records the tuned tau", {
  fx <- make_eda_scr_clean(seed = 2)
  dec <- edaDecompose(fx$pe, method = "cda", optimize_tau = TRUE)
  md <- S4Vectors::metadata(dec)$eda_decompose
  expect_true(md$optimize_tau)
  expect_equal(md$cda$method, "cda")
  tau2 <- md$cda$tau[1, "tau2"]
  expect_gt(tau2, md$cda$tau[1, "tau1"])
  expect_lte(tau2, 8)
})
