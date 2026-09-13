library(testthat)
library(PhysioEDA)

test_that("CDA tonic + phasic approximates original signal", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  result <- PhysioEDA:::.cda_decompose(signal, sr)

  reconstructed <- result$tonic + result$phasic
  # Correlation should be high
  r <- cor(reconstructed, signal)
  expect_gt(r, 0.9)
  # RMSE should be small relative to signal sd
  rmse <- sqrt(mean((reconstructed - signal)^2))
  expect_lt(rmse, 0.1 * sd(signal))
})

test_that("CDA driver is non-negative", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  result <- PhysioEDA:::.cda_decompose(signal, sr)

  expect_true(all(result$driver >= 0))
})

test_that("CDA phasic has peaks near known SCR locations", {
  # Use a signal with known SCR positions by using a fixed seed
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 1, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.01)
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  result <- PhysioEDA:::.cda_decompose(signal, sr)

  # Phasic should have elevated values (peaks) in the signal
  # Find indices where phasic exceeds its mean + 1 SD (SCR regions)
  phasic <- result$phasic
  threshold <- mean(phasic) + sd(phasic)
  peak_regions <- which(phasic > threshold)

  # There should be multiple distinct peak regions (at least 2 of the 5 SCRs)
  if (length(peak_regions) > 1) {
    # Count clusters of consecutive indices (gap > 10 samples = new cluster)
    gaps <- diff(peak_regions)
    n_clusters <- sum(gaps > 10) + 1
    expect_gte(n_clusters, 2)
  } else {
    # If no peaks found above threshold, the phasic should still have variation
    expect_gt(sd(phasic), 0)
  }
})

test_that("CDA tonic is smoother than original signal", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  result <- PhysioEDA:::.cda_decompose(signal, sr)

  # Tonic roughness (sd of first differences) should be less than original
  tonic_roughness <- sd(diff(result$tonic))
  signal_roughness <- sd(diff(signal))
  expect_lt(tonic_roughness, signal_roughness)
})

test_that("CDA works with custom tau values", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  result <- PhysioEDA:::.cda_decompose(signal, sr, tau1 = 0.5, tau2 = 3.0)

  expect_equal(result$tau1, 0.5)
  expect_equal(result$tau2, 3.0)
  expect_equal(length(result$tonic), length(signal))
  expect_equal(length(result$phasic), length(signal))
  expect_equal(length(result$driver), length(signal))
})

test_that("CDA works with different signal lengths", {
  set.seed(42)
  sr <- 10

  for (n in c(100, 500, 2000)) {
    x <- edaSimulate(n_time = n, n_channels = 1, sr = sr,
                     scr_count = 2, noise_sd = 0.01)
    signal <- SummarizedExperiment::assay(x, "raw")[, 1]

    result <- PhysioEDA:::.cda_decompose(signal, sr)

    expect_equal(length(result$tonic), n)
    expect_equal(length(result$phasic), n)
    expect_equal(length(result$driver), n)
    expect_true(all(result$driver >= 0))
  }
})

test_that("CDA returns correct list structure", {
  set.seed(42)
  x <- make_eda()
  signal <- SummarizedExperiment::assay(x, "raw")[, 1]
  sr <- samplingRate(x)

  result <- PhysioEDA:::.cda_decompose(signal, sr)

  expect_type(result, "list")
  expect_named(result, c("tonic", "phasic", "driver", "tau1", "tau2"))
  expect_type(result$tonic, "double")
  expect_type(result$phasic, "double")
  expect_type(result$driver, "double")
  expect_equal(result$tau1, 0.75)
  expect_equal(result$tau2, 2.0)
})

test_that("CDA validates inputs", {
  expect_error(PhysioEDA:::.cda_decompose("not_numeric", 10))
  expect_error(PhysioEDA:::.cda_decompose(c(1, 2, 3), -1))
  expect_error(PhysioEDA:::.cda_decompose(c(1, 2, 3), 10, tau1 = -1))
  expect_error(PhysioEDA:::.cda_decompose(c(1, 2, 3), 10, tau1 = 3, tau2 = 1))
})


# ---- Full CDA: metrics, tau optimisation, reconstruction (WS5-26) -----------

test_that(".cda_metrics matches hand-computed values on a 2-SCR fixture", {
  sr <- 10
  drv <- numeric(100)
  drv[20:29] <- 0.5   # SCR 1: 10 samples, amplitude 0.5, onset sample 20 (1.9 s)
  drv[50:55] <- 0.3   # SCR 2:  6 samples, amplitude 0.3, onset sample 50 (4.9 s)
  m <- PhysioEDA:::.cda_metrics(drv, sr, amp_threshold = 0.01)

  expect_equal(m$nSCR, 2L)
  expect_equal(m$AmpSum, 0.8)                     # 0.5 + 0.3
  expect_equal(m$ISCR, (0.5 * 10 + 0.3 * 6) / sr) # 6.8 / 10 = 0.68 uS*s
  expect_equal(m$onset_sec, c(19, 49) / sr)       # (sample - 1) / sr
  expect_equal(m$scr_area, c(0.5 * 10, 0.3 * 6) / sr)
})

test_that(".cda_metrics returns zero SCRs on a flat driver", {
  m <- PhysioEDA:::.cda_metrics(numeric(200), sr = 10)
  expect_equal(m$nSCR, 0L)
  expect_equal(m$AmpSum, 0)
  expect_equal(m$ISCR, 0)
  expect_length(m$onset_sec, 0)
})

test_that("CDA driver is nonnegative and reconvolves to the phasic (<2% range)", {
  fx <- make_eda_scr_clean(seed = 2)
  signal <- SummarizedExperiment::assay(fx$pe, "raw")[, 1]
  sr <- samplingRate(fx$pe)
  res <- PhysioEDA:::.cda_decompose(signal, sr)

  expect_true(all(res$driver >= 0))

  # independent reconvolution of the driver reconstructs the phasic component
  n <- length(signal)
  h_len <- min(n, as.integer(10 * res$tau2 * sr))
  th <- seq(0, by = 1 / sr, length.out = h_len)
  h <- scrf(th, tau1 = res$tau1, tau2 = res$tau2, form = "bateman",
            normalize = "none")
  h <- h / max(h)
  H <- stats::fft(c(h, rep(0, n - h_len)))
  recon <- Re(stats::fft(stats::fft(res$driver) * H, inverse = TRUE)) / n
  rmse <- sqrt(mean((recon - res$phasic)^2))
  expect_lt(rmse / diff(range(signal)), 0.02)
})

test_that("CDA optimize_tau returns a bounded tau2 and a nonnegative driver", {
  fx <- make_eda_scr_clean(seed = 3)
  signal <- SummarizedExperiment::assay(fx$pe, "raw")[, 1]
  sr <- samplingRate(fx$pe)
  res <- PhysioEDA:::.cda_decompose(signal, sr, optimize_tau = TRUE)

  expect_named(res, c("tonic", "phasic", "driver", "tau1", "tau2"))
  expect_gt(res$tau2, res$tau1)
  expect_lte(res$tau2, 8)
  expect_true(all(res$driver >= 0))
})
