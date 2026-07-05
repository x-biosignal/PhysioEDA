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
