library(testthat)
library(PhysioEDA)

# --- edaFilter tests ---

test_that("edaFilter lowpass removes high frequencies", {
  x <- make_eda()
  sr <- samplingRate(x)

  result <- edaFilter(x, type = "lowpass", cutoff = 0.5)

  expect_s4_class(result, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(result))
  expect_equal(dim(SummarizedExperiment::assay(result, "filtered")),
               dim(SummarizedExperiment::assay(x, "raw")))

  # Filtered signal should be smoother (lower variance of diff)
  raw_sig <- SummarizedExperiment::assay(x, "raw")[, 1]
  filt_sig <- SummarizedExperiment::assay(result, "filtered")[, 1]
  expect_true(var(diff(filt_sig)) < var(diff(raw_sig)))
})

test_that("edaFilter highpass removes DC offset", {
  x <- make_eda()

  result <- edaFilter(x, type = "highpass", cutoff = 0.05)

  filt_sig <- SummarizedExperiment::assay(result, "filtered")[, 1]
  raw_sig <- SummarizedExperiment::assay(x, "raw")[, 1]

  # Highpass should make mean closer to zero
  expect_true(abs(mean(filt_sig)) < abs(mean(raw_sig)))
})

test_that("edaFilter bandpass passes intermediate frequencies", {
  x <- make_eda()

  result <- edaFilter(x, type = "bandpass", cutoff = c(0.01, 1.0))

  expect_true("filtered" %in% SummarizedExperiment::assayNames(result))
  filt_sig <- SummarizedExperiment::assay(result, "filtered")[, 1]
  expect_equal(length(filt_sig), nrow(SummarizedExperiment::assay(x, "raw")))
})

test_that("edaFilter stores metadata", {
  x <- make_eda()

  result <- edaFilter(x, type = "lowpass", cutoff = 1.0, order = 3)

  meta <- S4Vectors::metadata(result)$eda_filter
  expect_equal(meta$type, "lowpass")
  expect_equal(meta$cutoff, 1.0)
  expect_equal(meta$order, 3)
})

test_that("edaFilter respects custom output_assay name", {
  x <- make_eda()

  result <- edaFilter(x, type = "lowpass", cutoff = 1.0,
                      output_assay = "my_filtered")

  expect_true("my_filtered" %in% SummarizedExperiment::assayNames(result))
})

test_that("edaFilter validates inputs", {
  x <- make_eda()

  # Cutoff above Nyquist
  expect_error(edaFilter(x, type = "lowpass", cutoff = 10),
               "Nyquist")

  # Bandpass needs length-2 cutoff
 expect_error(edaFilter(x, type = "bandpass", cutoff = 1.0))

  # Bandpass cutoff[1] must be < cutoff[2]
  expect_error(edaFilter(x, type = "bandpass", cutoff = c(2, 1)))

  # Non-PhysioExperiment input
  expect_error(edaFilter(data.frame(), type = "lowpass", cutoff = 1))
})

test_that("edaFilter lowpass preserves slow signal components", {
  # Create a signal with a known slow sine wave
  sr <- 10
  n <- 1000
  t_sec <- seq(0, (n - 1) / sr, length.out = n)

  # Slow wave (0.02 Hz) + fast noise
  slow_wave <- sin(2 * pi * 0.02 * t_sec)
  noisy <- slow_wave + 0.5 * sin(2 * pi * 2 * t_sec)

  mat <- matrix(noisy + 5, ncol = 1)  # offset to keep positive
  x <- PhysioExperiment(
    assays = list(raw = mat),
    colData = S4Vectors::DataFrame(label = "EDA1", type = "EDA", unit = "uS"),
    samplingRate = sr
  )

  result <- edaFilter(x, type = "lowpass", cutoff = 0.5)
  filt_sig <- SummarizedExperiment::assay(result, "filtered")[, 1]

  # The slow wave should be mostly preserved
  # Correlation between filtered and slow_wave + offset should be high
  expect_true(cor(filt_sig, slow_wave + 5) > 0.9)
})

# --- edaDownsample tests ---

test_that("edaDownsample reduces sampling rate and length", {
  x <- make_eda()  # sr=10, n_time=6000
  sr <- samplingRate(x)

  result <- edaDownsample(x, target_sr = 5)

  expect_s4_class(result, "PhysioExperiment")
  expect_equal(samplingRate(result), 5)
  # 6000 / 2 = 3000
  expect_equal(nrow(SummarizedExperiment::assay(result, "raw")), 3000)
})

test_that("edaDownsample preserves channel metadata", {
  x <- make_eda()

  result <- edaDownsample(x, target_sr = 5)

  cd_orig <- SummarizedExperiment::colData(x)
  cd_new <- SummarizedExperiment::colData(result)
  expect_equal(as.character(cd_new$label), as.character(cd_orig$label))
  expect_equal(as.character(cd_new$type), as.character(cd_orig$type))
})

test_that("edaDownsample preserves events", {
  x <- make_eda_event()

  result <- edaDownsample(x, target_sr = 5)

  events_orig <- getEvents(x)
  events_new <- getEvents(result)

  # Event onsets (in seconds) should be preserved
 expect_equal(events_new@events$onset, events_orig@events$onset)
})

test_that("edaDownsample stores downsampling metadata", {
  x <- make_eda()

  result <- edaDownsample(x, target_sr = 5)

  meta <- S4Vectors::metadata(result)$eda_downsample
  expect_equal(meta$original_sr, 10)
  expect_equal(meta$target_sr, 5)
  expect_equal(meta$factor, 2)
  expect_equal(meta$original_n_time, 6000)
  expect_equal(meta$new_n_time, 3000)
})

test_that("edaDownsample validates inputs", {
  x <- make_eda()  # sr=10

  # target_sr must be lower than current sr
  expect_error(edaDownsample(x, target_sr = 10), "lower")
  expect_error(edaDownsample(x, target_sr = 20), "lower")

  # Non-PhysioExperiment input
  expect_error(edaDownsample(data.frame(), target_sr = 5))
})

test_that("edaDownsample applies anti-aliasing filter", {
  # Create a signal with both slow and fast components
  sr <- 100
  n <- 10000
  t_sec <- seq(0, (n - 1) / sr, length.out = n)

  # Slow (1 Hz) + fast (40 Hz that would alias at sr=10)
  slow <- sin(2 * pi * 1 * t_sec)
  fast <- 0.5 * sin(2 * pi * 40 * t_sec)
  signal <- slow + fast + 5

  mat <- matrix(signal, ncol = 1)
  x <- PhysioExperiment(
    assays = list(raw = mat),
    colData = S4Vectors::DataFrame(label = "EDA1", type = "EDA", unit = "uS"),
    samplingRate = sr
  )

  result <- edaDownsample(x, target_sr = 10)

  # The 40 Hz component should be attenuated by the anti-alias filter
  ds_sig <- SummarizedExperiment::assay(result, "raw")[, 1]
  ds_t <- seq(0, by = 1 / 10, length.out = length(ds_sig))
  slow_ref <- sin(2 * pi * 1 * ds_t) + 5

  # The downsampled signal should correlate well with the slow component
  expect_true(cor(ds_sig, slow_ref) > 0.9)
})

test_that("edaDownsample with non-integer factor uses nearest integer", {
  # sr=10, target_sr=3 -> factor=round(10/3)=3, actual_sr=10/3=3.33
  x <- make_eda()

  result <- edaDownsample(x, target_sr = 3)

  # factor should be 3, actual_sr should be 10/3
  meta <- S4Vectors::metadata(result)$eda_downsample
  expect_equal(meta$factor, 3)
  expect_equal(samplingRate(result), 10 / 3)
  expect_equal(nrow(SummarizedExperiment::assay(result, "raw")), 2000)
})
