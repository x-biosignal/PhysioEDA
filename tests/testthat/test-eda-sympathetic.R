library(testthat)
library(PhysioEDA)

# single-channel EDA PhysioExperiment carrying a controlled oscillation
mk_pe <- function(sig, fs) {
  PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(eda = matrix(sig, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "EDA", type = "EDA"),
    samplingRate = fs)
}
osc_pe <- function(A, fs = 4, dur = 300, seed = 1, f = 0.1) {
  set.seed(seed); n <- dur * fs; t <- (0:(n - 1)) / fs
  mk_pe(5 + 0.05 * rnorm(n) + A * sin(2 * pi * f * t), fs)
}

test_that("EDASymp increases monotonically with sympathetic-band amplitude", {
  amps <- c(0, 0.1, 0.3, 0.6)
  vals <- vapply(amps, function(A) edaSympatheticIndex(osc_pe(A))$edasymp,
                 numeric(1))
  expect_true(all(diff(vals) > 0))                       # strictly increasing
  norms <- vapply(amps, function(A) edaSympatheticIndex(osc_pe(A))$edasymp_norm,
                  numeric(1))
  expect_true(all(norms >= 0 & norms <= 1))              # bounded
  expect_gt(norms[length(norms)], norms[1])
})

test_that("EDASymp returns one row per channel and validates inputs", {
  pe <- osc_pe(0.3)
  r <- edaSympatheticIndex(pe)
  expect_s3_class(r, "data.frame")
  expect_equal(nrow(r), 1L)
  expect_named(r, c("channel", "edasymp", "edasymp_norm", "total_power",
                    "band_low", "band_high"))
  expect_error(edaSympatheticIndex(pe, band = c(0.25, 0.045)), "low < high")
  expect_error(edaSympatheticIndex(pe, band = c(0.1, 5), resample_hz = NULL),
               "Nyquist")
  expect_error(edaSympatheticIndex(list()), "PhysioExperiment")
})

test_that("TVSymp trajectory peaks when a transient oscillation appears", {
  fs <- 4; dur <- 300; n <- dur * fs; t <- (0:(n - 1)) / fs
  set.seed(3)
  env <- as.numeric(t >= 120 & t <= 180)                 # transient 0.1 Hz burst
  sig <- 5 + 0.05 * rnorm(n) + 0.5 * env * sin(2 * pi * 0.1 * t)
  tv <- edaTVSymp(mk_pe(sig, fs))
  peak_t <- tv$time_sec[which.max(tv$tvsymp)]
  expect_gt(peak_t, 110)
  expect_lt(peak_t, 190)                                 # peak inside the burst
  in_win <- mean(tv$tvsymp[tv$time_sec >= 120 & tv$time_sec <= 180])
  out_win <- mean(tv$tvsymp[tv$time_sec < 100 | tv$time_sec > 200])
  expect_gt(in_win, 3 * out_win)
})

test_that("TVSymp is finite, bounded and reproducible", {
  pe <- osc_pe(0.4)
  tv1 <- edaTVSymp(pe); tv2 <- edaTVSymp(pe)
  expect_s3_class(tv1, "data.frame")
  expect_true(all(is.finite(tv1$tvsymp)))
  expect_true(all(is.finite(tv1$tvsymp_norm)))
  expect_identical(tv1$tvsymp, tv2$tvsymp)               # deterministic
  expect_equal(nrow(tv1), 600L)                          # resampled to 2 Hz
  expect_length(attr(tv1, "components")[[1]], 2L)
  expect_error(edaTVSymp(pe, band = c(0.3, 0.1)), "low < high")
})

# --- regression tests for adversarial-review findings (WS5-30) ----------------

test_that("short and single-sample signals do not crash the Welch estimator", {
  # a very short (few-second) record must not throw the opaque seq() 'by' error
  fs <- 4; t <- (0:39) / fs                         # 10 s -> resamples to ~20 @2Hz
  pe <- mk_pe(5 + 0.05 * sin(2 * pi * 0.1 * t), fs)
  expect_s3_class(edaSympatheticIndex(pe), "data.frame")
  expect_s3_class(edaTVSymp(pe), "data.frame")
  # even a sub-second record stays finite instead of erroring internally
  short <- mk_pe(c(5, 5.1, 5.0, 4.9, 5.05, 5.1, 5.0, 4.95), fs)
  r <- edaSympatheticIndex(short)
  expect_true(is.finite(r$edasymp))
  # a single time point is rejected with a clear message
  one <- mk_pe(5, fs)
  expect_error(edaSympatheticIndex(one), "time points")
  expect_error(edaTVSymp(one), "time points")
})

test_that("an empty (zero-channel) assay is rejected, not returned as NULL", {
  m <- matrix(numeric(0), nrow = 100, ncol = 0)
  pe <- PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(eda = m),
    colData = S4Vectors::DataFrame(label = character(0), type = character(0)),
    samplingRate = 4)
  expect_error(edaSympatheticIndex(pe), "channel")
  expect_error(edaTVSymp(pe), "channel")
})

test_that("a constant channel yields zero band power and NA normalised index", {
  pe <- mk_pe(rep(5, 400), 4)
  r <- edaSympatheticIndex(pe)
  expect_equal(r$edasymp, 0)
  expect_true(is.na(r$edasymp_norm))                # documented NA fallback
  expect_true(all(is.finite(edaTVSymp(pe)$tvsymp))) # constant -> finite (zero)
})

test_that("odd-length Welch segments recover the correct top-bin power", {
  # a unit sine near the top of an odd-length grid must not be halved
  fs <- 2; n <- 101; t <- (0:(n - 1)) / fs
  x <- sin(2 * pi * (fs * (n %/% 2) / n) * t)       # top odd-grid frequency
  sp <- PhysioEDA:::.eda_welch(x - mean(x), fs, seg_len = n)
  total <- PhysioEDA:::.eda_bandpower(sp$freq, sp$psd, sp$freq[2], max(sp$freq))
  expect_gt(total, 0.4)                             # ~0.5, not ~0.27 (halved)
})

test_that("both indices handle a multi-channel object", {
  fs <- 4; dur <- 200; n <- dur * fs; t <- (0:(n - 1)) / fs
  set.seed(5)
  m <- cbind(5 + 0.05 * rnorm(n) + 0.4 * sin(2 * pi * 0.1 * t),
             5 + 0.05 * rnorm(n))
  pe <- PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(eda = m),
    colData = S4Vectors::DataFrame(label = c("E1", "E2"), type = "EDA"),
    samplingRate = fs)
  r <- edaSympatheticIndex(pe)
  expect_equal(nrow(r), 2L)
  expect_gt(r$edasymp[1], r$edasymp[2])                  # ch1 carries the oscillation
  tv <- edaTVSymp(pe)
  expect_setequal(unique(tv$channel), c(1, 2))
})
