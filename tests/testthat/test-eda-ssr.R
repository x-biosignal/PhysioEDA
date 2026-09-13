library(testthat)
library(PhysioEDA)

mk_ssr_pe <- function(sig, sr = 10, labels = "EDA") {
  m <- if (is.matrix(sig)) sig else matrix(sig, ncol = 1)
  PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(eda = m),
    colData = S4Vectors::DataFrame(label = labels, type = "EDA"),
    samplingRate = sr)
}
# add a Bateman (biexponential) SCR with peak amplitude `amp` at onset `onset_s`
add_scr <- function(sig, onset_s, amp, sr = 10, tr = 0.75, td = 3) {
  t <- (seq_along(sig) - 1) / sr - onset_s
  k <- exp(-t / td) - exp(-t / tr); k[t <= 0] <- 0
  peak_norm <- max(exp(-(0:2000) / 10 / td) - exp(-(0:2000) / 10 / tr))
  sig + amp * k / peak_norm
}

test_that("SSR latency and amplitude match a synthetic fixture", {
  sr <- 10; n <- 1000; sig <- rep(5, n)
  stim <- c(10, 40, 70); lat <- 1.5; A <- 0.5
  for (st in stim) sig <- add_scr(sig, st + lat, A, sr)
  res <- edaSSR(mk_ssr_pe(sig, sr), stim, events_unit = "seconds")
  expect_s3_class(res, "eda_ssr")
  pr <- res$responses
  expect_true(all(pr$present))
  expect_true(all(abs(pr$latency - lat) <= 0.2))       # latency within +/-0.2 s
  expect_true(all(abs(pr$amplitude - A) / A < 0.05))   # amplitude within 5%
  expect_true(all(pr$morphology == "monophasic"))
})

test_that("a stimulus with no following SCR is classified absent", {
  sr <- 10; n <- 800; sig <- rep(5, n) + 0.02 * sin(2 * pi * 0.01 * (0:(n - 1)) / sr)
  sig <- add_scr(sig, 10 + 1.5, 0.5, sr)               # only the first stimulus
  res <- edaSSR(mk_ssr_pe(sig, sr), c(10, 40), events_unit = "seconds")
  expect_true(res$responses$present[1])
  expect_false(res$responses$present[2])
  expect_true(is.na(res$responses$latency[2]))
  expect_equal(res$summary$response_rate, 0.5)
  expect_false(res$summary$ssr_absent)
  # an entirely flat trace -> clinical absent-SSR flag
  flat <- edaSSR(mk_ssr_pe(rep(5, n), sr), c(10, 40), events_unit = "seconds")
  expect_true(flat$summary$ssr_absent)
  expect_equal(flat$summary$response_rate, 0)
})

test_that("habituation slope is negative for declining responses", {
  sr <- 10; n <- 1600; sig <- rep(5, n)
  amps <- c(0.7, 0.5, 0.35, 0.25); stim <- c(10, 40, 70, 100)
  for (i in seq_along(stim)) sig <- add_scr(sig, stim[i] + 1.5, amps[i], sr)
  res <- edaSSR(mk_ssr_pe(sig, sr), stim, events_unit = "seconds")
  expect_true(all(res$responses$present))
  expect_lt(res$summary$habituation_slope, 0)          # amplitudes decline
  expect_equal(res$summary$response_rate, 1)
})

test_that("a response that dips below baseline is detected as biphasic", {
  sr <- 10; n <- 600; sig <- rep(5, n)
  sig <- add_scr(sig, 10 + 1.5, 0.5, sr)               # positive phase
  # a later negative deflection dipping below baseline
  tb <- (seq_along(sig) - 1) / sr - (10 + 3.5)
  k <- exp(-tb / 1.5) - exp(-tb / 0.5); k[tb <= 0] <- 0
  sig <- sig - 0.7 * k / max(k)
  res <- edaSSR(mk_ssr_pe(sig, sr), 10, events_unit = "seconds")
  expect_true(res$responses$present)
  expect_equal(res$responses$morphology, "biphasic")
})

# --- regression tests for adversarial-review findings (WS5-31) ----------------

test_that("monotone tonic drift is NOT counted as an SSR (absent stays absent)", {
  sr <- 10; n <- 1000; t <- (0:(n - 1)) / sr
  drift <- 5 + 0.08 * t                                # slow rising SCL, no SCR
  res <- edaSSR(mk_ssr_pe(drift, sr), c(10, 40, 70), events_unit = "seconds")
  expect_false(any(res$responses$present))
  expect_true(res$summary$ssr_absent)                  # clinical flag correct
  expect_equal(res$summary$response_rate, 0)
})

test_that("a monophasic SCR on a declining tonic is not misclassified biphasic", {
  sr <- 10; n <- 800; t <- (0:(n - 1)) / sr
  sig <- 5 - 0.06 * t                                  # post-arousal decline
  sig <- add_scr(sig, 10 + 1.5, 0.5, sr)               # clean monophasic SCR
  res <- edaSSR(mk_ssr_pe(sig, sr), 10, events_unit = "seconds")
  expect_true(res$responses$present)
  expect_equal(res$responses$morphology, "monophasic")
  expect_lt(abs(res$responses$amplitude - 0.5) / 0.5, 0.05)
  expect_lte(abs(res$responses$latency - 1.5), 0.2)
})

test_that("habituation slope is chronological, independent of event order", {
  sr <- 10; n <- 1600; sig <- rep(5, n)
  amps <- c(0.7, 0.5, 0.35, 0.25); stim <- c(10, 40, 70, 100)
  for (i in seq_along(stim)) sig <- add_scr(sig, stim[i] + 1.5, amps[i], sr)
  in_order <- edaSSR(mk_ssr_pe(sig, sr), stim, events_unit = "seconds")
  shuffled <- edaSSR(mk_ssr_pe(sig, sr), c(100, 10, 70, 40),
                     events_unit = "seconds")
  expect_equal(in_order$summary$habituation_slope,
               shuffled$summary$habituation_slope)
  expect_lt(in_order$summary$habituation_slope, 0)
})

test_that("stimulus events resolving outside the record are rejected", {
  pe <- mk_ssr_pe(rep(5, 200), 10)
  expect_error(edaSSR(pe, 999, events_unit = "index"), "outside")
  expect_error(edaSSR(pe, 0, events_unit = "index"), "outside")
  expect_error(edaSSR(pe, 10000, events_unit = "seconds"), "outside")
})

test_that("edaSSR validates inputs and handles units and channels", {
  sr <- 10; n <- 800; sig <- add_scr(rep(5, n), 10 + 1.5, 0.5, sr)
  pe <- mk_ssr_pe(sig, sr)
  expect_error(edaSSR(pe, numeric(0)), "non-empty")
  expect_error(edaSSR(pe, c(10, NA)), "finite")
  expect_error(edaSSR(pe, 10, response_window = c(5, 1)), "start < end")
  expect_error(edaSSR(list(), 10), "PhysioExperiment")
  # index and seconds units agree
  r_idx <- edaSSR(pe, round(10 * sr) + 1L, events_unit = "index")
  r_sec <- edaSSR(pe, 10, events_unit = "seconds")
  expect_equal(r_idx$responses$present, r_sec$responses$present)
  expect_equal(r_idx$responses$latency, r_sec$responses$latency)
  # bilateral recording: the second channel is flat -> absent
  pe2 <- mk_ssr_pe(cbind(sig, rep(5, n)), sr, labels = c("L", "R"))
  res <- edaSSR(pe2, 10, events_unit = "seconds")
  expect_equal(nrow(res$summary), 2L)
  expect_true(res$summary$ssr_absent[2])
  expect_output(print(res), "eda_ssr")
})
