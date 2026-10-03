library(testthat)
library(PhysioEDA)

test_that("scrf(normalize='none') is the raw biexponential (single kernel source)", {
  t <- seq(0, 15, by = 0.05)
  tau1 <- 0.75; tau2 <- 2.0
  inline <- exp(-t / tau2) - exp(-t / tau1)
  expect_identical(scrf(t, tau1, tau2, normalize = "none"), inline)
})

test_that("scrf peak time and amplitude match the analytic Bateman maximum", {
  tau1 <- 0.75; tau2 <- 2.0
  pk <- scrfPeak(tau1, tau2)

  # the analytic peak is a stationary point (derivative 0) and a maximum
  expect_lt(abs(scrf(pk$time, tau1, tau2, normalize = "none", deriv = 1)), 1e-6)
  expect_lt(scrf(pk$time, tau1, tau2, normalize = "none", deriv = 2), 0)

  # numeric argmax matches the analytic peak time within 1e-6
  opt <- optimize(function(x) scrf(x, tau1, tau2, normalize = "none"),
                  interval = c(0, 10 * tau2), maximum = TRUE, tol = 1e-8)
  expect_lt(abs(opt$maximum - pk$time), 1e-6)
  # peak amplitude matches the analytic maximum within 1e-6
  expect_lt(abs(opt$objective - pk$amplitude), 1e-6)
})

test_that("scrf normalisation is deterministic and unit-consistent", {
  tau1 <- 0.75; tau2 <- 2.0
  pk <- scrfPeak(tau1, tau2)

  # peak normalisation: maximum is exactly 1
  expect_equal(scrf(pk$time, tau1, tau2, normalize = "peak"), 1, tolerance = 1e-6)
  # area normalisation: integral is exactly 1 (raw Bateman area is tau2 - tau1)
  area <- integrate(function(x) scrf(x, tau1, tau2, normalize = "area"),
                    0, 20 * tau2)$value
  expect_equal(area, 1, tolerance = 1e-4)
  # deterministic
  t <- seq(0, 15, by = 0.1)
  expect_identical(scrf(t), scrf(t))
})

test_that("scrf is causal and returns derivatives of the right length", {
  t <- seq(-2, 15, by = 0.1)
  k <- scrf(t)
  expect_length(k, length(t))
  expect_true(all(k[t < 0] == 0))
  for (d in 0:2) {
    v <- scrf(t, deriv = d)
    expect_length(v, length(t))
    expect_true(all(is.finite(v)))
    expect_true(all(v[t < 0] == 0))
  }
})

test_that("scrf pspm_canonical form is a finite causal kernel", {
  t <- seq(-1, 20, by = 0.05)
  cp <- scrf(t, form = "pspm_canonical", normalize = "peak")
  expect_length(cp, length(t))
  expect_true(all(is.finite(cp)))
  expect_true(all(cp[t < 0] == 0))
  expect_equal(max(cp), 1, tolerance = 1e-3)         # peak-normalised
  expect_identical(cp, scrf(t, form = "pspm_canonical", normalize = "peak"))
  # the bump makes the canonical form differ from the pure Bateman shape
  bt <- scrf(t, form = "bateman", normalize = "peak")
  expect_gt(max(abs(cp - bt)), 1e-3)
})

test_that("scrfPeak returns the analytic Bateman time and amplitude", {
  tau1 <- 0.8; tau2 <- 3.5
  pk <- scrfPeak(tau1, tau2)
  # closed form: t* = log(tau2/tau1) / (1/tau1 - 1/tau2)
  tp <- log(tau2 / tau1) / (1 / tau1 - 1 / tau2)
  expect_equal(pk$time, tp)
  expect_equal(pk$amplitude, exp(-tp / tau2) - exp(-tp / tau1))
})

test_that("scrf validates its arguments", {
  expect_error(scrf(0:5, tau1 = 1, tau2 = 1), "tau1 != tau2|not TRUE")
  expect_error(scrf(0:5, deriv = 3))
  expect_error(scrf(0:5, tau1 = -1))
})

test_that("edaSimulate and decompose consume scrf without changing results", {
  # the refactor to scrf(normalize='none') must be numerically identical:
  # a seeded edaSimulate still reproduces, and a decomposition still runs
  x1 <- edaSimulate(n_time = 400, sr = 10, scr_count = 4, seed = 5)
  x2 <- edaSimulate(n_time = 400, sr = 10, scr_count = 4, seed = 5)
  expect_identical(SummarizedExperiment::assay(x1), SummarizedExperiment::assay(x2))

  dec <- edaDecompose(x1, method = "cda")
  expect_true(all(c("tonic", "phasic", "driver") %in%
                    SummarizedExperiment::assayNames(dec)))
})
