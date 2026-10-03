library(testthat)
library(PhysioEDA)

test_that(".dda_decompose returns a discrete driver and the expected structure", {
  fx <- make_eda_scr_clean(seed = 1)
  signal <- SummarizedExperiment::assay(fx$pe, "raw")[, 1]
  sr <- samplingRate(fx$pe)
  res <- PhysioEDA:::.dda_decompose(signal, sr, amp_threshold = 0.1,
                                    min_interval_sec = 2)

  expect_type(res, "list")
  expect_true(all(c("tonic", "phasic", "driver", "tau1", "tau2",
                    "onsets", "nSCR") %in% names(res)))
  expect_length(res$tonic, length(signal))
  expect_length(res$phasic, length(signal))
  expect_length(res$driver, length(signal))
  expect_true(all(res$driver >= 0))
  # the driver is a sparse impulse train: far fewer nonzero samples than n
  expect_lt(sum(res$driver > 0), length(signal) / 4)
  expect_equal(res$nSCR, length(res$onsets))
})

test_that("DDA and CDA agree on nSCR and onset timing on clean data", {
  # matched-tau SCRs so both recover the same events (acceptance: agreement
  # within +/-0.5 s, and every true SCR is found)
  for (seed in 1:3) {
    fx <- make_eda_scr_clean(seed = seed)
    signal <- SummarizedExperiment::assay(fx$pe, "raw")[, 1]
    sr <- samplingRate(fx$pe)

    cda <- PhysioEDA:::.cda_decompose(signal, sr)
    thr <- 0.1 * max(pmax(cda$driver, 0))
    mc <- PhysioEDA:::.cda_metrics(cda$driver, sr, amp_threshold = thr,
                                   min_interval_sec = 2)
    dda <- PhysioEDA:::.dda_decompose(signal, sr, amp_threshold = thr,
                                      min_interval_sec = 2)

    # same number of SCRs
    expect_equal(mc$nSCR, dda$nSCR)
    # onset timing agrees within 0.5 s
    expect_true(all(abs(sort(mc$onset_sec) - sort(dda$onsets)) <= 0.5))
    # every true SCR is recovered by CDA
    matched <- vapply(fx$onsets,
                      function(o) any(abs(mc$onset_sec - o) <= 0.5), logical(1))
    expect_true(all(matched))
  }
})

test_that("edaDecompose(method='dda') adds assays and records driver metrics", {
  fx <- make_eda_scr_clean(seed = 2)
  dec <- edaDecompose(fx$pe, method = "dda")

  nms <- SummarizedExperiment::assayNames(dec)
  expect_true(all(c("tonic", "phasic", "driver") %in% nms))

  info <- S4Vectors::metadata(dec)$eda_decompose$cda
  expect_equal(info$method, "dda")
  expect_equal(nrow(info$tau), 1L)
  expect_true(is.list(info$metrics))
  expect_true(info$metrics[[1]]$nSCR >= 0)
})
