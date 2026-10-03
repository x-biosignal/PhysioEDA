test_that("ns_scr_freq counts only amplitude-qualifying non-specific SCRs (< scr_rate)", {
  x <- make_eda(n_time = 6000, n_channels = 1, sr = 10)
  lab <- SummarizedExperiment::colData(x)$label
  if (is.null(lab)) lab <- "Ch1"
  lab <- lab[1]

  peaks <- data.frame(
    channel = lab,
    onset_sec = c(10, 60, 120, 200, 300),
    amplitude = c(0.10, 0.03, 0.20, 0.04, 0.50),  # 3 >= 0.05, 2 < 0.05
    stringsAsFactors = FALSE
  )
  feat <- edaFeatures(x, peaks = peaks)
  row <- feat[feat$channel == lab, ]

  duration_min <- (6000 - 1) / 10 / 60
  expect_equal(row$scr_count, 5)
  expect_equal(row$ns_scr_count, 3)
  expect_equal(row$ns_scr_freq, 3 / duration_min, tolerance = 1e-9)
  # non-specific frequency must be strictly below the overall SCR rate here
  expect_lt(row$ns_scr_freq, row$scr_rate_per_min)
  # regression guard against the old bug: ns_scr_freq == scr_rate_per_min
  expect_false(isTRUE(all.equal(row$ns_scr_freq, row$scr_rate_per_min)))
})
