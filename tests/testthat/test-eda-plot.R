library(testthat)
library(PhysioEDA)

test_that("plotEda runs without error for single channel", {
  x <- make_eda(n_channels = 1)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotEda(x))
})

test_that("plotEda runs without error for multiple channels", {
  x <- make_eda(n_channels = 3)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotEda(x))
})

test_that("plotEda respects time_range parameter", {
  x <- make_eda(n_channels = 1)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotEda(x, time_range = c(10, 50)))
})

test_that("plotEda displays events when present", {
  x <- make_eda_event(n_channels = 1)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotEda(x, show_events = TRUE))
})

test_that("plotEda respects channel selection by index", {
  x <- make_eda(n_channels = 3)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotEda(x, channels = c(1, 3)))
})

test_that("plotEda respects channel selection by name", {
  x <- make_eda(n_channels = 3)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotEda(x, channels = "EDA2"))
})

test_that("plotEda with custom colors", {
  x <- make_eda(n_channels = 2)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotEda(x, col = c("darkblue", "darkred")))
})

test_that("plotDecompose runs without error", {
  x <- make_eda(n_channels = 1)
  x <- edaDecompose(x)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotDecompose(x))
})

test_that("plotDecompose errors without tonic/phasic assays", {
  x <- make_eda(n_channels = 1)
  expect_error(plotDecompose(x), "tonic.*phasic")
})

test_that("plotDecompose respects time_range", {
  x <- make_eda(n_channels = 1)
  x <- edaDecompose(x)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotDecompose(x, time_range = c(20, 80)))
})

test_that("plotDecompose works with channel by name", {
  x <- make_eda(n_channels = 2)
  x <- edaDecompose(x)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotDecompose(x, channel = "EDA2"))
})

test_that("plotPeaks runs without error", {
  x <- make_eda(n_channels = 1)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotPeaks(x))
})

test_that("plotPeaks works with pre-computed peaks", {
  x <- make_eda(n_channels = 1)
  peaks <- edaPeaks(x)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotPeaks(x, peaks = peaks))
})

test_that("plotPeaks respects time_range", {
  x <- make_eda(n_channels = 1)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotPeaks(x, time_range = c(50, 200)))
})

test_that("plotPeaks works after decomposition", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 1, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.001)
  x <- edaDecompose(x)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotPeaks(x))
})

test_that("plotPeaks works with multi-channel selecting channel 2", {
  set.seed(42)
  x <- edaSimulate(n_time = 6000, n_channels = 2, sr = 10,
                   scr_count = 5, scr_amplitude = 1.0, noise_sd = 0.001)
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expect_no_error(plotPeaks(x, channel = 2))
})

test_that("all plot functions return invisible NULL", {
  x <- make_eda(n_channels = 1)
  x_dec <- edaDecompose(x)

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off())

  result1 <- plotEda(x)
  result2 <- plotDecompose(x_dec)
  result3 <- plotPeaks(x)

  expect_null(result1)
  expect_null(result2)
  expect_null(result3)
})
