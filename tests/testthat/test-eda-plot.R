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

# --- ER-SCR grand-average plot (WS5-32) ---------------------------------------

test_that("grand-average mean equals the arithmetic mean of the epochs", {
  x <- make_eda_with_scr()
  ga <- plotErscrGrandAverage(x, plot = FALSE)
  expect_s3_class(ga, "erscr_grand_average")
  g <- ga$grand_average[ga$grand_average$type == "stimulus", ]
  ep <- ga$epochs[["stimulus"]]
  expect_equal(g$mean, rowMeans(ep), tolerance = 1e-6)   # AC1
  expect_equal(nrow(g), length(ga$time))
})

test_that("CI band width scales as 1/sqrt(n_trials)", {
  x <- make_eda_with_scr()
  ga <- plotErscrGrandAverage(x, plot = FALSE, ci_level = 0.95)
  g <- ga$grand_average[ga$grand_average$type == "stimulus", ]
  z <- stats::qnorm(0.975); n <- g$n_trials[1]
  expect_equal(g$sem, g$sd / sqrt(n))                    # SEM = sd / sqrt(n)
  expect_equal(g$ci_upper - g$ci_lower, 2 * z * g$sd / sqrt(n))  # AC2 width ~ 1/sqrt(n)
  # more trials -> proportionally narrower band for the same sd
  ga2 <- plotErscrGrandAverage(make_eda_with_scr(n_events = 16), plot = FALSE)
  g2 <- ga2$grand_average[ga2$grand_average$type == "stimulus", ]
  expect_equal(g2$sem, g2$sd / sqrt(g2$n_trials[1]))
  expect_equal(g2$n_trials[1], 16)
})

test_that("plotErscrGrandAverage draws and returns an object without error", {
  x <- make_eda_with_scr()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(res <- plotErscrGrandAverage(x))         # AC3
  expect_s3_class(res, "erscr_grand_average")
  expect_output(print(res), "erscr_grand_average")
})

test_that("multiple event types are averaged and overlaid", {
  x <- make_eda_with_scr(n_time = 4000, sr = 10, n_events = 6)
  onsets <- seq(40, 360, length.out = 6)
  x <- PhysioCore::setEvents(x, PhysioCore::PhysioEvents(
    onset = onsets, duration = rep(0.5, 6),
    type = rep(c("A", "B"), 3), value = paste0("E", seq_len(6))))
  ga <- plotErscrGrandAverage(x, plot = FALSE)
  expect_setequal(unique(ga$grand_average$type), c("A", "B"))
  expect_equal(ga$grand_average$n_trials[ga$grand_average$type == "A"][1], 3)
  # selecting a single type
  gb <- plotErscrGrandAverage(x, event_types = "B", plot = FALSE)
  expect_equal(unique(gb$grand_average$type), "B")
})

test_that("plotErscrGrandAverage validates its inputs", {
  x <- make_eda_with_scr()
  expect_error(plotErscrGrandAverage(x, window = c(6, -1)), "pre < post")
  expect_error(plotErscrGrandAverage(x, baseline = c(-5, 0)), "within")
  expect_error(plotErscrGrandAverage(x, channel = 99), "channel")
  expect_error(plotErscrGrandAverage(x, event_types = "nope"), "present")
  expect_error(plotErscrGrandAverage(list()), "PhysioExperiment")
  # no events attached
  flat <- edaSimulate(n_time = 600, sr = 10, scr_count = 0, seed = 1)
  expect_error(plotErscrGrandAverage(flat), "no events")
})

test_that("epochs extending beyond the recording are dropped, keeping mean/CI consistent", {
  # an event so near the end that its [-1, 6] s window runs off the recording
  x <- make_eda_with_scr(n_time = 700, sr = 10, n_events = 3)
  onsets <- c(20, 40, 69)                                # 69 s + 6 s > 70 s record
  x <- PhysioCore::setEvents(x, PhysioCore::PhysioEvents(
    onset = onsets, duration = rep(0.5, 3),
    type = rep("stimulus", 3), value = paste0("E", seq_len(3))))
  expect_message(ga <- plotErscrGrandAverage(x, plot = FALSE), "dropped 1")
  g <- ga$grand_average
  expect_equal(g$n_trials[1], 2)                        # only the 2 complete epochs
  expect_false(anyNA(g$mean))                           # no ragged NA latencies
  expect_equal(g$mean, rowMeans(ga$epochs[["stimulus"]]), tolerance = 1e-6)
})

test_that("an empty baseline sub-window errors instead of silently skipping correction", {
  x <- make_eda_with_scr(n_time = 400, sr = 4, n_events = 4)
  # a narrow baseline that contains no sample on the 4 Hz grid
  expect_error(plotErscrGrandAverage(x, baseline = c(-0.15, -0.05), plot = FALSE),
               "no samples")
  # baseline = NULL is the explicit no-correction path (no error)
  expect_s3_class(plotErscrGrandAverage(x, baseline = NULL, plot = FALSE),
                  "erscr_grand_average")
})

test_that("plot label overrides via ... do not collide with fixed arguments", {
  x <- make_eda_with_scr()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plotErscrGrandAverage(x, main = "Custom", xlab = "t", ylab = "y"))
  expect_s3_class(plotErscrGrandAverage(x, main = "Custom"), "erscr_grand_average")
})

test_that("a non-integer channel is rejected rather than silently truncated", {
  x <- make_eda_with_scr(n_channels = 2)
  expect_error(plotErscrGrandAverage(x, channel = 1.5, plot = FALSE), "integer")
})
