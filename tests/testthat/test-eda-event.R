library(testthat)
library(PhysioEDA)

test_that("edaErscr returns data.frame with correct columns", {
  set.seed(42)
  x <- make_eda_event()
  x <- edaDecompose(x)
  result <- edaErscr(x)

  expect_s3_class(result, "data.frame")
  expected_cols <- c("event_index", "event_onset", "channel",
                     "scr_present", "scr_amplitude", "scr_latency",
                     "scr_rise_time", "scr_recovery_time")
  expect_equal(colnames(result), expected_cols)
})

test_that("edaErscr has one row per event per channel", {
  set.seed(42)
  x <- make_eda_event(n_channels = 2, n_events = 4)
  x <- edaDecompose(x)
  result <- edaErscr(x)

  expect_equal(nrow(result), 8)
})

test_that("edaErscr scr_present is logical", {
  set.seed(42)
  x <- make_eda_event()
  x <- edaDecompose(x)
  result <- edaErscr(x)

  expect_type(result$scr_present, "logical")
})

test_that("edaErscr errors when no events found", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x)

  expect_error(edaErscr(x), "No events")
})

test_that("edaErscr uses all events when event_type is NULL", {
  set.seed(42)
  x <- make_eda_event(n_channels = 1, n_events = 4)
  x <- edaDecompose(x)
  result <- edaErscr(x)

  expect_equal(nrow(result), 4)
})
