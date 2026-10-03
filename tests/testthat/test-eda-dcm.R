library(testthat)
library(PhysioEDA)

test_that("edaDCM per-event amplitudes correlate r > 0.9 with the ground truth", {
  # acceptance: DCM amplitude estimates correlate r > 0.9 with truth
  for (seed in 1:4) {
    fx <- make_eda_dcm_fixture(seed)
    d <- edaDCM(fx$pe, fx$events, refine_latency = FALSE)
    expect_gt(cor(d$estimates$amplitude, fx$amps), 0.9)
  }
})

test_that("edaDCM returns non-negative per-event estimates and a good fit", {
  fx <- make_eda_dcm_fixture(1)
  d <- edaDCM(fx$pe, fx$events)
  expect_s3_class(d, "eda_dcm")
  expect_s3_class(d$estimates, "data.frame")
  expect_setequal(names(d$estimates),
                  c("event", "onset", "onset_refined", "amplitude"))
  expect_equal(nrow(d$estimates), nrow(fx$events))
  expect_true(all(d$estimates$amplitude >= 0))          # latent input is non-negative
  expect_gt(d$r_squared, 0.9)
})

test_that("edaDCM latency refinement recovers shifted onsets", {
  # shift the true onsets; the recorded onsets are the shifted ones, and
  # refinement should move onset_refined back towards the true generators
  fx <- make_eda_dcm_fixture(3)
  true_onsets <- fx$events$onset
  shifted <- data.frame(onset = true_onsets + 0.8)      # 0.8 s late
  d <- edaDCM(fx$pe, shifted, refine_latency = TRUE, latency_window = 1.5)
  err_before <- mean(abs(shifted$onset - true_onsets))
  err_after <- mean(abs(d$estimates$onset_refined - true_onsets))
  expect_lt(err_after, err_before)
  expect_gt(cor(d$estimates$amplitude, fx$amps), 0.9)
})

test_that("edaDCM carries a condition column through when present", {
  fx <- make_eda_dcm_fixture(1)
  fx$events$condition <- rep(c("A", "B"), length.out = nrow(fx$events))
  d <- edaDCM(fx$pe, fx$events, refine_latency = FALSE)
  expect_true("condition" %in% names(d$estimates))
})

test_that("edaDCM validates its inputs", {
  fx <- make_eda_dcm_fixture(1)
  expect_error(edaDCM(fx$pe, data.frame(x = 1)), "onset")
})
