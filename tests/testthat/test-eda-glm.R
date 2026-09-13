library(testthat)
library(PhysioEDA)

test_that("edaGLM recovers per-condition amplitude ratios within 10%", {
  # acceptance: GLM betas recover the known amplitude ratios within 10%
  for (seed in 1:4) {
    fx <- make_eda_glm_fixture(seed)
    g <- edaGLM(fx$pe, fx$events)
    b <- g$coefficients$beta
    true <- fx$conds$amp
    rb <- b / b[1]
    rt <- true / true[1]
    expect_lt(max(abs(rb - rt) / rt), 0.10)
  }
})

test_that("edaGLM returns a full coefficient table and diagnostics", {
  fx <- make_eda_glm_fixture(1)
  g <- edaGLM(fx$pe, fx$events)
  expect_s3_class(g, "eda_glm")
  expect_s3_class(g$coefficients, "data.frame")
  expect_setequal(names(g$coefficients), c("condition", "beta", "se", "t", "p"))
  expect_equal(g$coefficients$condition, fx$conds$condition)
  # standard errors and t statistics are finite and positive SE
  expect_true(all(is.finite(g$coefficients$se) & g$coefficients$se > 0))
  expect_true(all(is.finite(g$coefficients$t)))
  # design is full rank on well-separated events
  expect_false(g$collinear)
  expect_equal(g$rank, g$n_regressors)
  expect_gt(g$r_squared, 0.95)
})

test_that("edaGLM warns and flags collinear (rank-deficient) designs", {
  fx <- make_eda_glm_fixture(1)
  # two conditions sharing identical onsets are perfectly collinear
  ev <- data.frame(onset = c(20, 60, 100, 20, 60, 100),
                   condition = rep(c("X", "Y"), each = 3))
  expect_warning(g <- edaGLM(fx$pe, ev), "rank-deficient")
  expect_true(g$collinear)
  expect_lt(g$rank, g$n_regressors)
})

test_that("edaGLM derivative basis adds per-condition derivative regressors", {
  fx <- make_eda_glm_fixture(2)
  g <- edaGLM(fx$pe, fx$events, basis = "derivative")
  expect_false(is.null(g$derivative))
  expect_equal(nrow(g$derivative), nrow(fx$conds))
  # canonical amplitude ratios still recovered
  b <- g$coefficients$beta
  rb <- b / b[1]; rt <- fx$conds$amp / fx$conds$amp[1]
  expect_lt(max(abs(rb - rt) / rt), 0.12)
})

test_that("edaGLM validates its inputs", {
  fx <- make_eda_glm_fixture(1)
  expect_error(edaGLM(fx$pe, data.frame(x = 1)), "onset")
  expect_error(edaGLM(fx$pe, fx$events, channel = 99), "out of range")
})
