library(testthat)
library(PhysioEDA)

test_that("edaTransform log produces expected output", {
  set.seed(42)
  x <- make_eda()
  x <- edaTransform(x, method = "log")

  raw <- SummarizedExperiment::assay(x, "raw")[, 1]
  transformed <- SummarizedExperiment::assay(x, "log_transformed")[, 1]

  expect_equal(transformed, log1p(raw))
  expect_true(all(is.finite(transformed)))
})

test_that("edaTransform sqrt clamps negatives and transforms", {
  set.seed(42)
  x <- make_eda()
  x <- edaTransform(x, method = "sqrt")

  raw <- SummarizedExperiment::assay(x, "raw")[, 1]
  transformed <- SummarizedExperiment::assay(x, "sqrt_transformed")[, 1]

  expect_equal(transformed, sqrt(pmax(raw, 0)))
  expect_true(all(transformed >= 0))
})

test_that("edaTransform zscore has mean ~0 and sd ~1 per channel", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  x <- edaTransform(x, method = "zscore")

  transformed <- SummarizedExperiment::assay(x, "zscore_transformed")

  for (ch in seq_len(ncol(transformed))) {
    expect_equal(mean(transformed[, ch]), 0, tolerance = 1e-10)
    expect_equal(sd(transformed[, ch]), 1, tolerance = 1e-10)
  }
})

test_that("edaTransform range produces values in [0,1]", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  x <- edaTransform(x, method = "range")

  transformed <- SummarizedExperiment::assay(x, "range_transformed")

  expect_true(all(transformed >= 0))
  expect_true(all(transformed <= 1))
  for (ch in seq_len(ncol(transformed))) {
    expect_equal(min(transformed[, ch]), 0)
    expect_equal(max(transformed[, ch]), 1)
  }
})

test_that("edaTransform stores metadata correctly", {
  set.seed(42)
  x <- make_eda()
  x <- edaTransform(x, method = "zscore")

  meta <- S4Vectors::metadata(x)$eda_transform
  expect_equal(meta$method, "zscore")
  expect_true(!is.null(meta$means))
  expect_true(!is.null(meta$sds))
  expect_equal(length(meta$means), 1)
  expect_equal(length(meta$sds), 1)
})

test_that("edaTransform preserves dimensions", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  raw_dim <- dim(SummarizedExperiment::assay(x, "raw"))

  for (method in c("log", "sqrt", "zscore", "range")) {
    x2 <- edaTransform(x, method = method)
    out_name <- paste0(method, "_transformed")
    expect_equal(dim(SummarizedExperiment::assay(x2, out_name)), raw_dim,
                 info = paste("dimensions preserved for", method))
  }
})

test_that("edaUntransform recovers original for all methods", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  raw <- SummarizedExperiment::assay(x, "raw")

  for (method in c("log", "sqrt", "zscore", "range")) {
    x2 <- edaTransform(x, method = method)
    x2 <- edaUntransform(x2)
    recovered <- SummarizedExperiment::assay(x2, "untransformed")
    expect_equal(recovered, raw, tolerance = 1e-10,
                 info = paste("untransform recovers original for", method))
  }
})

test_that("edaUntransform errors when no transform metadata exists", {
  set.seed(42)
  x <- make_eda()
  expect_error(edaUntransform(x), "No transform metadata found")
})

test_that("edaTransform supports custom output assay name", {
  set.seed(42)
  x <- make_eda()
  x <- edaTransform(x, method = "log", output_assay = "my_log")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("my_log" %in% nms)
})

test_that("edaTransform rejects invalid method", {
  set.seed(42)
  x <- make_eda()
  expect_error(edaTransform(x, method = "invalid"))
})

test_that("edaTransform range metadata stores mins and maxs", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  x <- edaTransform(x, method = "range")

  meta <- S4Vectors::metadata(x)$eda_transform
  expect_equal(meta$method, "range")
  expect_equal(length(meta$mins), 2)
  expect_equal(length(meta$maxs), 2)
})
