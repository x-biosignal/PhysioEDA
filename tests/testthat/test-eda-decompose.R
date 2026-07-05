library(testthat)
library(PhysioEDA)

test_that("edaDecompose with median method adds tonic and phasic assays", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("tonic" %in% nms)
  expect_true("phasic" %in% nms)
})

test_that("edaDecompose preserves dimensions", {
  set.seed(42)
  x <- make_eda(n_channels = 2)
  x <- edaDecompose(x, method = "median")

  raw_dim <- dim(SummarizedExperiment::assay(x, "raw"))
  tonic_dim <- dim(SummarizedExperiment::assay(x, "tonic"))
  phasic_dim <- dim(SummarizedExperiment::assay(x, "phasic"))

  expect_equal(tonic_dim, raw_dim)
  expect_equal(phasic_dim, raw_dim)
})

test_that("edaDecompose median produces smoother tonic than raw", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  raw <- SummarizedExperiment::assay(x, "raw")[, 1]
  tonic <- SummarizedExperiment::assay(x, "tonic")[, 1]

  expect_lt(sd(tonic), sd(raw))
})

test_that("edaDecompose tonic + phasic approximately equals raw", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  raw <- SummarizedExperiment::assay(x, "raw")
  tonic <- SummarizedExperiment::assay(x, "tonic")
  phasic <- SummarizedExperiment::assay(x, "phasic")

  expect_equal(tonic + phasic, raw, tolerance = 1e-10)
})

test_that("edaDecompose highpass method works", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "highpass")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("tonic" %in% nms)
  expect_true("phasic" %in% nms)
})

test_that("edaDecompose stores method in metadata", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median")

  meta <- S4Vectors::metadata(x)
  expect_equal(meta$eda_decompose$method, "median")
})

test_that("edaDecompose supports custom output assay names", {
  set.seed(42)
  x <- make_eda()
  x <- edaDecompose(x, method = "median",
                    output_tonic = "my_tonic",
                    output_phasic = "my_phasic")

  nms <- SummarizedExperiment::assayNames(x)
  expect_true("my_tonic" %in% nms)
  expect_true("my_phasic" %in% nms)
})
