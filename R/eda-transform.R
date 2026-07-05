#' Transform EDA Data
#'
#' Applies a mathematical transformation to EDA signal data. Common
#' transformations include log and square root (to reduce skewness of SCR
#' amplitudes), z-score standardization, and range normalization. Transformation
#' parameters are stored in metadata so that \code{\link{edaUntransform}} can
#' reverse the operation.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param method Transformation method: \code{"log"}, \code{"sqrt"},
#'   \code{"zscore"}, or \code{"range"}.
#' @param assay_name Name of the input assay. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_assay Name for the output assay. If NULL, defaults to
#'   \code{"{method}_transformed"}.
#' @return A modified PhysioExperiment with the transformed assay added and
#'   transformation parameters stored in \code{metadata(x)$eda_transform}.
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{edaUntransform}} for reversing the transformation,
#'   \code{\link{edaFeatures}} for feature extraction,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition
#' @export
edaTransform <- function(x, method = c("log", "sqrt", "zscore", "range"),
                         assay_name = NULL, output_assay = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }
  if (is.null(output_assay)) {
    output_assay <- paste0(method, "_transformed")
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  n_channels <- ncol(data)

  params <- list(method = method, assay_name = assay_name,
                 output_assay = output_assay)

  if (method == "log") {
    result <- log1p(data)

  } else if (method == "sqrt") {
    result <- sqrt(pmax(data, 0))

  } else if (method == "zscore") {
    means <- numeric(n_channels)
    sds <- numeric(n_channels)
    result <- matrix(NA_real_, nrow = nrow(data), ncol = n_channels)
    for (ch in seq_len(n_channels)) {
      means[ch] <- mean(data[, ch])
      sds[ch] <- sd(data[, ch])
      if (sds[ch] == 0) {
        result[, ch] <- 0
      } else {
        result[, ch] <- (data[, ch] - means[ch]) / sds[ch]
      }
    }
    params$means <- means
    params$sds <- sds

  } else if (method == "range") {
    mins <- numeric(n_channels)
    maxs <- numeric(n_channels)
    result <- matrix(NA_real_, nrow = nrow(data), ncol = n_channels)
    for (ch in seq_len(n_channels)) {
      mins[ch] <- min(data[, ch])
      maxs[ch] <- max(data[, ch])
      rng <- maxs[ch] - mins[ch]
      if (rng == 0) {
        result[, ch] <- 0
      } else {
        result[, ch] <- (data[, ch] - mins[ch]) / rng
      }
    }
    params$mins <- mins
    params$maxs <- maxs
  }

  current_assays <- SummarizedExperiment::assays(x)
  current_assays[[output_assay]] <- result
  SummarizedExperiment::assays(x) <- current_assays

  meta <- S4Vectors::metadata(x)
  meta$eda_transform <- params
  S4Vectors::metadata(x) <- meta

  x
}

#' Reverse EDA Transformation
#'
#' Reverses a transformation previously applied by \code{\link{edaTransform}},
#' using parameters stored in \code{metadata(x)$eda_transform}.
#'
#' @param x A PhysioExperiment object that has been transformed via
#'   \code{edaTransform}.
#' @param assay_name Name of the transformed assay to reverse. If NULL, uses
#'   the output_assay recorded in the transform metadata.
#' @param output_assay Name for the untransformed output assay. If NULL,
#'   defaults to \code{"untransformed"}.
#' @return A modified PhysioExperiment with the untransformed assay added.
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{edaTransform}} for applying the initial transformation,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition
#' @export
edaUntransform <- function(x, assay_name = NULL, output_assay = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  meta <- S4Vectors::metadata(x)
  params <- meta$eda_transform
  if (is.null(params)) {
    stop("No transform metadata found. Run edaTransform() first.")
  }

  if (is.null(assay_name)) {
    assay_name <- params$output_assay
  }
  if (is.null(output_assay)) {
    output_assay <- "untransformed"
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  method <- params$method
  n_channels <- ncol(data)

  if (method == "log") {
    result <- expm1(data)

  } else if (method == "sqrt") {
    result <- data^2

  } else if (method == "zscore") {
    means <- params$means
    sds <- params$sds
    result <- matrix(NA_real_, nrow = nrow(data), ncol = n_channels)
    for (ch in seq_len(n_channels)) {
      result[, ch] <- data[, ch] * sds[ch] + means[ch]
    }

  } else if (method == "range") {
    mins <- params$mins
    maxs <- params$maxs
    result <- matrix(NA_real_, nrow = nrow(data), ncol = n_channels)
    for (ch in seq_len(n_channels)) {
      rng <- maxs[ch] - mins[ch]
      result[, ch] <- data[, ch] * rng + mins[ch]
    }
  }

  current_assays <- SummarizedExperiment::assays(x)
  current_assays[[output_assay]] <- result
  SummarizedExperiment::assays(x) <- current_assays

  x
}
