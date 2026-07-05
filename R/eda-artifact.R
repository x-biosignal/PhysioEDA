#' Detect and Correct EDA Artifacts
#'
#' Identifies artifacts in electrodermal activity (EDA) signals using one or
#' more detection methods (threshold, gradient, flatline) and optionally
#' corrects them via interpolation or NA replacement.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param methods Character vector of detection methods to apply. Any
#'   combination of \code{"threshold"}, \code{"gradient"}, and
#'   \code{"flatline"}. Default is all three.
#' @param threshold_range Numeric vector of length 2 giving the acceptable
#'   range of EDA values in microsiemens (default: \code{c(0.001, 60)}).
#'   Values outside this range are flagged as artifacts.
#' @param gradient_max Maximum allowable absolute gradient in microsiemens
#'   per sample. If NULL (default), computed as \code{10 / samplingRate(x)}
#'   (equivalent to 10 uS/sec).
#' @param flatline_sec Minimum duration in seconds of a constant-value
#'   segment to be flagged as a flatline artifact (default: 5).
#' @param correct Correction strategy: \code{"interpolate"} (linear
#'   interpolation across artifact regions), \code{"na"} (replace with
#'   NA), or \code{"none"} (detection only, no correction). Default is
#'   \code{"interpolate"}.
#' @param assay_name Name of the input assay. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_assay Name for the corrected output assay (default:
#'   \code{"cleaned"}). Only used when \code{correct != "none"}.
#' @return A modified \code{\link[PhysioCore]{PhysioExperiment}} with artifact
#'   information stored in \code{metadata(x)$eda_artifacts}, a list containing:
#'   \describe{
#'     \item{mask}{Logical matrix (time x channels) where \code{TRUE} indicates
#'       an artifact sample.}
#'     \item{summary}{A \code{data.frame} with columns \code{channel},
#'       \code{method}, \code{n_artifacts}, and \code{pct} giving artifact
#'       counts per channel per detection method.}
#'   }
#'   If \code{correct != "none"}, the corrected signal is stored in the
#'   \code{output_assay}.
#' @references
#' Kleckner, I.R., et al. (2018). "Simple, transparent, and flexible automated
#' quality assessment procedures for ambulatory electrodermal activity data."
#' \emph{IEEE Transactions on Biomedical Engineering}, 65(7), 1460-1467.
#' \doi{10.1109/TBME.2017.2758643}
#'
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{edaQuality}} for signal quality assessment,
#'   \code{\link{edaFilter}} for frequency-domain filtering,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition
#' @export
edaArtifact <- function(x, methods = c("threshold", "gradient", "flatline"),
                        threshold_range = c(0.001, 60),
                        gradient_max = NULL,
                        flatline_sec = 5,
                        correct = "interpolate",
                        assay_name = NULL,
                        output_assay = "cleaned") {
  stopifnot(inherits(x, "PhysioExperiment"))
  methods <- match.arg(methods, c("threshold", "gradient", "flatline"),
                       several.ok = TRUE)
  stopifnot(is.numeric(threshold_range), length(threshold_range) == 2,
            threshold_range[1] < threshold_range[2])
  stopifnot(is.numeric(flatline_sec), flatline_sec > 0)
  correct <- match.arg(correct, c("interpolate", "na", "none"))

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  if (is.null(gradient_max)) {
    gradient_max <- 10 / sr
  }
  stopifnot(is.numeric(gradient_max), gradient_max > 0)

  # Get channel labels
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% colnames(cd)) {
    ch_labels <- as.character(cd$label)
  } else {
    ch_labels <- paste0("Ch", seq_len(n_channels))
  }

  # Combined artifact mask (TRUE = artifact)
  mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)

  # Per-method masks for summary
  method_masks <- list()

  # --- Threshold detection ---
  if ("threshold" %in% methods) {
    thresh_mask <- data < threshold_range[1] | data > threshold_range[2]
    method_masks[["threshold"]] <- thresh_mask
    mask <- mask | thresh_mask
  }

  # --- Gradient detection ---
  if ("gradient" %in% methods) {
    grad_mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)
    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      d <- abs(diff(sig))
      bad <- which(d > gradient_max)
      # Flag both the sample before and after each large gradient
      bad_idx <- unique(c(bad, bad + 1L))
      bad_idx <- bad_idx[bad_idx >= 1L & bad_idx <= n_time]
      grad_mask[bad_idx, ch] <- TRUE
    }
    method_masks[["gradient"]] <- grad_mask
    mask <- mask | grad_mask
  }

  # --- Flatline detection ---
  if ("flatline" %in% methods) {
    flat_mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)
    min_run <- as.integer(ceiling(flatline_sec * sr))
    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      d <- diff(sig)
      is_zero <- d == 0
      # Find runs of zero-diff
      # is_zero[i] corresponds to diff between sig[i] and sig[i+1]
      run_len <- 0L
      for (i in seq_along(is_zero)) {
        if (is_zero[i]) {
          run_len <- run_len + 1L
        } else {
          if (run_len >= min_run) {
            # Zero-diff run covers is_zero indices (i-run_len) to (i-1)
            # In original signal: samples (i-run_len) to i
            start_idx <- i - run_len
            end_idx <- i
            flat_mask[start_idx:end_idx, ch] <- TRUE
          }
          run_len <- 0L
        }
      }
      # Handle run that extends to end of signal
      if (run_len >= min_run) {
        n_diff <- length(is_zero)
        start_idx <- n_diff - run_len + 1L
        end_idx <- n_time
        flat_mask[start_idx:end_idx, ch] <- TRUE
      }
    }
    method_masks[["flatline"]] <- flat_mask
    mask <- mask | flat_mask
  }

  # --- Build summary ---
  summary_rows <- vector("list", n_channels * length(methods))
  idx <- 1L
  for (ch in seq_len(n_channels)) {
    for (m in methods) {
      count <- sum(method_masks[[m]][, ch])
      summary_rows[[idx]] <- data.frame(
        channel = ch_labels[ch],
        method = m,
        n_artifacts = count,
        pct = 100 * count / n_time,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  summary_df <- do.call(rbind, summary_rows)

  # Store artifact info in metadata
  meta <- S4Vectors::metadata(x)
  meta$eda_artifacts <- list(
    mask = mask,
    summary = summary_df
  )
  S4Vectors::metadata(x) <- meta

  # --- Correction ---
  if (correct != "none") {
    corrected <- data
    for (ch in seq_len(n_channels)) {
      ch_mask <- mask[, ch]
      if (!any(ch_mask)) next

      if (correct == "na") {
        corrected[ch_mask, ch] <- NA_real_
      } else {
        # Linear interpolation across artifact regions
        good_idx <- which(!ch_mask)
        if (length(good_idx) >= 2) {
          bad_idx <- which(ch_mask)
          corrected[bad_idx, ch] <- stats::approx(
            x = good_idx,
            y = data[good_idx, ch],
            xout = bad_idx,
            rule = 2  # extend to boundary
          )$y
        } else {
          # Not enough good samples to interpolate
          corrected[ch_mask, ch] <- NA_real_
        }
      }
    }

    current_assays <- SummarizedExperiment::assays(x)
    current_assays[[output_assay]] <- corrected
    SummarizedExperiment::assays(x) <- current_assays
  }

  x
}
