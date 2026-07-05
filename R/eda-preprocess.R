#' Filter EDA Signal
#'
#' Applies a frequency-domain (FFT-based) filter to an EDA signal. Supports
#' lowpass, highpass, and bandpass filter types using a smooth quadratic
#' transition in the frequency domain.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param type Filter type: \code{"lowpass"}, \code{"highpass"}, or
#'   \code{"bandpass"}.
#' @param cutoff Cutoff frequency in Hz. A single numeric value for lowpass or
#'   highpass, or a numeric vector of length 2 (\code{c(low, high)}) for
#'   bandpass.
#' @param order Filter steepness parameter controlling the smoothness of the
#'   frequency transition (default: 2). Higher values produce sharper rolloff.
#' @param assay_name Name of the input assay. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_assay Name for the output assay (default: "filtered").
#' @return A modified \code{\link[PhysioCore]{PhysioExperiment}} with filtered
#'   data stored in the \code{output_assay} and filter parameters recorded in
#'   \code{metadata(x)$eda_filter} (a list with \code{type}, \code{cutoff},
#'   \code{order}, \code{assay_name}, and \code{output_assay}).
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{edaDownsample}} for decimation with anti-aliasing,
#'   \code{\link{edaArtifact}} for artifact detection and correction,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition
#' @export
edaFilter <- function(x, type = c("lowpass", "highpass", "bandpass"),
                      cutoff, order = 2,
                      assay_name = NULL, output_assay = "filtered") {
  stopifnot(inherits(x, "PhysioExperiment"))
  type <- match.arg(type)

  # Validate cutoff
  if (type == "bandpass") {
    stopifnot(is.numeric(cutoff), length(cutoff) == 2, cutoff[1] < cutoff[2])
  } else {
    stopifnot(is.numeric(cutoff), length(cutoff) == 1, cutoff > 0)
  }
  stopifnot(is.numeric(order), length(order) == 1, order > 0)

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)
  nyquist <- sr / 2

  # Validate cutoff against Nyquist
  max_cutoff <- if (type == "bandpass") max(cutoff) else cutoff
  if (max_cutoff >= nyquist) {
    stop(sprintf(
      "Cutoff frequency (%.2f Hz) must be below Nyquist frequency (%.2f Hz).",
      max_cutoff, nyquist
    ), call. = FALSE)
  }

  # Build frequency vector for all FFT bins
  freq_bins <- (seq_len(n_time) - 1) * sr / n_time
  # Mirror frequencies above Nyquist
  freq_bins[freq_bins > nyquist] <- sr - freq_bins[freq_bins > nyquist]

  # Build frequency-domain mask
  mask <- .build_filter_mask(freq_bins, type, cutoff, order)

  # Apply filter channel by channel
  filtered_mat <- matrix(NA_real_, nrow = n_time, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    ft <- stats::fft(data[, ch])
    filtered_ft <- ft * mask
    filtered_mat[, ch] <- Re(stats::fft(filtered_ft, inverse = TRUE)) / n_time
  }

  # Store results as new assay
  current_assays <- SummarizedExperiment::assays(x)
  current_assays[[output_assay]] <- filtered_mat
  SummarizedExperiment::assays(x) <- current_assays

  # Store filter parameters in metadata
  meta <- S4Vectors::metadata(x)
  meta$eda_filter <- list(
    type = type,
    cutoff = cutoff,
    order = order,
    assay_name = assay_name,
    output_assay = output_assay
  )
  S4Vectors::metadata(x) <- meta

  x
}

#' Build FFT filter mask
#'
#' Constructs a smooth frequency-domain mask for lowpass, highpass, or bandpass
#' filtering using a power-law (raised cosine-like) transition.
#'
#' @param freq_bins Numeric vector of frequencies for each FFT bin.
#' @param type Filter type.
#' @param cutoff Cutoff frequency or frequencies.
#' @param order Steepness parameter.
#' @return Numeric vector of mask values (0 to 1) for each FFT bin.
#' @keywords internal
.build_filter_mask <- function(freq_bins, type, cutoff, order) {
  n <- length(freq_bins)
  mask <- numeric(n)

  if (type == "lowpass") {
    for (i in seq_len(n)) {
      f <- freq_bins[i]
      if (f <= cutoff) {
        mask[i] <- 1.0
      } else {
        mask[i] <- 1.0 / (1.0 + (f / cutoff)^(2 * order))
      }
    }
  } else if (type == "highpass") {
    for (i in seq_len(n)) {
      f <- freq_bins[i]
      if (f >= cutoff) {
        mask[i] <- 1.0
      } else if (f == 0) {
        mask[i] <- 0.0
      } else {
        mask[i] <- 1.0 / (1.0 + (cutoff / f)^(2 * order))
      }
    }
  } else {
    # bandpass: product of highpass at cutoff[1] and lowpass at cutoff[2]
    low <- cutoff[1]
    high <- cutoff[2]
    for (i in seq_len(n)) {
      f <- freq_bins[i]
      if (f == 0) {
        mask[i] <- 0.0
      } else {
        hp <- 1.0 / (1.0 + (low / f)^(2 * order))
        lp <- 1.0 / (1.0 + (f / high)^(2 * order))
        mask[i] <- hp * lp
      }
    }
  }

  mask
}

#' Downsample EDA Signal
#'
#' Reduces the sampling rate of an EDA signal by first applying an anti-aliasing
#' lowpass filter at half the target sampling rate, then decimating the signal.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param target_sr Target sampling rate in Hz. Must be lower than the current
#'   sampling rate.
#' @param assay_name Name of the input assay. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @return A new \code{\link[PhysioCore]{PhysioExperiment}} with the downsampled
#'   signal in the \code{"raw"} assay and the sampling rate set to the actual
#'   achieved rate. Events are preserved. Downsampling parameters are stored in
#'   \code{metadata(x)$eda_downsample} (a list with \code{original_sr},
#'   \code{target_sr}, \code{factor}, \code{original_n_time}, and
#'   \code{new_n_time}).
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{edaFilter}} for frequency-domain filtering,
#'   \code{\link{edaArtifact}} for artifact detection and correction
#' @export
edaDownsample <- function(x, target_sr, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.numeric(target_sr), length(target_sr) == 1, target_sr > 0)

  sr <- samplingRate(x)

  if (target_sr >= sr) {
    stop(sprintf(
      "target_sr (%.1f Hz) must be lower than current sampling rate (%.1f Hz).",
      target_sr, sr
    ), call. = FALSE)
  }

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Compute decimation factor
  factor <- round(sr / target_sr)
  actual_sr <- sr / factor

  # Anti-aliasing lowpass filter at target_sr / 2
  aa_cutoff <- actual_sr / 2
  nyquist <- sr / 2

  # Build frequency vector and lowpass mask
  freq_bins <- (seq_len(n_time) - 1) * sr / n_time
  freq_bins[freq_bins > nyquist] <- sr - freq_bins[freq_bins > nyquist]
  mask <- .build_filter_mask(freq_bins, "lowpass", aa_cutoff, order = 4)

  # Apply anti-aliasing filter and decimate
  keep_indices <- seq(1, n_time, by = factor)
  n_new <- length(keep_indices)

  downsampled_mat <- matrix(NA_real_, nrow = n_new, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    ft <- stats::fft(data[, ch])
    filtered <- Re(stats::fft(ft * mask, inverse = TRUE)) / n_time
    downsampled_mat[, ch] <- filtered[keep_indices]
  }

  # Get channel metadata
  cd <- SummarizedExperiment::colData(x)

  # Create new PhysioExperiment with downsampled data
  result <- PhysioExperiment(
    assays = list(raw = downsampled_mat),
    colData = cd,
    samplingRate = actual_sr
  )

  # Preserve events if they exist
  meta <- S4Vectors::metadata(x)
  if (!is.null(meta$events)) {
    result <- setEvents(result, getEvents(x))
  }

  # Store downsampling info in metadata
  result_meta <- S4Vectors::metadata(result)
  result_meta$eda_downsample <- list(
    original_sr = sr,
    target_sr = actual_sr,
    factor = factor,
    original_n_time = n_time,
    new_n_time = n_new
  )
  S4Vectors::metadata(result) <- result_meta

  result
}
