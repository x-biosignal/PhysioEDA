#' Detect Skin Conductance Response (SCR) Peaks
#'
#' Identifies SCR peaks in an EDA signal using either a gradient-based
#' zero-crossing method or an amplitude threshold method. Returns onset,
#' peak, amplitude, rise time, and recovery time for each detected SCR.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param method Detection method: \code{"gradient"} (first-derivative
#'   zero-crossing) or \code{"threshold"} (amplitude threshold). Default
#'   is \code{"gradient"}.
#' @param amplitude_min Minimum SCR amplitude in microsiemens (default: 0.01).
#' @param rise_time_min Minimum rise time in seconds (default: 0.1).
#' @param rise_time_max Maximum rise time in seconds (default: 5.0).
#' @param assay_name Name of the input assay. If NULL, uses \code{"phasic"}
#'   if available, otherwise \code{defaultAssay(x)}.
#' @return A \code{data.frame} with one row per detected SCR and the following
#'   columns:
#'   \describe{
#'     \item{channel}{Character channel label.}
#'     \item{onset_sample}{Integer sample index of SCR onset.}
#'     \item{onset_sec}{Numeric onset time in seconds.}
#'     \item{peak_sample}{Integer sample index of SCR peak.}
#'     \item{peak_sec}{Numeric peak time in seconds.}
#'     \item{amplitude}{Numeric SCR amplitude in microsiemens (peak minus onset).}
#'     \item{rise_time}{Numeric rise time from onset to peak in seconds.}
#'     \item{recovery_time}{Numeric 50\% recovery time in seconds, or \code{NA}.}
#'   }
#'   Returns an empty \code{data.frame} with the same columns if no peaks are
#'   found.
#' @references
#' Bach, D.R., et al. (2010). "Modelling event-related skin conductance
#' responses." \emph{International Journal of Psychophysiology}, 75(3),
#' 349-356. \doi{10.1016/j.ijpsycho.2010.01.005}
#'
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#' @seealso \code{\link{edaDecompose}} for tonic/phasic decomposition (run first),
#'   \code{\link{edaFeatures}} for summary feature extraction,
#'   \code{\link{plotPeaks}} for peak visualization,
#'   \code{\link{edaErscr}} for event-related SCR analysis
#' @export
edaPeaks <- function(x, method = c("gradient", "threshold"),
                     amplitude_min = 0.01, rise_time_min = 0.1,
                     rise_time_max = 5.0, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  # Auto-resolve assay: prefer "phasic" if available

  if (is.null(assay_name)) {
    anames <- SummarizedExperiment::assayNames(x)
    if ("phasic" %in% anames) {
      assay_name <- "phasic"
    } else {
      assay_name <- defaultAssay(x)
    }
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_channels <- ncol(data)

  # Column template for empty results
  empty_df <- data.frame(
    channel = character(0),
    onset_sample = integer(0),
    onset_sec = numeric(0),
    peak_sample = integer(0),
    peak_sec = numeric(0),
    amplitude = numeric(0),
    rise_time = numeric(0),
    recovery_time = numeric(0),
    stringsAsFactors = FALSE
  )

  # Get channel labels
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% colnames(cd)) {
    ch_labels <- as.character(cd$label)
  } else {
    ch_labels <- paste0("Ch", seq_len(n_channels))
  }

  results <- vector("list", n_channels)

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    if (method == "gradient") {
      peaks_df <- .detect_peaks_gradient(sig, sr, amplitude_min,
                                         rise_time_min, rise_time_max)
    } else {
      peaks_df <- .detect_peaks_threshold(sig, sr, amplitude_min,
                                          rise_time_min, rise_time_max)
    }

    if (nrow(peaks_df) > 0) {
      peaks_df$channel <- ch_labels[ch]
      results[[ch]] <- peaks_df
    }
  }

  result <- do.call(rbind, results)

  if (is.null(result) || nrow(result) == 0) {
    return(empty_df)
  }

  # Reorder columns to put channel first
  result <- result[, c("channel", "onset_sample", "onset_sec",
                        "peak_sample", "peak_sec", "amplitude",
                        "rise_time", "recovery_time")]
  rownames(result) <- NULL
  result
}

#' Detect SCR peaks using first-derivative zero-crossing
#' @keywords internal
.detect_peaks_gradient <- function(sig, sr, amplitude_min,
                                   rise_time_min, rise_time_max) {
  n <- length(sig)
  if (n < 3) {
    return(data.frame(
      onset_sample = integer(0), onset_sec = numeric(0),
      peak_sample = integer(0), peak_sec = numeric(0),
      amplitude = numeric(0), rise_time = numeric(0),
      recovery_time = numeric(0), stringsAsFactors = FALSE
    ))
  }

  dsig <- diff(sig)

  # Find peaks: derivative goes from positive to non-positive
  peak_indices <- integer(0)
  for (i in seq_len(length(dsig) - 1)) {
    if (dsig[i] > 0 && dsig[i + 1] <= 0) {
      # Peak is at index i+1 in the original signal
      peak_indices <- c(peak_indices, i + 1L)
    }
  }

  if (length(peak_indices) == 0) {
    return(data.frame(
      onset_sample = integer(0), onset_sec = numeric(0),
      peak_sample = integer(0), peak_sec = numeric(0),
      amplitude = numeric(0), rise_time = numeric(0),
      recovery_time = numeric(0), stringsAsFactors = FALSE
    ))
  }

  onset_samples <- integer(0)
  peak_samples <- integer(0)
  amplitudes <- numeric(0)
  rise_times <- numeric(0)
  recovery_times <- numeric(0)

  for (pk in peak_indices) {
    # Walk backward from peak to find onset (where signal stops decreasing)
    onset <- pk
    for (j in seq(pk - 1, 1, by = -1)) {
      if (sig[j] <= sig[j + 1]) {
        onset <- j
      } else {
        break
      }
    }
    if (onset == pk) next

    amp <- sig[pk] - sig[onset]
    rt <- (pk - onset) / sr

    # Filter by criteria
    if (amp < amplitude_min) next
    if (rt < rise_time_min || rt > rise_time_max) next

    # Recovery time: time for signal to drop to 50% of amplitude after peak
    half_amp <- sig[pk] - amp * 0.5
    rec_time <- NA_real_
    for (j in seq(pk + 1, n)) {
      if (sig[j] <= half_amp) {
        rec_time <- (j - pk) / sr
        break
      }
    }

    onset_samples <- c(onset_samples, onset)
    peak_samples <- c(peak_samples, pk)
    amplitudes <- c(amplitudes, amp)
    rise_times <- c(rise_times, rt)
    recovery_times <- c(recovery_times, rec_time)
  }

  data.frame(
    onset_sample = onset_samples,
    onset_sec = (onset_samples - 1) / sr,
    peak_sample = peak_samples,
    peak_sec = (peak_samples - 1) / sr,
    amplitude = amplitudes,
    rise_time = rise_times,
    recovery_time = recovery_times,
    stringsAsFactors = FALSE
  )
}

#' Detect SCR peaks using amplitude threshold
#' @keywords internal
.detect_peaks_threshold <- function(sig, sr, amplitude_min,
                                    rise_time_min, rise_time_max) {
  n <- length(sig)
  if (n < 3) {
    return(data.frame(
      onset_sample = integer(0), onset_sec = numeric(0),
      peak_sample = integer(0), peak_sec = numeric(0),
      amplitude = numeric(0), rise_time = numeric(0),
      recovery_time = numeric(0), stringsAsFactors = FALSE
    ))
  }

  # Find segments where signal exceeds amplitude_min
  above <- sig > amplitude_min

  onset_samples <- integer(0)
  peak_samples <- integer(0)
  amplitudes <- numeric(0)
  rise_times <- numeric(0)
  recovery_times <- numeric(0)

  in_segment <- FALSE
  seg_start <- 1L

  for (i in seq_len(n)) {
    if (above[i] && !in_segment) {
      in_segment <- TRUE
      seg_start <- i
    } else if (!above[i] && in_segment) {
      in_segment <- FALSE
      seg_end <- i - 1L

      # Find peak within segment
      seg_range <- seg_start:seg_end
      pk_local <- which.max(sig[seg_range])
      pk <- seg_start + pk_local - 1L

      # Onset is the start of the segment
      onset <- seg_start

      amp <- sig[pk] - sig[onset]
      rt <- (pk - onset) / sr

      # Filter by criteria
      if (amp < amplitude_min) next
      if (rt < rise_time_min || rt > rise_time_max) next

      # Recovery time: time for signal to drop to 50% of amplitude after peak
      half_amp <- sig[pk] - amp * 0.5
      rec_time <- NA_real_
      for (j in seq(pk + 1, n)) {
        if (sig[j] <= half_amp) {
          rec_time <- (j - pk) / sr
          break
        }
      }

      onset_samples <- c(onset_samples, onset)
      peak_samples <- c(peak_samples, pk)
      amplitudes <- c(amplitudes, amp)
      rise_times <- c(rise_times, rt)
      recovery_times <- c(recovery_times, rec_time)
    }
  }

  # Handle segment that extends to end of signal
  if (in_segment) {
    seg_end <- n
    seg_range <- seg_start:seg_end
    pk_local <- which.max(sig[seg_range])
    pk <- seg_start + pk_local - 1L
    onset <- seg_start

    amp <- sig[pk] - sig[onset]
    rt <- (pk - onset) / sr

    if (amp >= amplitude_min && rt >= rise_time_min && rt <= rise_time_max) {
      half_amp <- sig[pk] - amp * 0.5
      rec_time <- NA_real_
      for (j in seq(pk + 1, n)) {
        if (sig[j] <= half_amp) {
          rec_time <- (j - pk) / sr
          break
        }
      }

      onset_samples <- c(onset_samples, onset)
      peak_samples <- c(peak_samples, pk)
      amplitudes <- c(amplitudes, amp)
      rise_times <- c(rise_times, rt)
      recovery_times <- c(recovery_times, rec_time)
    }
  }

  data.frame(
    onset_sample = onset_samples,
    onset_sec = (onset_samples - 1) / sr,
    peak_sample = peak_samples,
    peak_sec = (peak_samples - 1) / sr,
    amplitude = amplitudes,
    rise_time = rise_times,
    recovery_time = recovery_times,
    stringsAsFactors = FALSE
  )
}
