#' Extract EDA Features per Channel
#'
#' Computes summary features from electrodermal activity data, including
#' skin conductance response (SCR) statistics and skin conductance level
#' (SCL) measures. Requires tonic/phasic decomposition to have been
#' performed (via \code{\link{edaDecompose}}).
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param peaks Optional pre-computed peaks data.frame from
#'   \code{\link{edaPeaks}}. If NULL, peaks are computed automatically
#'   using the gradient method.
#' @param window Optional numeric vector of length 2 giving the time window
#'   in seconds as \code{c(start_sec, end_sec)} to restrict analysis.
#' @param assay_name Name of the input assay for peak detection. If NULL,
#'   uses \code{defaultAssay(x)}.
#' @return A \code{data.frame} with one row per channel and the following
#'   columns:
#'   \describe{
#'     \item{channel}{Character channel label.}
#'     \item{scr_count}{Integer number of detected SCR peaks.}
#'     \item{scr_rate_per_min}{Numeric SCR count per minute.}
#'     \item{mean_amplitude}{Numeric mean SCR peak amplitude in microsiemens
#'       (0 if no peaks detected).}
#'     \item{mean_scl}{Numeric mean tonic skin conductance level.}
#'     \item{scl_sd}{Numeric standard deviation of the tonic signal.}
#'     \item{auc_phasic}{Numeric area under the phasic curve (positive values
#'       only, in uS*s), or \code{NA} if no phasic assay exists.}
#'     \item{ns_scr_freq}{Numeric non-specific SCR frequency (per minute).}
#'   }
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#'
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#' @seealso \code{\link{edaDecompose}} for tonic/phasic decomposition (run first),
#'   \code{\link{edaPeaks}} for SCR peak detection,
#'   \code{\link{edaErscr}} for event-related SCR analysis
#' @export
edaFeatures <- function(x, peaks = NULL, window = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  sr <- samplingRate(x)

  # Determine assay names for tonic and phasic
  assay_names <- SummarizedExperiment::assayNames(x)
  tonic_assay <- if ("tonic" %in% assay_names) "tonic" else defaultAssay(x)
  has_phasic <- "phasic" %in% assay_names

  tonic_data <- SummarizedExperiment::assay(x, tonic_assay)
  phasic_data <- if (has_phasic) SummarizedExperiment::assay(x, "phasic") else NULL

  n_time <- nrow(tonic_data)
  n_channels <- ncol(tonic_data)

  # Subset data to window if specified
  if (!is.null(window)) {
    stopifnot(is.numeric(window), length(window) == 2, window[1] < window[2])
    start_idx <- max(1L, as.integer(floor(window[1] * sr)) + 1L)
    end_idx <- min(n_time, as.integer(floor(window[2] * sr)) + 1L)
    tonic_data <- tonic_data[start_idx:end_idx, , drop = FALSE]
    if (has_phasic) {
      phasic_data <- phasic_data[start_idx:end_idx, , drop = FALSE]
    }
    duration_sec <- (end_idx - start_idx) / sr
  } else {
    duration_sec <- (n_time - 1) / sr
  }

  duration_min <- duration_sec / 60

  # Auto-compute peaks if not provided
  if (is.null(peaks)) {
    peaks <- edaPeaks(x, method = "gradient", assay_name = assay_name)
  }

  # Filter peaks to window if specified
  if (!is.null(window) && nrow(peaks) > 0) {
    peaks <- peaks[peaks$onset_sec >= window[1] & peaks$onset_sec <= window[2], ]
  }

  # Get channel labels
  ch_labels <- SummarizedExperiment::colData(x)$label
  if (is.null(ch_labels)) {
    ch_labels <- paste0("Ch", seq_len(n_channels))
  }

  # Compute features per channel
  results <- vector("list", n_channels)
  for (ch in seq_len(n_channels)) {
    ch_label <- ch_labels[ch]
    ch_peaks <- peaks[peaks$channel == ch_label, ]

    scr_count <- nrow(ch_peaks)
    scr_rate <- if (duration_min > 0) scr_count / duration_min else 0
    mean_amp <- if (scr_count > 0) mean(ch_peaks$amplitude) else 0

    tonic_ch <- tonic_data[, ch]
    mean_scl <- mean(tonic_ch)
    scl_sd_val <- stats::sd(tonic_ch)

    auc <- NA_real_
    if (has_phasic) {
      phasic_ch <- phasic_data[, ch]
      positive_vals <- phasic_ch[phasic_ch > 0]
      auc <- sum(positive_vals) / sr
    }

    results[[ch]] <- data.frame(
      channel = ch_label,
      scr_count = scr_count,
      scr_rate_per_min = scr_rate,
      mean_amplitude = mean_amp,
      mean_scl = mean_scl,
      scl_sd = scl_sd_val,
      auc_phasic = auc,
      ns_scr_freq = scr_rate,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}
