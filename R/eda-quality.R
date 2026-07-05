#' Assess EDA Signal Quality
#'
#' Computes per-channel signal quality metrics for electrodermal activity data,
#' including basic statistics, flatline detection, artifact estimation, and an
#' overall quality score.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param assay_name Name of the input assay. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @return A data.frame with one row per channel and the following columns:
#'   \describe{
#'     \item{channel}{Channel label}
#'     \item{mean_sc}{Mean skin conductance}
#'     \item{sd_sc}{Standard deviation of skin conductance}
#'     \item{min_sc}{Minimum skin conductance}
#'     \item{max_sc}{Maximum skin conductance}
#'     \item{pct_negative}{Percentage of samples <= 0}
#'     \item{pct_flatline}{Percentage of signal in flatline segments (runs of
#'       consecutive near-zero differences longer than 1 second)}
#'     \item{pct_artifact}{Percentage of gradient-based artifact samples}
#'     \item{snr_db}{Estimated signal-to-noise ratio in dB (capped at 60)}
#'     \item{quality_score}{Overall quality score from 0 to 100}
#'     \item{quality_label}{"good" (>= 70), "acceptable" (>= 40), or "poor"}
#'   }
#' @references
#' Kleckner, I.R., et al. (2018). "Simple, transparent, and flexible automated
#' quality assessment procedures for ambulatory electrodermal activity data."
#' \emph{IEEE Transactions on Biomedical Engineering}, 65(7), 1460-1467.
#' \doi{10.1109/TBME.2017.2758643}
#' @seealso \code{\link{edaArtifact}} for artifact detection and correction,
#'   \code{\link{edaFilter}} for frequency-domain filtering,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition
#' @export
edaQuality <- function(x, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Get channel labels
  ch_labels <- SummarizedExperiment::colData(x)$label
  if (is.null(ch_labels)) {
    ch_labels <- paste0("Ch", seq_len(n_channels))
  }

  results <- vector("list", n_channels)
  for (ch in seq_len(n_channels)) {
    signal <- data[, ch]
    n <- length(signal)

    # Basic statistics
    mean_sc <- mean(signal)
    sd_sc <- stats::sd(signal)
    min_sc <- min(signal)
    max_sc <- max(signal)

    # Percentage of negative/zero samples
    pct_negative <- 100 * sum(signal <= 0) / n

    # Flatline detection: runs of abs(diff) < 1e-6 for > 1*sr consecutive
    d <- abs(diff(signal))
    flat_threshold <- 1e-6
    min_run <- as.integer(ceiling(sr))  # 1 second worth of samples
    flat_samples <- 0L
    run_len <- 0L
    for (i in seq_along(d)) {
      if (d[i] < flat_threshold) {
        run_len <- run_len + 1L
      } else {
        if (run_len > min_run) {
          flat_samples <- flat_samples + run_len
        }
        run_len <- 0L
      }
    }
    # Handle trailing run
    if (run_len > min_run) {
      flat_samples <- flat_samples + run_len
    }
    pct_flatline <- 100 * flat_samples / n

    # Artifact detection: gradient-based
    gradient_threshold <- 10 / sr
    pct_artifact <- 100 * sum(abs(diff(signal)) > gradient_threshold) / (n - 1)

    # SNR estimation using Donoho's estimator
    signal_power <- stats::var(signal)
    dd <- diff(diff(signal))
    noise_power <- stats::var(dd) / 6

    if (noise_power < .Machine$double.eps) {
      snr_db <- 60
    } else {
      snr_db <- 10 * log10(signal_power / noise_power)
      snr_db <- min(snr_db, 60)
    }

    # Quality score: start at 100, subtract penalties
    quality_score <- 100
    if (pct_negative > 5) quality_score <- quality_score - 30
    if (pct_flatline > 10) quality_score <- quality_score - 20
    if (pct_artifact > 5) quality_score <- quality_score - 20
    if (mean_sc < 0.05 || mean_sc > 50) quality_score <- quality_score - 15
    if (snr_db < 10) quality_score <- quality_score - 15
    quality_score <- max(0, min(100, quality_score))

    # Quality label
    quality_label <- if (quality_score >= 70) {
      "good"
    } else if (quality_score >= 40) {
      "acceptable"
    } else {
      "poor"
    }

    results[[ch]] <- data.frame(
      channel = ch_labels[ch],
      mean_sc = mean_sc,
      sd_sc = sd_sc,
      min_sc = min_sc,
      max_sc = max_sc,
      pct_negative = pct_negative,
      pct_flatline = pct_flatline,
      pct_artifact = pct_artifact,
      snr_db = snr_db,
      quality_score = quality_score,
      quality_label = quality_label,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}
