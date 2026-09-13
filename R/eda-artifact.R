#' Detect and Correct EDA Artifacts
#'
#' Identifies artifacts in electrodermal activity (EDA) signals using one or
#' more detection methods (threshold, gradient, flatline, and epoch-based
#' wavelet / machine-learning methods) and optionally corrects them via
#' interpolation or NA replacement. An accelerometer channel can additionally
#' gate out high-motion epochs (Taylor et al., 2015).
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param methods Character vector of detection methods to apply. Any
#'   combination of \code{"threshold"}, \code{"gradient"}, \code{"flatline"},
#'   \code{"wavelet"} (per-5 s-epoch stationary-wavelet fine-scale energy
#'   outliers), and \code{"ml"} (a shipped ridge-logistic classifier on
#'   wavelet/derivative epoch features). Default is the first three.
#' @param epoch_sec Epoch length in seconds for the \code{"wavelet"} /
#'   \code{"ml"} methods and accelerometer gating (default: 5).
#' @param wavelet_k Robust-outlier multiplier for the \code{"wavelet"} method;
#'   an epoch is flagged when its peak level-1 detail coefficient exceeds
#'   \code{median + wavelet_k * MAD} across epochs (default: 5).
#' @param ml_threshold Probability threshold for the \code{"ml"} method's binary
#'   label (default: 0.5).
#' @param acc_channel Optional channel index or label (or several) identifying
#'   accelerometer channel(s). When supplied, per-epoch motion energy (variance
#'   of the accelerometer magnitude) gates the EDA channels: epochs with energy
#'   above \code{acc_threshold} are flagged. The accelerometer channels are
#'   excluded from EDA artifact detection.
#' @param acc_threshold Motion-energy threshold for accelerometer gating. If
#'   \code{NULL} (default), computed robustly as \code{median + 3 * MAD} of the
#'   per-epoch motion energies.
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
#'     \item{epochs}{(only when \code{"wavelet"}/\code{"ml"} or accelerometer
#'       gating is used) A \code{data.frame}, one row per EDA channel per epoch,
#'       with the epoch bounds and, as applicable, \code{wavelet_flag},
#'       \code{ml_prob}, \code{ml_flag}, \code{motion_energy}, and \code{acc_flag}.}
#'   }
#'   If \code{correct != "none"}, the corrected signal is stored in the
#'   \code{output_assay}.
#' @references
#' Taylor, S., et al. (2015). "Automatic identification of artifacts in
#' electrodermal activity data." \emph{IEEE EMBC}, 1934-1937.
#' \doi{10.1109/EMBC.2015.7318762}
#'
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
                        epoch_sec = 5,
                        wavelet_k = 5,
                        ml_threshold = 0.5,
                        acc_channel = NULL,
                        acc_threshold = NULL,
                        correct = "interpolate",
                        assay_name = NULL,
                        output_assay = "cleaned") {
  stopifnot(inherits(x, "PhysioExperiment"))
  methods <- match.arg(methods,
                       c("threshold", "gradient", "flatline", "wavelet", "ml"),
                       several.ok = TRUE)
  stopifnot(is.numeric(threshold_range), length(threshold_range) == 2,
            threshold_range[1] < threshold_range[2])
  stopifnot(is.numeric(flatline_sec), flatline_sec > 0)
  stopifnot(is.numeric(epoch_sec), epoch_sec > 0)
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

  # Resolve accelerometer channel(s); they are excluded from EDA detection.
  acc_idx <- integer(0)
  if (!is.null(acc_channel)) {
    acc_idx <- if (is.character(acc_channel)) {
      match(acc_channel, ch_labels)
    } else {
      as.integer(acc_channel)
    }
    if (anyNA(acc_idx) || any(acc_idx < 1 | acc_idx > n_channels)) {
      stop("acc_channel does not match a channel in x", call. = FALSE)
    }
  }
  eda_channels <- setdiff(seq_len(n_channels), acc_idx)

  # Combined artifact mask (TRUE = artifact)
  mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)

  # Per-method masks for summary
  method_masks <- list()

  # Per-epoch record (populated by wavelet / ml / acc gating)
  epoch_rows <- list()

  # --- Threshold detection ---
  if ("threshold" %in% methods) {
    thresh_mask <- data < threshold_range[1] | data > threshold_range[2]
    if (length(acc_idx)) thresh_mask[, acc_idx] <- FALSE
    method_masks[["threshold"]] <- thresh_mask
    mask <- mask | thresh_mask
  }

  # --- Gradient detection ---
  if ("gradient" %in% methods) {
    grad_mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)
    for (ch in eda_channels) {
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
    for (ch in eda_channels) {
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

  # --- Epoch-based methods (wavelet / ml) and accelerometer gating ---
  use_wavelet <- "wavelet" %in% methods
  use_ml <- "ml" %in% methods
  use_acc <- length(acc_idx) > 0
  if (use_wavelet || use_ml || use_acc) {
    # Per-epoch accelerometer motion energy + flag (shared across EDA channels).
    acc_energy <- NULL
    acc_flag <- NULL
    if (use_acc) {
      ae <- .acc_motion_energy(data[, acc_idx, drop = FALSE], sr, epoch_sec)
      acc_energy <- ae$motion_energy
      acc_thr <- if (is.null(acc_threshold)) {
        stats::median(acc_energy) + 3 * stats::mad(acc_energy)
      } else acc_threshold
      acc_flag <- acc_energy > acc_thr
    }

    wave_mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)
    ml_mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)
    acc_mask <- matrix(FALSE, nrow = n_time, ncol = n_channels)
    ml_model <- if (use_ml) .load_eda_artifact_model() else NULL

    for (ch in eda_channels) {
      feat <- .eda_epoch_features(data[, ch], sr, epoch_sec)
      n_ep <- nrow(feat)
      w_flag <- if (use_wavelet) .eda_wavelet_flag(feat, wavelet_k) else rep(NA, n_ep)
      m_prob <- if (use_ml) .eda_ml_predict(feat, ml_model) else rep(NA_real_, n_ep)
      m_flag <- if (use_ml) m_prob > ml_threshold else rep(NA, n_ep)

      expand_epochs <- function(flag) {
        col <- logical(n_time)
        for (e in which(flag %in% TRUE)) col[feat$start[e]:feat$end[e]] <- TRUE
        col
      }
      if (use_wavelet) wave_mask[, ch] <- expand_epochs(w_flag)
      if (use_ml) ml_mask[, ch] <- expand_epochs(m_flag)
      if (use_acc) acc_mask[, ch] <- expand_epochs(acc_flag)

      epoch_rows[[length(epoch_rows) + 1L]] <- data.frame(
        channel = ch_labels[ch],
        epoch = feat$epoch,
        start_sec = feat$start_sec,
        end_sec = feat$end_sec,
        wavelet_flag = w_flag,
        ml_prob = m_prob,
        ml_flag = m_flag,
        motion_energy = if (use_acc) acc_energy else NA_real_,
        acc_flag = if (use_acc) acc_flag else NA,
        stringsAsFactors = FALSE
      )
    }

    if (use_wavelet) { method_masks[["wavelet"]] <- wave_mask; mask <- mask | wave_mask }
    if (use_ml) { method_masks[["ml"]] <- ml_mask; mask <- mask | ml_mask }
    if (use_acc) { method_masks[["acc"]] <- acc_mask; mask <- mask | acc_mask }
  }

  # --- Build summary (over EDA channels and every method actually run) ---
  run_methods <- names(method_masks)
  summary_rows <- vector("list", length(eda_channels) * length(run_methods))
  idx <- 1L
  for (ch in eda_channels) {
    for (m in run_methods) {
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
  epochs_df <- if (length(epoch_rows)) do.call(rbind, epoch_rows) else NULL
  meta <- S4Vectors::metadata(x)
  meta$eda_artifacts <- list(
    mask = mask,
    summary = summary_df,
    epochs = epochs_df
  )
  S4Vectors::metadata(x) <- meta

  # --- Correction (EDA channels only) ---
  if (correct != "none") {
    corrected <- data
    for (ch in eda_channels) {
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
