# Wavelet / ML / accelerometer helpers for edaArtifact() (WS5-28).

#' Undecimated (Stationary) Haar Wavelet Transform
#'
#' A dependency-free maximal-overlap ("a trous") Haar wavelet transform. Every
#' returned component has the same length as \code{x} (undecimated), so the
#' detail coefficients are time-aligned with the samples - ideal for locating
#' high-frequency artifacts. The normalisation matches
#' \code{waveslim::modwt(x, "haar")}: the sum of squared coefficients over all
#' levels plus the final smooth equals the energy of \code{x}.
#'
#' @param x Numeric signal vector.
#' @param J Number of decomposition levels (default: 4).
#' @return A named list \code{d1..dJ} (level detail coefficients) and \code{sJ}
#'   (level-J smooth), each the same length as \code{x}.
#' @keywords internal
.swt_haar <- function(x, J = 4) {
  n <- length(x)
  s <- x
  out <- vector("list", J + 1)
  for (j in seq_len(J)) {
    lag <- 2^(j - 1)
    idx <- ((seq_len(n) - 1 - lag) %% n) + 1     # circular shift by 2^(j-1)
    d <- (s - s[idx]) / 2                          # level-j detail (MODWT norm)
    s <- (s + s[idx]) / 2                          # level-j smooth
    out[[j]] <- d
  }
  out[[J + 1]] <- s
  names(out) <- c(paste0("d", seq_len(J)), paste0("s", J))
  out
}

#' Non-overlapping Epoch Boundaries
#'
#' @param n Number of samples.
#' @param sr Sampling rate in Hz.
#' @param epoch_sec Epoch length in seconds.
#' @return A \code{data.frame} with columns \code{epoch}, \code{start},
#'   \code{end} (sample indices), \code{start_sec}, \code{end_sec}.
#' @keywords internal
.eda_epoch_bounds <- function(n, sr, epoch_sec = 5) {
  L <- max(2L, as.integer(round(epoch_sec * sr)))
  starts <- seq(1L, n, by = L)
  ends <- pmin(starts + L - 1L, n)
  data.frame(epoch = seq_along(starts), start = starts, end = ends,
             start_sec = (starts - 1) / sr, end_sec = (ends - 1) / sr)
}

#' Per-epoch Artifact Features for EDA
#'
#' Computes, for each non-overlapping epoch, the derivative statistics and the
#' undecimated-Haar wavelet detail statistics used by the wavelet rule and the
#' logistic model (EDA-Explorer-style features; Taylor et al., 2015).
#'
#' @param sig Numeric EDA signal (single channel).
#' @param sr Sampling rate in Hz.
#' @param epoch_sec Epoch length in seconds (default: 5).
#' @param J Wavelet levels (default: 4).
#' @return A \code{data.frame} with one row per epoch: epoch bounds plus
#'   \code{amp}, \code{maxabsd}, \code{avgabsd}, \code{maxabs2d},
#'   \code{avgabs2d}, \code{wav_max_d1}, \code{wav_std_d1}, \code{wav_max_d2},
#'   \code{wav_std_d2}, \code{wav_std_d3}.
#' @keywords internal
.eda_epoch_features <- function(sig, sr, epoch_sec = 5, J = 4) {
  n <- length(sig)
  eb <- .eda_epoch_bounds(n, sr, epoch_sec)
  sw <- .swt_haar(sig, J)
  d1 <- sw$d1; d2 <- sw$d2; d3 <- sw$d3
  feats <- lapply(seq_len(nrow(eb)), function(i) {
    rng <- eb$start[i]:eb$end[i]
    se <- sig[rng]
    dd <- if (length(se) > 1) diff(se) else 0
    dd2 <- if (length(se) > 2) diff(se, differences = 2) else 0
    sd0 <- function(v) if (length(v) > 1) stats::sd(v) else 0
    data.frame(
      amp = mean(se),
      maxabsd = max(abs(dd)), avgabsd = mean(abs(dd)),
      maxabs2d = max(abs(dd2)), avgabs2d = mean(abs(dd2)),
      wav_max_d1 = max(abs(d1[rng])), wav_std_d1 = sd0(d1[rng]),
      wav_max_d2 = max(abs(d2[rng])), wav_std_d2 = sd0(d2[rng]),
      wav_std_d3 = sd0(d3[rng])
    )
  })
  cbind(eb, do.call(rbind, feats))
}

#' Wavelet-rule Artifact Flag per Epoch
#'
#' Flags an epoch as an artifact when its fine-scale wavelet energy
#' (\code{wav_max_d1}, the peak level-1 detail coefficient) is a robust outlier:
#' greater than \code{median + k * MAD} across epochs. Motion and contact
#' artifacts inject high-frequency energy that genuine (slow) skin conductance
#' does not.
#'
#' @param features Epoch features from \code{\link{.eda_epoch_features}}.
#' @param k Robust-outlier multiplier (default: 5).
#' @return Logical vector, one per epoch (\code{TRUE} = artifact).
#' @keywords internal
.eda_wavelet_flag <- function(features, k = 5) {
  stat <- features$wav_max_d1
  thr <- stats::median(stat) + k * stats::mad(stat)
  stat > thr
}

# Package-level cache for the shipped logistic model.
.eda_model_cache <- new.env(parent = emptyenv())

#' Load the Shipped EDA Artifact Logistic Model
#'
#' Reads the lightweight ridge-logistic artifact model shipped in
#' \code{inst/extdata/eda_artifact_model.csv} (columns \code{feature},
#' \code{mean}, \code{sd}, \code{coef}); the result is cached.
#'
#' @return A \code{data.frame} with the model coefficients and standardisation.
#' @keywords internal
.load_eda_artifact_model <- function() {
  if (!is.null(.eda_model_cache$model)) return(.eda_model_cache$model)
  path <- system.file("extdata", "eda_artifact_model.csv", package = "PhysioEDA")
  if (path == "") {
    stop("eda_artifact_model.csv not found; is PhysioEDA installed correctly?",
         call. = FALSE)
  }
  m <- utils::read.csv(path, stringsAsFactors = FALSE)
  .eda_model_cache$model <- m
  m
}

#' Per-epoch Artifact Probability from the Logistic Model
#'
#' Applies the shipped standardised ridge-logistic model to the epoch features.
#'
#' @param features Epoch features from \code{\link{.eda_epoch_features}}.
#' @param model The model \code{data.frame} from
#'   \code{\link{.load_eda_artifact_model}}.
#' @return Numeric vector of artifact probabilities, one per epoch.
#' @keywords internal
.eda_ml_predict <- function(features, model) {
  is_int <- model$feature == "(Intercept)"
  eta <- model$coef[is_int]
  fm <- model[!is_int, , drop = FALSE]
  for (i in seq_len(nrow(fm))) {
    fname <- fm$feature[i]
    z <- (features[[fname]] - fm$mean[i]) / fm$sd[i]
    eta <- eta + fm$coef[i] * z
  }
  stats::plogis(eta)
}

#' Per-epoch Accelerometer Motion Energy
#'
#' Motion energy per epoch: the variance of the accelerometer magnitude within
#' the epoch. A still sensor sees near-constant magnitude (low variance); motion
#' raises it.
#'
#' @param acc Numeric matrix (time x axes) of accelerometer channels.
#' @param sr Sampling rate in Hz.
#' @param epoch_sec Epoch length in seconds (default: 5).
#' @return A \code{data.frame} of epoch bounds plus \code{motion_energy}.
#' @keywords internal
.acc_motion_energy <- function(acc, sr, epoch_sec = 5) {
  acc <- as.matrix(acc)
  mag <- if (ncol(acc) > 1) sqrt(rowSums(acc^2)) else acc[, 1]
  eb <- .eda_epoch_bounds(length(mag), sr, epoch_sec)
  energy <- vapply(seq_len(nrow(eb)), function(i) {
    seg <- mag[eb$start[i]:eb$end[i]]
    if (length(seg) > 1) stats::var(seg) else 0
  }, numeric(1))
  eb$motion_energy <- energy
  eb
}
