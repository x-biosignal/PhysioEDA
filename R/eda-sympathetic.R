# Electrodermal sympathetic-tone indices: EDASymp (spectral power in the
# sympathetic band) and TVSymp (time-varying sympathetic index via complex
# demodulation), after Posada-Quintero et al. (2016).

#' Validate an optional positive finite scalar (rate / length parameter)
#' @keywords internal
#' @noRd
.eda_pos_scalar <- function(v, name, allow_null = TRUE) {
  if (is.null(v)) {
    if (allow_null) return(invisible())
    stop(sprintf("`%s` must be a single positive number.", name), call. = FALSE)
  }
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
    stop(sprintf("`%s` must be a single positive number.", name), call. = FALSE)
  }
  invisible()
}

#' Validate the sympathetic band and the electrodermal signal matrix
#' @keywords internal
#' @noRd
.eda_check_band <- function(band) {
  if (!is.numeric(band) || length(band) != 2L || any(!is.finite(band)) ||
      band[1] < 0 || band[2] <= band[1]) {
    stop("`band` must be c(low, high) with 0 <= low < high.", call. = FALSE)
  }
}

#' Resample a signal to a target rate by linear interpolation
#' @keywords internal
#' @noRd
.eda_resample <- function(x, fs_in, fs_out) {
  if (is.null(fs_out) || fs_out >= fs_in) {
    return(list(x = x, fs = fs_in))
  }
  n_in <- length(x)
  dur <- (n_in - 1) / fs_in
  n_out <- max(2L, as.integer(floor(dur * fs_out)) + 1L)
  t_in <- (seq_len(n_in) - 1) / fs_in
  t_out <- (seq_len(n_out) - 1) / fs_out
  list(x = stats::approx(t_in, x, xout = t_out, rule = 2)$y, fs = fs_out)
}

#' Welch power spectral density (Hann taper, 50% overlap)
#' @keywords internal
#' @noRd
.eda_welch <- function(x, fs, seg_len = NULL) {
  n <- length(x)
  if (n < 2L) return(list(freq = 0, psd = 0))
  if (is.null(seg_len) || seg_len > n) seg_len <- n
  seg_len <- max(2L, min(as.integer(seg_len), n))     # never exceed the record
  step <- max(1L, seg_len %/% 2L)
  win <- 0.5 - 0.5 * cos(2 * pi * (seq_len(seg_len) - 1) / (seg_len - 1))  # Hann
  u <- sum(win^2)                                     # window power
  if (u <= 0) { win <- rep(1, seg_len); u <- seg_len }  # boxcar fallback
  starts <- seq(1L, n - seg_len + 1L, by = step)
  if (!length(starts)) starts <- 1L
  nfreq <- seg_len %/% 2L + 1L
  # bins to double for the one-sided spectrum: all above DC except a genuine
  # unique Nyquist bin (present only for even seg_len)
  hi <- if (seg_len %% 2L == 0L) nfreq - 1L else nfreq
  psd <- numeric(nfreq)
  for (s in starts) {
    seg <- x[s:(s + seg_len - 1L)]
    seg <- (seg - mean(seg)) * win
    X <- stats::fft(seg)[seq_len(nfreq)]
    p <- (Mod(X)^2) / (fs * u)
    if (hi >= 2L) p[2:hi] <- 2 * p[2:hi]
    psd <- psd + p
  }
  psd <- psd / length(starts)
  freq <- (seq_len(nfreq) - 1) * fs / seg_len
  list(freq = freq, psd = psd)
}

#' Trapezoidal band power of a PSD over [lo, hi]
#' @keywords internal
#' @noRd
.eda_bandpower <- function(freq, psd, lo, hi) {
  keep <- freq >= lo & freq <= hi
  if (sum(keep) < 2L) return(0)
  f <- freq[keep]; p <- psd[keep]
  sum(diff(f) * (utils::head(p, -1) + utils::tail(p, -1)) / 2)
}

#' Electrodermal sympathetic index (EDASymp)
#'
#' Spectral power of the electrodermal signal in the sympathetic band
#' (0.045-0.25 Hz by default), following Posada-Quintero et al. (2016). The
#' signal is optionally downsampled, its Welch power spectral density is
#' estimated, and the band power is reported both in absolute units and
#' normalised by the total spectral power.
#'
#' @param x A PhysioExperiment object with an electrodermal assay.
#' @param band Numeric `c(low, high)` sympathetic band in Hz (default
#'   `c(0.045, 0.25)`).
#' @param assay_name Input assay (default: `defaultAssay(x)`).
#' @param resample_hz Rate (Hz) the signal is downsampled to before spectral
#'   estimation (default 2; `NULL` keeps the native rate). Ignored if it is not
#'   below the native sampling rate.
#' @param seg_sec Welch segment length in seconds (`NULL` = one whole-signal
#'   Hann window).
#'
#' @return A data frame with one row per channel: `channel`, `edasymp`
#'   (absolute band power), `edasymp_norm` (band power / total power, in
#'   `[0, 1]`, or `NA` for a constant / zero-power channel), `total_power`, and
#'   the `band_low`/`band_high` used.
#' @references Posada-Quintero, H.F. et al. (2016). Power spectral density
#'   analysis of electrodermal activity for sympathetic function assessment.
#'   Ann Biomed Eng 44(10):3124-3135.
#' @seealso [edaTVSymp()]
#' @export
#' @examples
#' pe <- edaSimulate(n_time = 1200, sr = 4, scr_count = 6, seed = 1)
#' edaSympatheticIndex(pe)
edaSympatheticIndex <- function(x, band = c(0.045, 0.25), assay_name = NULL,
                                resample_hz = 2, seg_sec = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  .eda_check_band(band)
  .eda_pos_scalar(resample_hz, "resample_hz")
  .eda_pos_scalar(seg_sec, "seg_sec")
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (!is.numeric(data) || any(!is.finite(data))) {
    stop("the electrodermal signal must be a finite numeric matrix.",
         call. = FALSE)
  }
  if (nrow(data) < 2L || ncol(data) < 1L) {
    stop("the electrodermal assay must have >= 2 time points and >= 1 channel.",
         call. = FALSE)
  }
  sr <- samplingRate(x)
  eff_fs <- if (is.null(resample_hz) || resample_hz >= sr) sr else resample_hz
  if (band[2] >= eff_fs / 2) {
    stop("`band` upper edge must be below the (effective) Nyquist frequency.",
         call. = FALSE)
  }

  rows <- lapply(seq_len(ncol(data)), function(ch) {
    rs <- .eda_resample(as.numeric(data[, ch]), sr, resample_hz)
    seg_len <- if (is.null(seg_sec)) NULL else as.integer(round(seg_sec * rs$fs))
    sp <- .eda_welch(rs$x, rs$fs, seg_len)
    band_p <- .eda_bandpower(sp$freq, sp$psd, band[1], band[2])
    total_p <- .eda_bandpower(sp$freq, sp$psd, sp$freq[2], max(sp$freq))
    data.frame(
      channel = ch, edasymp = band_p,
      edasymp_norm = if (total_p > 0) band_p / total_p else NA_real_,
      total_power = total_p, band_low = band[1], band_high = band[2],
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Reflect-padded low-pass by a Hann-weighted moving average
#' @keywords internal
#' @noRd
.eda_lowpass <- function(v, len) {
  len <- max(3L, as.integer(len))
  if (len %% 2L == 0L) len <- len + 1L
  k <- 0.5 - 0.5 * cos(2 * pi * (seq_len(len) - 1) / (len - 1))
  k <- k / sum(k)
  h <- len %/% 2L
  n <- length(v)
  if (n < 2L) return(v)
  pad <- c(rev(v[2:min(h + 1L, n)]), v, rev(v[max(1L, n - h):(n - 1L)]))
  # pad may be short if h >= n; fall back to edge replication
  need <- n + 2L * h
  if (length(pad) < need) {
    pad <- c(rep(v[1], h), v, rep(v[n], h))
  } else if (length(pad) > need) {
    pad <- pad[seq_len(need)]
  }
  as.numeric(stats::filter(pad, k, sides = 2L))[(h + 1L):(h + n)]
}

#' Instantaneous amplitude of the band component at fc by complex demodulation
#' @keywords internal
#' @noRd
.eda_cdm_amp <- function(x, fs, fc, lp_hz) {
  t <- (seq_along(x) - 1) / fs
  z <- x * exp(-2i * pi * fc * t)               # shift fc to baseband
  len <- max(3L, round(fs / lp_hz))
  re <- .eda_lowpass(Re(z), len)
  im <- .eda_lowpass(Im(z), len)
  2 * sqrt(re^2 + im^2)                          # instantaneous amplitude
}

#' Time-varying electrodermal sympathetic index (TVSymp)
#'
#' A time-resolved sympathetic-tone index built by complex demodulation of the
#' dominant components of the electrodermal signal in the sympathetic band
#' (Posada-Quintero et al. 2016). The two strongest spectral components in the
#' band are demodulated to their instantaneous amplitudes, combined, and
#' normalised to a reproducible trajectory that rises when sympathetic
#' oscillations appear.
#'
#' @param x A PhysioExperiment object with an electrodermal assay.
#' @param band Numeric `c(low, high)` sympathetic band in Hz (default
#'   `c(0.08, 0.24)`).
#' @param assay_name Input assay (default: `defaultAssay(x)`).
#' @param resample_hz Rate (Hz) the signal is downsampled to (default 2; `NULL`
#'   keeps the native rate).
#' @param n_components Number of dominant band components to demodulate
#'   (default 2).
#' @param lp_hz Low-pass cutoff (Hz) of the demodulation envelope (default 0.5
#'   of the band width).
#'
#' @return A data frame with one row per time sample per channel: `channel`,
#'   `time_sec`, `tvsymp` (combined instantaneous amplitude) and `tvsymp_norm`
#'   (unit-variance normalised). The dominant component frequencies are attached
#'   as the `"components"` attribute.
#' @references Posada-Quintero, H.F. et al. (2016). Time-varying analysis of
#'   electrodermal activity during exercise. PLoS ONE 11(8):e0161553.
#' @seealso [edaSympatheticIndex()]
#' @export
#' @examples
#' pe <- edaSimulate(n_time = 1200, sr = 4, scr_count = 6, seed = 1)
#' tv <- edaTVSymp(pe)
#' head(tv)
edaTVSymp <- function(x, band = c(0.08, 0.24), assay_name = NULL,
                      resample_hz = 2, n_components = 2, lp_hz = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  .eda_check_band(band)
  stopifnot(is.numeric(n_components), length(n_components) == 1L,
            is.finite(n_components), n_components >= 1)
  .eda_pos_scalar(resample_hz, "resample_hz")
  .eda_pos_scalar(lp_hz, "lp_hz")
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (!is.numeric(data) || any(!is.finite(data))) {
    stop("the electrodermal signal must be a finite numeric matrix.",
         call. = FALSE)
  }
  if (nrow(data) < 2L || ncol(data) < 1L) {
    stop("the electrodermal assay must have >= 2 time points and >= 1 channel.",
         call. = FALSE)
  }
  sr <- samplingRate(x)
  eff_fs <- if (is.null(resample_hz) || resample_hz >= sr) sr else resample_hz
  if (band[2] >= eff_fs / 2) {
    stop("`band` upper edge must be below the (effective) Nyquist frequency.",
         call. = FALSE)
  }
  if (is.null(lp_hz)) lp_hz <- 0.5 * (band[2] - band[1])

  comp_all <- vector("list", ncol(data))
  rows <- lapply(seq_len(ncol(data)), function(ch) {
    rs <- .eda_resample(as.numeric(data[, ch]), sr, resample_hz)
    sig <- rs$x - mean(rs$x)
    sp <- .eda_welch(sig, rs$fs)
    inband <- which(sp$freq >= band[1] & sp$freq <= band[2])
    # dominant in-band component frequencies (fall back to band centre)
    if (length(inband)) {
      ord <- inband[order(sp$psd[inband], decreasing = TRUE)]
      fcs <- sp$freq[utils::head(ord, as.integer(n_components))]
    } else {
      fcs <- mean(band)
    }
    comp_all[[ch]] <<- fcs
    amp <- vapply(fcs, function(fc) .eda_cdm_amp(sig, rs$fs, fc, lp_hz),
                  numeric(length(sig)))
    if (is.null(dim(amp))) amp <- matrix(amp, ncol = 1)
    tv <- sqrt(rowSums(amp^2))
    sdv <- stats::sd(tv)
    data.frame(
      channel = ch, time_sec = (seq_along(tv) - 1) / rs$fs,
      tvsymp = tv,
      tvsymp_norm = if (is.finite(sdv) && sdv > 0) tv / sdv else tv * 0,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "components") <- comp_all
  out
}
