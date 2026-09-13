# Sympathetic skin response (SSR): a clinical autonomic test in which a
# stimulus (electrical, acoustic or inspiratory) evokes a transient skin
# conductance deflection. edaSSR() extracts per-stimulus latency, amplitude and
# morphology, classifies present/absent, and summarises habituation and the
# clinical absent-SSR flag (Vetrugno et al. 2003).

#' Detrend a response window by the pre-stimulus tonic trend
#'
#' Fits the tonic slope + level from the `base_n` samples immediately before the
#' stimulus and subtracts that (extrapolated) linear trend from the response
#' window, leaving the phasic residual (0 at the tonic baseline).
#' @keywords internal
#' @noRd
.ssr_detrend <- function(sig, s, s_start, s_end, base_n) {
  pb_start <- max(1L, s - base_n)
  pb_end <- max(pb_start, s - 1L)
  win <- s_start:s_end
  if (pb_end > pb_start) {
    pt <- (pb_start:pb_end) - s                 # centred on the stimulus
    y <- sig[pb_start:pb_end]
    vt <- stats::var(pt)
    b <- if (vt > 0) stats::cov(pt, y) / vt else 0
    a <- mean(y) - b * mean(pt)
  } else {
    a <- sig[pb_start]; b <- 0
  }
  sig[win] - (a + b * (win - s))
}

#' Resolve stimulus events to sample indices
#' @keywords internal
#' @noRd
.ssr_event_samples <- function(events, n, sr, unit) {
  events <- as.numeric(events)
  if (length(events) == 0L || any(!is.finite(events))) {
    stop("`stimulus_events` must be a non-empty finite numeric vector.",
         call. = FALSE)
  }
  is_index <- switch(unit,
    index = TRUE, seconds = FALSE,
    auto = all(events == round(events)) && all(events >= 1) && max(events) <= n)
  idx <- if (is_index) as.integer(round(events)) else as.integer(round(events * sr) + 1L)
  idx
}

#' Sympathetic skin response (SSR) analysis
#'
#' Extracts the sympathetic skin response evoked by each stimulus from a skin
#' conductance recording: the response onset latency, the baseline-to-peak
#' amplitude, the waveform morphology (monophasic / biphasic) and a
#' present/absent classification. A per-channel summary reports the response
#' rate, mean latency and amplitude, an amplitude-habituation slope across the
#' repeated stimuli, and the clinical absent-SSR flag.
#'
#' The tonic trend estimated from the pre-stimulus baseline is subtracted from
#' each response window, so slow skin-conductance drift is neither counted as a
#' response nor mistaken for a biphasic deflection. In the detrended window the
#' phasic peak is located; a response is accepted only if it exceeds
#' `amplitude_min` and recovers (rises then falls), and the latency is measured
#' from the stimulus to the onset of that phasic deflection (Vetrugno et al.
#' 2003).
#'
#' @param x A PhysioExperiment object with a skin conductance assay.
#' @param stimulus_events Stimulus times, as sample indices or times in seconds
#'   (see `events_unit`).
#' @param assay_name Input assay (default: `defaultAssay(x)`).
#' @param response_window Numeric `c(start, end)` seconds after each stimulus in
#'   which the SSR is sought (default `c(0.5, 8)`).
#' @param baseline_sec Length in seconds of the pre-stimulus window used to
#'   estimate and remove the tonic trend (default 2).
#' @param amplitude_min Minimum detrended peak amplitude (uS) for a response to
#'   count as present (default 0.05).
#' @param recovery_frac Fraction of the peak amplitude the trace must fall back
#'   after the peak for the response to be accepted as phasic (rather than a
#'   monotone drift; default 0.1).
#' @param onset_frac Fraction of the SSR amplitude above baseline at which the
#'   response onset (and hence latency) is marked (default 0.05).
#' @param morphology_frac Fraction of the SSR amplitude an opposite-going
#'   deflection must reach for the response to be classed biphasic (default 0.3).
#' @param events_unit Unit of `stimulus_events`: `"auto"` (default), `"index"`
#'   or `"seconds"`.
#'
#' @return An `eda_ssr` object with `responses` (one row per stimulus per
#'   channel: `channel`, `stimulus`, `stim_time`, `present`, `latency`,
#'   `amplitude`, `peak_time`, `morphology`) and `summary` (one row per channel:
#'   `channel`, `n_stimuli`, `n_present`, `response_rate`, `mean_latency`,
#'   `mean_amplitude`, `habituation_slope`, `ssr_absent`).
#' @references Vetrugno, R. et al. (2003). Sympathetic skin response: basic
#'   mechanisms and clinical applications. Clin Auton Res 13(4):256-270.
#' @seealso [edaSympatheticIndex()]
#' @export
#' @examples
#' # 3 stimuli at 10, 40 and 70 s, each followed by an SCR
#' pe <- edaSimulate(n_time = 900, sr = 10, scr_count = 0, seed = 1)
#' sig <- as.numeric(SummarizedExperiment::assay(pe))
#' sr <- 10
#' for (st in c(10, 40, 70)) {
#'   t <- (seq_along(sig) - 1) / sr - st
#'   sig <- sig + 0.5 * ifelse(t > 0, (exp(-t / 3) - exp(-t / 0.75)), 0)
#' }
#' SummarizedExperiment::assay(pe) <- matrix(sig, ncol = 1)
#' res <- edaSSR(pe, stimulus_events = c(10, 40, 70), events_unit = "seconds")
#' res$summary
edaSSR <- function(x, stimulus_events, assay_name = NULL,
                   response_window = c(0.5, 8), baseline_sec = 2,
                   amplitude_min = 0.05, recovery_frac = 0.1,
                   onset_frac = 0.05, morphology_frac = 0.3,
                   events_unit = c("auto", "index", "seconds")) {
  stopifnot(inherits(x, "PhysioExperiment"))
  events_unit <- match.arg(events_unit)
  if (!is.numeric(response_window) || length(response_window) != 2L ||
      any(!is.finite(response_window)) ||
      response_window[1] < 0 || response_window[2] <= response_window[1]) {
    stop("`response_window` must be c(start, end) with 0 <= start < end.",
         call. = FALSE)
  }
  stopifnot(is.numeric(baseline_sec), length(baseline_sec) == 1L,
            is.finite(baseline_sec), baseline_sec > 0)
  stopifnot(is.numeric(amplitude_min), length(amplitude_min) == 1L,
            is.finite(amplitude_min), amplitude_min >= 0)
  stopifnot(is.numeric(recovery_frac), length(recovery_frac) == 1L,
            is.finite(recovery_frac), recovery_frac >= 0, recovery_frac < 1)
  stopifnot(is.numeric(onset_frac), length(onset_frac) == 1L,
            is.finite(onset_frac), onset_frac >= 0, onset_frac < 1)
  stopifnot(is.numeric(morphology_frac), length(morphology_frac) == 1L,
            is.finite(morphology_frac), morphology_frac > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  if (!is.numeric(data) || any(!is.finite(data)) || nrow(data) < 2L ||
      ncol(data) < 1L) {
    stop("the skin conductance assay must be a finite numeric matrix with ",
         ">= 2 time points and >= 1 channel.", call. = FALSE)
  }
  sr <- samplingRate(x)
  n_time <- nrow(data); n_channels <- ncol(data)

  stim <- .ssr_event_samples(stimulus_events, n_time, sr, events_unit)
  if (any(stim < 1L | stim > n_time)) {
    stop(sprintf(paste0("`stimulus_events` resolve to sample indices outside ",
                        "[1, %d]; check `events_unit`."), n_time), call. = FALSE)
  }
  n_stim <- length(stim)
  ch_info <- SummarizedExperiment::colData(x)
  ch_names <- if ("label" %in% colnames(ch_info)) {
    as.character(ch_info$label)
  } else {
    paste0("Ch", seq_len(n_channels))
  }
  w0 <- round(response_window[1] * sr)
  w1 <- round(response_window[2] * sr)
  base_n <- max(2L, as.integer(round(baseline_sec * sr)))

  resp <- list()
  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]
    for (k in seq_len(n_stim)) {
      s <- stim[k]
      s_start <- s + w0
      s_end <- s + w1
      # do not let the window run into the next stimulus
      later <- stim[stim > s]
      if (length(later)) s_end <- min(s_end, min(later) - 1L)
      s_start <- max(1L, s_start); s_end <- min(n_time, s_end)

      present <- FALSE; latency <- NA_real_; amplitude <- NA_real_
      peak_time <- NA_real_; morphology <- NA_character_
      if (s_end - s_start >= 1L && s >= 1L && s <= n_time) {
        # remove the tonic trend estimated from the pre-stimulus baseline, so
        # slow drift is not mistaken for (or masked as) a phasic response
        detr <- .ssr_detrend(sig, s, s_start, s_end, base_n)
        peak_local <- which.max(detr)
        peak_sample <- s_start + peak_local - 1L
        amplitude <- detr[peak_local]              # relative to the tonic trend
        post <- detr[peak_local:length(detr)]
        recovery <- amplitude - min(post)          # a real SSR rises then falls
        present <- amplitude >= amplitude_min &&
          recovery >= recovery_frac * amplitude
        if (present) {
          # onset = last pre-peak sample still at/below onset_frac of the amplitude
          thr <- onset_frac * amplitude
          pre <- detr[seq_len(peak_local)]
          below <- which(pre <= thr)
          onset_sample <- s_start + (if (length(below)) max(below) else 1L) - 1L
          latency <- (onset_sample - s) / sr
          peak_time <- (peak_sample - s) / sr
          # biphasic if the detrended trace dips a fraction of the amplitude
          # below the tonic baseline after the peak
          morphology <- if (min(post) < -morphology_frac * amplitude) {
            "biphasic"
          } else {
            "monophasic"
          }
        }
      }
      resp[[length(resp) + 1L]] <- data.frame(
        channel = ch_names[ch], stimulus = k, stim_time = (s - 1) / sr,
        present = present, latency = latency, amplitude = amplitude,
        peak_time = peak_time, morphology = morphology,
        stringsAsFactors = FALSE)
    }
  }
  responses <- do.call(rbind, resp)

  # per-channel clinical summary. Subset by row position (responses are built
  # channel-outer, stimulus-inner) rather than by label, so channels that share
  # a label are never merged.
  summ <- lapply(seq_len(n_channels), function(ch) {
    r <- responses[((ch - 1L) * n_stim + 1L):(ch * n_stim), , drop = FALSE]
    r <- r[order(r$stim_time), ]                 # chronological, not listed, order
    pres <- r$present
    hab <- NA_real_
    if (sum(pres) >= 2L) {
      amp <- r$amplitude[pres]; ord <- which(pres)
      hab <- unname(stats::coef(stats::lm(amp ~ ord))[2])   # amplitude trend
    }
    data.frame(
      channel = ch_names[ch], n_stimuli = nrow(r), n_present = sum(pres),
      response_rate = mean(pres),
      mean_latency = if (any(pres)) mean(r$latency[pres]) else NA_real_,
      mean_amplitude = if (any(pres)) mean(r$amplitude[pres]) else NA_real_,
      habituation_slope = hab, ssr_absent = sum(pres) == 0L,
      stringsAsFactors = FALSE)
  })
  out <- list(responses = responses, summary = do.call(rbind, summ))
  class(out) <- "eda_ssr"
  out
}

#' @export
print.eda_ssr <- function(x, ...) {
  cat(sprintf("<eda_ssr> %d channel(s), %d stimulus/stimuli\n",
              nrow(x$summary), max(x$responses$stimulus)))
  for (i in seq_len(nrow(x$summary))) {
    s <- x$summary[i, ]
    cat(sprintf("  %s: %d/%d responses (%.0f%%)%s",
                s$channel, s$n_present, s$n_stimuli, 100 * s$response_rate,
                if (s$ssr_absent) "  [ABSENT SSR]" else ""))
    if (s$n_present > 0) {
      cat(sprintf(", latency %.2f s, amplitude %.3f uS",
                  s$mean_latency, s$mean_amplitude))
    }
    cat("\n")
  }
  invisible(x)
}
