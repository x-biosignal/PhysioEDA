#' Event-Related SCR Analysis
#'
#' Performs event-related skin conductance response (ER-SCR) analysis by
#' extracting SCR features time-locked to experimental events. For each event,
#' the function searches for an SCR onset within a specified window and, if
#' found, computes amplitude, latency, rise time, and 50 percent recovery time.
#'
#' @param x A PhysioExperiment object containing EDA data and events.
#' @param event_type Character string specifying the event type to analyze.
#'   If NULL (default), all events are used.
#' @param onset_window Numeric vector of length 2 giving the valid SCR onset
#'   latency window in seconds relative to event onset (default: c(1, 4)).
#' @param peak_window Numeric vector of length 2 giving the minimum and maximum
#'   time from SCR onset to peak in seconds (default: c(0.5, 5)).
#' @param amplitude_min Minimum SCR amplitude in microsiemens to be considered
#'   a valid response (default: 0.01).
#' @param assay_name Name of the assay to use. If NULL, uses "phasic" if
#'   available, otherwise the default assay.
#' @return A data.frame with one row per event per channel containing columns:
#'   \describe{
#'     \item{event_index}{Integer index of the event.}
#'     \item{event_onset}{Numeric onset time of the event in seconds.}
#'     \item{channel}{Character name of the channel.}
#'     \item{scr_present}{Logical indicating whether a valid SCR was detected.}
#'     \item{scr_amplitude}{Numeric SCR amplitude (peak minus onset), or NA.}
#'     \item{scr_latency}{Numeric latency from event to SCR onset in seconds, or NA.}
#'     \item{scr_rise_time}{Numeric time from SCR onset to peak in seconds, or NA.}
#'     \item{scr_recovery_time}{Numeric time for 50 percent recovery from peak, or NA.}
#'   }
#' @references
#' Bach, D.R., et al. (2010). "Modelling event-related skin conductance
#' responses." \emph{International Journal of Psychophysiology}, 75(3),
#' 349-356. \doi{10.1016/j.ijpsycho.2010.01.005}
#' @seealso \code{\link{edaPeaks}} for general SCR peak detection,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition (run first),
#'   \code{\link{edaFeatures}} for summary feature extraction,
#'   \code{\link{plotPeaks}} for visualizing SCR peaks
#' @export
edaErscr <- function(x, event_type = NULL, onset_window = c(1, 4),
                     peak_window = c(0.5, 5), amplitude_min = 0.01,
                     assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  # Get events
  events <- getEvents(x, type = event_type)
  event_df <- events@events
  if (nrow(event_df) == 0) {
    stop("No events found", call. = FALSE)
  }

  # Resolve assay name
  if (is.null(assay_name)) {
    assay_names <- SummarizedExperiment::assayNames(x)
    if ("phasic" %in% assay_names) {
      assay_name <- "phasic"
    } else {
      assay_name <- defaultAssay(x)
    }
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Channel names
  ch_info <- SummarizedExperiment::colData(x)
  if ("label" %in% colnames(ch_info)) {
    ch_names <- as.character(ch_info$label)
  } else {
    ch_names <- paste0("Ch", seq_len(n_channels))
  }

  n_events <- nrow(event_df)

  # Pre-allocate result vectors
  total_rows <- n_events * n_channels
  res_event_index <- integer(total_rows)
  res_event_onset <- numeric(total_rows)
  res_channel <- character(total_rows)
  res_scr_present <- logical(total_rows)
  res_scr_amplitude <- rep(NA_real_, total_rows)
  res_scr_latency <- rep(NA_real_, total_rows)
  res_scr_rise_time <- rep(NA_real_, total_rows)
  res_scr_recovery_time <- rep(NA_real_, total_rows)

  row_idx <- 0L

  for (ei in seq_len(n_events)) {
    ev_onset <- as.numeric(event_df$onset[ei])

    for (ch in seq_len(n_channels)) {
      row_idx <- row_idx + 1L
      res_event_index[row_idx] <- ei
      res_event_onset[row_idx] <- ev_onset
      res_channel[row_idx] <- ch_names[ch]

      signal <- data[, ch]

      # Define search window in seconds
      win_start_sec <- ev_onset + onset_window[1]
      win_end_sec <- ev_onset + onset_window[2] + peak_window[2]

      # Convert to sample indices
      s_start <- max(1L, round(win_start_sec * sr) + 1L)
      s_end <- min(n_time, round(win_end_sec * sr) + 1L)

      if (s_start >= s_end || s_start > n_time) {
        res_scr_present[row_idx] <- FALSE
        next
      }

      # Find peak in the search window
      win_signal <- signal[s_start:s_end]
      peak_local <- which.max(win_signal)
      peak_sample <- s_start + peak_local - 1L
      peak_value <- signal[peak_sample]

      # Walk backward from peak to find SCR onset (local minimum)
      onset_sample <- peak_sample
      for (s in seq(peak_sample - 1L, s_start, by = -1L)) {
        if (s < 1L) break
        if (signal[s] <= signal[onset_sample]) {
          onset_sample <- s
        } else {
          break
        }
      }
      onset_value <- signal[onset_sample]

      # Compute features
      amplitude <- peak_value - onset_value
      onset_sec <- (onset_sample - 1L) / sr
      latency <- onset_sec - ev_onset
      rise_time <- (peak_sample - onset_sample) / sr

      # Validate SCR
      valid <- amplitude >= amplitude_min &&
        latency >= onset_window[1] &&
        latency <= onset_window[2]

      res_scr_present[row_idx] <- valid

      if (valid) {
        res_scr_amplitude[row_idx] <- amplitude
        res_scr_latency[row_idx] <- latency
        res_scr_rise_time[row_idx] <- rise_time

        # Compute 50% recovery time
        half_recovery_level <- peak_value - amplitude * 0.5
        recovery_time <- NA_real_
        if (peak_sample < n_time) {
          for (s in seq(peak_sample + 1L, n_time)) {
            if (signal[s] <= half_recovery_level) {
              recovery_time <- (s - peak_sample) / sr
              break
            }
          }
        }
        res_scr_recovery_time[row_idx] <- recovery_time
      }
    }
  }

  data.frame(
    event_index = res_event_index,
    event_onset = res_event_onset,
    channel = res_channel,
    scr_present = res_scr_present,
    scr_amplitude = res_scr_amplitude,
    scr_latency = res_scr_latency,
    scr_rise_time = res_scr_rise_time,
    scr_recovery_time = res_scr_recovery_time,
    stringsAsFactors = FALSE
  )
}
