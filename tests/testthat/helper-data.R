#' Create test EDA PhysioExperiment
#' @param n_time Number of time points
#' @param n_channels Number of EDA channels
#' @param sr Sampling rate in Hz
#' @return PhysioExperiment with simulated EDA data
make_eda <- function(n_time = 6000, n_channels = 1, sr = 10) {
  edaSimulate(n_time = n_time, n_channels = n_channels, sr = sr,
              scr_count = 5, scl_level = 5.0, scr_amplitude = 0.5,
              noise_sd = 0.01, seed = 42)
}

#' Create test EDA PhysioExperiment with stimulus events
#' @param n_time Number of time points
#' @param n_channels Number of EDA channels
#' @param sr Sampling rate in Hz
#' @param n_events Number of stimulus events
#' @param event_interval Interval between events in seconds
#' @return PhysioExperiment with simulated EDA data and events
make_eda_event <- function(n_time = 6000, n_channels = 1, sr = 10,
                           n_events = 4, event_interval = 10) {
  x <- edaSimulate(n_time = n_time, n_channels = n_channels, sr = sr,
                   scr_count = n_events, scl_level = 5.0,
                   scr_amplitude = 0.8, noise_sd = 0.01, seed = 42)

  # Add stimulus events at regular intervals
  duration_sec <- (n_time - 1) / sr
  margin <- duration_sec * 0.1
  event_onsets <- seq(margin, duration_sec - margin,
                      length.out = n_events)

  x <- setEvents(x, PhysioEvents(
    onset = event_onsets,
    duration = rep(0.5, n_events),
    type = rep("stimulus", n_events),
    value = paste0("S", seq_len(n_events))
  ))

  x
}
