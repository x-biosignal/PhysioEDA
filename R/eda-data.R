#' Create a Simulated EDA PhysioExperiment
#'
#' Generates a synthetic electrodermal activity (EDA) PhysioExperiment object
#' with known tonic (SCL) and phasic (SCR) components. Useful for testing,
#' demonstrations, and vignettes.
#'
#' @param n_time Number of time points (default: 6000).
#' @param n_channels Number of EDA channels (default: 1).
#' @param sr Sampling rate in Hz (default: 10).
#' @return A \code{\link[PhysioCore]{PhysioExperiment}} object with a single
#'   \code{"raw"} assay containing simulated EDA data. Channel metadata has
#'   \code{type = "EDA"} and \code{unit = "uS"}.
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{make_eda_with_scr}} for EDA with stimulus events,
#'   \code{\link{edaSimulate}} for the underlying simulation function,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition
#' @export
#' @examples
#' x <- make_eda()
#' x
make_eda <- function(n_time = 6000, n_channels = 1, sr = 10) {
  edaSimulate(n_time = n_time, n_channels = n_channels, sr = sr,
              scr_count = 5, scl_level = 5.0, scr_amplitude = 0.5,
              noise_sd = 0.01, seed = 42)
}

#' Create a Simulated EDA PhysioExperiment with SCR Events
#'
#' Generates a synthetic electrodermal activity (EDA) PhysioExperiment object
#' with known SCR events time-locked to stimulus markers. The object includes
#' regularly spaced stimulus events stored via
#' \code{\link[PhysioCore]{PhysioEvents}}, making it suitable for testing
#' event-related SCR analysis workflows.
#'
#' @param n_time Number of time points (default: 6000).
#' @param n_channels Number of EDA channels (default: 1).
#' @param sr Sampling rate in Hz (default: 10).
#' @param n_events Number of stimulus events to embed (default: 4).
#' @param event_interval Interval between events in seconds (default: 10).
#'   Note: events are spaced evenly within the usable signal duration rather
#'   than at exact fixed intervals.
#' @return A \code{\link[PhysioCore]{PhysioExperiment}} object with a single
#'   \code{"raw"} assay containing simulated EDA data and stimulus events
#'   accessible via \code{\link[PhysioCore]{getEvents}}. Each event has
#'   \code{type = "stimulus"} and \code{duration = 0.5} seconds.
#' @references
#' Bach, D.R., et al. (2010). "Modelling event-related skin conductance
#' responses." \emph{International Journal of Psychophysiology}, 75(3),
#' 349-356. \doi{10.1016/j.ijpsycho.2010.01.005}
#' @seealso \code{\link{make_eda}} for basic EDA without events,
#'   \code{\link{edaErscr}} for event-related SCR analysis,
#'   \code{\link{edaSimulate}} for the underlying simulation function
#' @export
#' @examples
#' x <- make_eda_with_scr()
#' x
make_eda_with_scr <- function(n_time = 6000, n_channels = 1, sr = 10,
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
