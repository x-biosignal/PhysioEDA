# Event-Related SCR Analysis

Performs event-related skin conductance response (ER-SCR) analysis by
extracting SCR features time-locked to experimental events. For each
event, the function searches for an SCR onset within a specified window
and, if found, computes amplitude, latency, rise time, and 50 percent
recovery time.

## Usage

``` r
edaErscr(
  x,
  event_type = NULL,
  onset_window = c(1, 4),
  peak_window = c(0.5, 5),
  amplitude_min = 0.01,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data and events.

- event_type:

  Character string specifying the event type to analyze. If NULL
  (default), all events are used.

- onset_window:

  Numeric vector of length 2 giving the valid SCR onset latency window
  in seconds relative to event onset (default: c(1, 4)).

- peak_window:

  Numeric vector of length 2 giving the minimum and maximum time from
  SCR onset to peak in seconds (default: c(0.5, 5)).

- amplitude_min:

  Minimum SCR amplitude in microsiemens to be considered a valid
  response (default: 0.01).

- assay_name:

  Name of the assay to use. If NULL, uses "phasic" if available,
  otherwise the default assay.

## Value

A data.frame with one row per event per channel containing columns:

- event_index:

  Integer index of the event.

- event_onset:

  Numeric onset time of the event in seconds.

- channel:

  Character name of the channel.

- scr_present:

  Logical indicating whether a valid SCR was detected.

- scr_amplitude:

  Numeric SCR amplitude (peak minus onset), or NA.

- scr_latency:

  Numeric latency from event to SCR onset in seconds, or NA.

- scr_rise_time:

  Numeric time from SCR onset to peak in seconds, or NA.

- scr_recovery_time:

  Numeric time for 50 percent recovery from peak, or NA.

## References

Bach, D.R., et al. (2010). "Modelling event-related skin conductance
responses." *International Journal of Psychophysiology*, 75(3), 349-356.
[doi:10.1016/j.ijpsycho.2010.01.005](https://doi.org/10.1016/j.ijpsycho.2010.01.005)

## See also

[`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md)
for general SCR peak detection,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition (run first),
[`edaFeatures`](https://x-biosignal.github.io/PhysioEDA/reference/edaFeatures.md)
for summary feature extraction,
[`plotPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/plotPeaks.md)
for visualizing SCR peaks
