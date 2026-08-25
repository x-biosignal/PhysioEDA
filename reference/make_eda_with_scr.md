# Create a Simulated EDA PhysioExperiment with SCR Events

Generates a synthetic electrodermal activity (EDA) PhysioExperiment
object with known SCR events time-locked to stimulus markers. The object
includes regularly spaced stimulus events stored via
[`PhysioEvents`](https://x-biosignal.github.io/PhysioCore//reference/PhysioEvents.html),
making it suitable for testing event-related SCR analysis workflows.

## Usage

``` r
make_eda_with_scr(
  n_time = 6000,
  n_channels = 1,
  sr = 10,
  n_events = 4,
  event_interval = 10
)
```

## Arguments

- n_time:

  Number of time points (default: 6000).

- n_channels:

  Number of EDA channels (default: 1).

- sr:

  Sampling rate in Hz (default: 10).

- n_events:

  Number of stimulus events to embed (default: 4).

- event_interval:

  Interval between events in seconds (default: 10). Note: events are
  spaced evenly within the usable signal duration rather than at exact
  fixed intervals.

## Value

A
[`PhysioExperiment`](https://x-biosignal.github.io/PhysioCore//reference/PhysioExperiment.html)
object with a single `"raw"` assay containing simulated EDA data and
stimulus events accessible via
[`getEvents`](https://x-biosignal.github.io/PhysioCore//reference/getEvents.html).
Each event has `type = "stimulus"` and `duration = 0.5` seconds.

## References

Bach, D.R., et al. (2010). "Modelling event-related skin conductance
responses." *International Journal of Psychophysiology*, 75(3), 349-356.
[doi:10.1016/j.ijpsycho.2010.01.005](https://doi.org/10.1016/j.ijpsycho.2010.01.005)

## See also

[`make_eda`](https://x-biosignal.github.io/PhysioEDA/reference/make_eda.md)
for basic EDA without events,
[`edaErscr`](https://x-biosignal.github.io/PhysioEDA/reference/edaErscr.md)
for event-related SCR analysis,
[`edaSimulate`](https://x-biosignal.github.io/PhysioEDA/reference/edaSimulate.md)
for the underlying simulation function

## Examples

``` r
x <- make_eda_with_scr()
x
#> class: PhysioExperiment
#> dim: 6000 x 1 
#> assays(1): raw
#> samplingRate: 10 Hz
#> channels(1): EDA1
#> colData names(3): label, type, unit
#> events: 4 
```
