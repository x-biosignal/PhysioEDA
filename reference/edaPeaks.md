# Detect Skin Conductance Response (SCR) Peaks

Identifies SCR peaks in an EDA signal using either a gradient-based
zero-crossing method or an amplitude threshold method. Returns onset,
peak, amplitude, rise time, and recovery time for each detected SCR.

## Usage

``` r
edaPeaks(
  x,
  method = c("gradient", "threshold"),
  amplitude_min = 0.01,
  rise_time_min = 0.1,
  rise_time_max = 5,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- method:

  Detection method: `"gradient"` (first-derivative zero-crossing) or
  `"threshold"` (amplitude threshold). Default is `"gradient"`.

- amplitude_min:

  Minimum SCR amplitude in microsiemens (default: 0.01).

- rise_time_min:

  Minimum rise time in seconds (default: 0.1).

- rise_time_max:

  Maximum rise time in seconds (default: 5.0).

- assay_name:

  Name of the input assay. If NULL, uses `"phasic"` if available,
  otherwise `defaultAssay(x)`.

## Value

A `data.frame` with one row per detected SCR and the following columns:

- channel:

  Character channel label.

- onset_sample:

  Integer sample index of SCR onset.

- onset_sec:

  Numeric onset time in seconds.

- peak_sample:

  Integer sample index of SCR peak.

- peak_sec:

  Numeric peak time in seconds.

- amplitude:

  Numeric SCR amplitude in microsiemens (peak minus onset).

- rise_time:

  Numeric rise time from onset to peak in seconds.

- recovery_time:

  Numeric 50% recovery time in seconds, or `NA`.

Returns an empty `data.frame` with the same columns if no peaks are
found.

## References

Bach, D.R., et al. (2010). "Modelling event-related skin conductance
responses." *International Journal of Psychophysiology*, 75(3), 349-356.
[doi:10.1016/j.ijpsycho.2010.01.005](https://doi.org/10.1016/j.ijpsycho.2010.01.005)

Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
electrodermal activity." *Journal of Neuroscience Methods*, 190(1),
80-91.
[doi:10.1016/j.jneumeth.2010.04.028](https://doi.org/10.1016/j.jneumeth.2010.04.028)

## See also

[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition (run first),
[`edaFeatures`](https://x-biosignal.github.io/PhysioEDA/reference/edaFeatures.md)
for summary feature extraction,
[`plotPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/plotPeaks.md)
for peak visualization,
[`edaErscr`](https://x-biosignal.github.io/PhysioEDA/reference/edaErscr.md)
for event-related SCR analysis
