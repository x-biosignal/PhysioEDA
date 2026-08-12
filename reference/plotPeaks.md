# Plot EDA Signal with SCR Peaks

Plots an EDA signal with detected SCR peaks, onsets, and amplitude lines
annotated. If no peaks data.frame is supplied, peaks are computed
automatically via
[`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md).

## Usage

``` r
plotPeaks(
  x,
  peaks = NULL,
  channel = 1,
  time_range = NULL,
  assay_name = NULL,
  main = "SCR Peaks",
  ...
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- peaks:

  A data.frame as returned by
  [`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md),
  or NULL to compute peaks automatically.

- channel:

  Integer or character specifying which channel to plot (default: 1).

- time_range:

  Numeric vector of length 2 giving start and end times in seconds. If
  NULL, the full signal is plotted.

- assay_name:

  Name of the assay to plot. If NULL, uses "phasic" if available,
  otherwise `defaultAssay(x)`.

- main:

  Character string for the plot title (default: "SCR Peaks").

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisible NULL. Called for its side effect of producing a plot.

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

## See also

[`plotEda`](https://x-biosignal.github.io/PhysioEDA/reference/plotEda.md)
for basic EDA time series plotting,
[`plotDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/plotDecompose.md)
for decomposition visualization,
[`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md)
for SCR peak detection
