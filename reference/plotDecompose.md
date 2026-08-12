# Plot EDA Decomposition

Displays the tonic/phasic decomposition of an EDA signal in a
multi-panel layout. Requires that
[`edaDecompose()`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
has been run first so that "tonic" and "phasic" assays exist.

## Usage

``` r
plotDecompose(
  x,
  channel = 1,
  time_range = NULL,
  main = "EDA Decomposition",
  ...
)
```

## Arguments

- x:

  A PhysioExperiment object with tonic and phasic assays.

- channel:

  Integer or character specifying which channel to plot (default: 1).

- time_range:

  Numeric vector of length 2 giving start and end times in seconds. If
  NULL, the full signal is plotted.

- main:

  Character string for the overall plot title (default: "EDA
  Decomposition").

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
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for performing the decomposition,
[`plotPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/plotPeaks.md)
for SCR peak visualization
