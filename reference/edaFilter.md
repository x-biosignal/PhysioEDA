# Filter EDA Signal

Applies a frequency-domain (FFT-based) filter to an EDA signal. Supports
lowpass, highpass, and bandpass filter types using a smooth quadratic
transition in the frequency domain.

## Usage

``` r
edaFilter(
  x,
  type = c("lowpass", "highpass", "bandpass"),
  cutoff,
  order = 2,
  assay_name = NULL,
  output_assay = "filtered"
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- type:

  Filter type: `"lowpass"`, `"highpass"`, or `"bandpass"`.

- cutoff:

  Cutoff frequency in Hz. A single numeric value for lowpass or
  highpass, or a numeric vector of length 2 (`c(low, high)`) for
  bandpass.

- order:

  Filter steepness parameter controlling the smoothness of the frequency
  transition (default: 2). Higher values produce sharper rolloff.

- assay_name:

  Name of the input assay. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name for the output assay (default: "filtered").

## Value

A modified
[`PhysioExperiment`](https://rdrr.io/pkg/PhysioCore/man/PhysioExperiment.html)
with filtered data stored in the `output_assay` and filter parameters
recorded in `metadata(x)$eda_filter` (a list with `type`, `cutoff`,
`order`, `assay_name`, and `output_assay`).

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

## See also

[`edaDownsample`](https://x-biosignal.github.io/PhysioEDA/reference/edaDownsample.md)
for decimation with anti-aliasing,
[`edaArtifact`](https://x-biosignal.github.io/PhysioEDA/reference/edaArtifact.md)
for artifact detection and correction,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition
