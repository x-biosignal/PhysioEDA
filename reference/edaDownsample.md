# Downsample EDA Signal

Reduces the sampling rate of an EDA signal by first applying an
anti-aliasing lowpass filter at half the target sampling rate, then
decimating the signal.

## Usage

``` r
edaDownsample(x, target_sr, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- target_sr:

  Target sampling rate in Hz. Must be lower than the current sampling
  rate.

- assay_name:

  Name of the input assay. If NULL, uses `defaultAssay(x)`.

## Value

A new
[`PhysioExperiment`](https://x-biosignal.github.io/PhysioCore//reference/PhysioExperiment.html)
with the downsampled signal in the `"raw"` assay and the sampling rate
set to the actual achieved rate. Events are preserved. Downsampling
parameters are stored in `metadata(x)$eda_downsample` (a list with
`original_sr`, `target_sr`, `factor`, `original_n_time`, and
`new_n_time`).

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

## See also

[`edaFilter`](https://x-biosignal.github.io/PhysioEDA/reference/edaFilter.md)
for frequency-domain filtering,
[`edaArtifact`](https://x-biosignal.github.io/PhysioEDA/reference/edaArtifact.md)
for artifact detection and correction
