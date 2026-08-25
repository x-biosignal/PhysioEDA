# Create a Simulated EDA PhysioExperiment

Generates a synthetic electrodermal activity (EDA) PhysioExperiment
object with known tonic (SCL) and phasic (SCR) components. Useful for
testing, demonstrations, and vignettes.

## Usage

``` r
make_eda(n_time = 6000, n_channels = 1, sr = 10)
```

## Arguments

- n_time:

  Number of time points (default: 6000).

- n_channels:

  Number of EDA channels (default: 1).

- sr:

  Sampling rate in Hz (default: 10).

## Value

A
[`PhysioExperiment`](https://x-biosignal.github.io/PhysioCore//reference/PhysioExperiment.html)
object with a single `"raw"` assay containing simulated EDA data.
Channel metadata has `type = "EDA"` and `unit = "uS"`.

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

## See also

[`make_eda_with_scr`](https://x-biosignal.github.io/PhysioEDA/reference/make_eda_with_scr.md)
for EDA with stimulus events,
[`edaSimulate`](https://x-biosignal.github.io/PhysioEDA/reference/edaSimulate.md)
for the underlying simulation function,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition

## Examples

``` r
x <- make_eda()
x
#> class: PhysioExperiment
#> dim: 6000 x 1 
#> assays(1): raw
#> samplingRate: 10 Hz
#> channels(1): EDA1
#> colData names(3): label, type, unit
```
