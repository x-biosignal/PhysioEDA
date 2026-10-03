# Simulate Synthetic EDA Signals

Generates synthetic electrodermal activity (EDA) signals with known
tonic (SCL) and phasic (SCR) components for testing and demonstration.
SCRs are modeled as biexponential impulse responses (Bateman function).

## Usage

``` r
edaSimulate(
  n_time = 6000,
  n_channels = 1,
  sr = 10,
  scr_count = 5,
  scl_level = 5,
  scr_amplitude = 0.5,
  noise_sd = 0.01,
  seed = NULL
)
```

## Arguments

- n_time:

  Number of time points (default: 6000).

- n_channels:

  Number of EDA channels (default: 1).

- sr:

  Sampling rate in Hz (default: 10).

- scr_count:

  Number of SCRs to embed (default: 5).

- scl_level:

  Baseline skin conductance level in microsiemens (default: 5.0).

- scr_amplitude:

  Mean SCR amplitude in microsiemens (default: 0.5).

- noise_sd:

  Standard deviation of Gaussian noise (default: 0.01).

- seed:

  Random seed for reproducibility (default: NULL).

## Value

A
[`PhysioExperiment`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/PhysioExperiment.html)
object with a single `"raw"` assay containing the simulated EDA signal
(time x channels matrix). Channel metadata has `type = "EDA"` and
`unit = "uS"`. The sampling rate is set to `sr`. The ground-truth
components used to build the signal are stored in
`metadata(x)$eda_truth` as a list with:

- tonic:

  time x channels matrix of the true tonic (SCL) component.

- phasic:

  time x channels matrix of the true phasic (SCR) component.

- onsets:

  a `data.frame` of the true SCR onsets with columns `channel`,
  `sample`, `time_sec`, and `amplitude`.

This lets `edaSimulate()` serve as a ground-truth oracle for testing
decomposition (e.g.
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
with `method="cvxeda"`).

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
electrodermal activity." *Journal of Neuroscience Methods*, 190(1),
80-91.
[doi:10.1016/j.jneumeth.2010.04.028](https://doi.org/10.1016/j.jneumeth.2010.04.028)

## See also

[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition,
[`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md)
for SCR detection,
[`edaFilter`](https://x-biosignal.github.io/PhysioEDA/reference/edaFilter.md)
for signal filtering

## Examples

``` r
x <- edaSimulate(n_time = 600, sr = 10, seed = 1)
x
#> class: PhysioExperiment
#> dim: 600 x 1 
#> assays(1): raw
#> samplingRate: 10 Hz
#> channels(1): EDA1
#> colData names(3): label, type, unit
str(S4Vectors::metadata(x)$eda_truth, max.level = 1)
#> List of 3
#>  $ tonic : num [1:600, 1] 5.1 5.1 5.1 5.11 5.11 ...
#>  $ phasic: num [1:600, 1] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ onsets:'data.frame':  5 obs. of  4 variables:
```
