# Time-varying electrodermal sympathetic index (TVSymp)

A time-resolved sympathetic-tone index built by complex demodulation of
the dominant components of the electrodermal signal in the sympathetic
band (Posada-Quintero et al. 2016). The two strongest spectral
components in the band are demodulated to their instantaneous
amplitudes, combined, and normalised to a reproducible trajectory that
rises when sympathetic oscillations appear.

## Usage

``` r
edaTVSymp(
  x,
  band = c(0.08, 0.24),
  assay_name = NULL,
  resample_hz = 2,
  n_components = 2,
  lp_hz = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with an electrodermal assay.

- band:

  Numeric `c(low, high)` sympathetic band in Hz (default
  `c(0.08, 0.24)`).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- resample_hz:

  Rate (Hz) the signal is downsampled to (default 2; `NULL` keeps the
  native rate).

- n_components:

  Number of dominant band components to demodulate (default 2).

- lp_hz:

  Low-pass cutoff (Hz) of the demodulation envelope (default 0.5 of the
  band width).

## Value

A data frame with one row per time sample per channel: `channel`,
`time_sec`, `tvsymp` (combined instantaneous amplitude) and
`tvsymp_norm` (unit-variance normalised). The dominant component
frequencies are attached as the `"components"` attribute.

## References

Posada-Quintero, H.F. et al. (2016). Time-varying analysis of
electrodermal activity during exercise. PLoS ONE 11(8):e0161553.

## See also

[`edaSympatheticIndex()`](https://x-biosignal.github.io/PhysioEDA/reference/edaSympatheticIndex.md)

## Examples

``` r
pe <- edaSimulate(n_time = 1200, sr = 4, scr_count = 6, seed = 1)
tv <- edaTVSymp(pe)
head(tv)
#>   channel time_sec    tvsymp tvsymp_norm
#> 1       1      0.0 0.1596698   1.2560419
#> 2       1      0.5 0.1565969   1.2318688
#> 3       1      1.0 0.1480711   1.1648007
#> 4       1      1.5 0.1362271   1.0716296
#> 5       1      2.0 0.1244556   0.9790289
#> 6       1      2.5 0.1166002   0.9172346
```
