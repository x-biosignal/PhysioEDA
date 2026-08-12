# Electrodermal sympathetic index (EDASymp)

Spectral power of the electrodermal signal in the sympathetic band
(0.045-0.25 Hz by default), following Posada-Quintero et al. (2016). The
signal is optionally downsampled, its Welch power spectral density is
estimated, and the band power is reported both in absolute units and
normalised by the total spectral power.

## Usage

``` r
edaSympatheticIndex(
  x,
  band = c(0.045, 0.25),
  assay_name = NULL,
  resample_hz = 2,
  seg_sec = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with an electrodermal assay.

- band:

  Numeric `c(low, high)` sympathetic band in Hz (default
  `c(0.045, 0.25)`).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- resample_hz:

  Rate (Hz) the signal is downsampled to before spectral estimation
  (default 2; `NULL` keeps the native rate). Ignored if it is not below
  the native sampling rate.

- seg_sec:

  Welch segment length in seconds (`NULL` = one whole-signal Hann
  window).

## Value

A data frame with one row per channel: `channel`, `edasymp` (absolute
band power), `edasymp_norm` (band power / total power, in `[0, 1]`, or
`NA` for a constant / zero-power channel), `total_power`, and the
`band_low`/`band_high` used.

## References

Posada-Quintero, H.F. et al. (2016). Power spectral density analysis of
electrodermal activity for sympathetic function assessment. Ann Biomed
Eng 44(10):3124-3135.

## See also

[`edaTVSymp()`](https://x-biosignal.github.io/PhysioEDA/reference/edaTVSymp.md)

## Examples

``` r
pe <- edaSimulate(n_time = 1200, sr = 4, scr_count = 6, seed = 1)
edaSympatheticIndex(pe)
#>   channel     edasymp edasymp_norm total_power band_low band_high
#> 1       1 0.006002253    0.1195644  0.05020099    0.045      0.25
```
