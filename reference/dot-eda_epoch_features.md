# Per-epoch Artifact Features for EDA

Computes, for each non-overlapping epoch, the derivative statistics and
the undecimated-Haar wavelet detail statistics used by the wavelet rule
and the logistic model (EDA-Explorer-style features; Taylor et al.,
2015).

## Usage

``` r
.eda_epoch_features(sig, sr, epoch_sec = 5, J = 4)
```

## Arguments

- sig:

  Numeric EDA signal (single channel).

- sr:

  Sampling rate in Hz.

- epoch_sec:

  Epoch length in seconds (default: 5).

- J:

  Wavelet levels (default: 4).

## Value

A `data.frame` with one row per epoch: epoch bounds plus `amp`,
`maxabsd`, `avgabsd`, `maxabs2d`, `avgabs2d`, `wav_max_d1`,
`wav_std_d1`, `wav_max_d2`, `wav_std_d2`, `wav_std_d3`.
