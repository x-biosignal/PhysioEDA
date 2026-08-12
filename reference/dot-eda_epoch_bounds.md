# Non-overlapping Epoch Boundaries

Non-overlapping Epoch Boundaries

## Usage

``` r
.eda_epoch_bounds(n, sr, epoch_sec = 5)
```

## Arguments

- n:

  Number of samples.

- sr:

  Sampling rate in Hz.

- epoch_sec:

  Epoch length in seconds.

## Value

A `data.frame` with columns `epoch`, `start`, `end` (sample indices),
`start_sec`, `end_sec`.
