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

## Examples

``` r
PhysioEDA:::.eda_epoch_bounds(n = 100, sr = 10, epoch_sec = 2)
#>   epoch start end start_sec end_sec
#> 1     1     1  20         0     1.9
#> 2     2    21  40         2     3.9
#> 3     3    41  60         4     5.9
#> 4     4    61  80         6     7.9
#> 5     5    81 100         8     9.9
```
