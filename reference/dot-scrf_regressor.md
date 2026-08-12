# Build an Event Regressor

Places the basis `kern` at each onset and sums (a stick function
convolved with the SCRF).

## Usage

``` r
.scrf_regressor(onsets, n, sr, kern)
```

## Arguments

- onsets:

  Numeric event onsets in seconds.

- n:

  Number of samples.

- sr:

  Sampling rate in Hz.

- kern:

  Basis kernel from
  [`.scrf_kernel`](https://x-biosignal.github.io/PhysioEDA/reference/dot-scrf_kernel.md).

## Value

Numeric regressor vector of length `n`.
