# Low-frequency Drift (High-pass) Nuisance Basis

A discrete-cosine drift basis (intercept plus low-frequency cosines)
that absorbs the slow tonic level, equivalent to a high-pass filter in
the GLM.

## Usage

``` r
.drift_basis(n, k = 6)
```

## Arguments

- n:

  Number of samples.

- k:

  Number of cosine terms (default: 6).

## Value

An `n x (k + 1)` matrix.
