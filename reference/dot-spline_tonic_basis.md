# Triangular (order-1) B-spline Tonic Basis for cvxEDA

Builds the triangular spline regressor matrix `B` whose columns are
overlapping tent functions spaced every `delta_knot` seconds, following
Greco et al. (2016). The tonic (SCL) component is modelled as
`B %*% l + C %*% d` where `C` is a linear trend.

## Usage

``` r
.spline_tonic_basis(n, dt, delta_knot = 10)
```

## Arguments

- n:

  Integer number of samples.

- dt:

  Sampling interval in seconds.

- delta_knot:

  Knot spacing in seconds (default: 10).

## Value

An `n x nB` matrix of spline regressors.
