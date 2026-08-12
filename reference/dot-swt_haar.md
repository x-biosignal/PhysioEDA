# Undecimated (Stationary) Haar Wavelet Transform

A dependency-free maximal-overlap ("a trous") Haar wavelet transform.
Every returned component has the same length as `x` (undecimated), so
the detail coefficients are time-aligned with the samples - ideal for
locating high-frequency artifacts. The normalisation matches
`waveslim::modwt(x, "haar")`: the sum of squared coefficients over all
levels plus the final smooth equals the energy of `x`.

## Usage

``` r
.swt_haar(x, J = 4)
```

## Arguments

- x:

  Numeric signal vector.

- J:

  Number of decomposition levels (default: 4).

## Value

A named list `d1..dJ` (level detail coefficients) and `sJ` (level-J
smooth), each the same length as `x`.
