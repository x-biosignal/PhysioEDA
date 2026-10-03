# Running Mean Smoother

Computes a symmetric running mean with edge handling via partial
windows.

## Usage

``` r
.running_mean(x, window)
```

## Arguments

- x:

  Numeric vector.

- window:

  Integer window size (must be odd).

## Value

Smoothed numeric vector of same length.

## Examples

``` r
PhysioEDA:::.running_mean(c(1, 2, 3, 4, 5, 6), window = 3)
#> [1] 1.5 2.0 3.0 4.0 5.0 5.5
```
