# Raw biexponential (Bateman) kernel and its analytic time-derivatives

Raw biexponential (Bateman) kernel and its analytic time-derivatives

## Usage

``` r
.scrf_bateman(x, tau1, tau2, deriv = 0)
```

## Arguments

- x:

  Numeric time vector (seconds).

- tau1:

  Rise time constant in seconds.

- tau2:

  Decay time constant in seconds.

- deriv:

  Derivative order: 0 (kernel), 1, or 2.

## Value

Numeric vector of the (`deriv`-th derivative of the) raw kernel. Values
are not zeroed for `x < 0`; the caller does that.

## Examples

``` r
PhysioEDA:::.scrf_bateman(seq(0, 10, by = 0.5), tau1 = 0.75, tau2 = 2.0)
#>  [1] 0.000000000 0.265383664 0.342933522 0.337031270 0.298395990 0.250830804
#>  [7] 0.204814521 0.164370381 0.130507333 0.102920472 0.080812365 0.063274469
#> [13] 0.049451606 0.038601976 0.030108956 0.023472346 0.018292330 0.014252267
#> [19] 0.011102852 0.008648541 0.006736327
```
