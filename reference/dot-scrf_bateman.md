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
