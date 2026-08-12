# Analytic Bateman peak time

Closed-form time of the maximum of the biexponential (Bateman) response
`h(t) = exp(-t / tau2) - exp(-t / tau1)`, obtained from the stationarity
condition `dh/dt = 0`.

## Usage

``` r
.scrf_peak_time(tau1, tau2)
```

## Arguments

- tau1:

  Rise time constant in seconds.

- tau2:

  Decay time constant in seconds (`tau2 > tau1`).

## Value

The peak time in seconds.
