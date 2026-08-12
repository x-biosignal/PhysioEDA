# Analytic peak of the Bateman SCRF

Returns the closed-form time and amplitude of the maximum of the raw
biexponential (Bateman) response used by
[`scrf`](https://x-biosignal.github.io/PhysioEDA/reference/scrf.md).

## Usage

``` r
scrfPeak(tau1 = 0.75, tau2 = 2)
```

## Arguments

- tau1:

  SCR rise time constant in seconds (default: 0.75).

- tau2:

  SCR decay time constant in seconds (default: 2.0).

## Value

A list with `time` (seconds) and `amplitude` (the raw Bateman value at
the peak).

## See also

[`scrf`](https://x-biosignal.github.io/PhysioEDA/reference/scrf.md)

## Examples

``` r
scrfPeak()
#> $time
#> [1] 1.176995
#> 
#> $amplitude
#> [1] 0.3469755
#> 
```
