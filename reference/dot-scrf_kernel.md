# Sampled SCRF Basis Kernel

Samples the skin-conductance response function
([`scrf`](https://x-biosignal.github.io/PhysioEDA/reference/scrf.md)) on
the signal grid for use as a GLM/DCM basis regressor kernel.

## Usage

``` r
.scrf_kernel(sr, tau1 = 0.75, tau2 = 2, deriv = 0, len_sec = NULL)
```

## Arguments

- sr:

  Sampling rate in Hz.

- tau1:

  SCR rise time constant in seconds.

- tau2:

  SCR decay time constant in seconds.

- deriv:

  Derivative order (0 canonical, 1 temporal derivative).

- len_sec:

  Kernel length in seconds (default: `10 * tau2`).

## Value

Numeric kernel vector, scaled to unit peak absolute value.
