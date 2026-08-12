# Bateman ARMA Matrices for cvxEDA

Builds the banded auto-regressive (`A`) and moving-average (`M`)
matrices of the discretised biexponential (Bateman) system used by the
exact cvxEDA formulation (Greco et al., 2016). The latent variable `q`
relates to the SMNA driver by `p = A %*% q` and to the phasic component
by `r = M %*% q`.

## Usage

``` r
.bateman_arma(n, dt, tau_rise, tau_decay)
```

## Arguments

- n:

  Integer number of samples.

- dt:

  Sampling interval in seconds (`1 / sr`).

- tau_rise:

  SCR rise time constant in seconds.

- tau_decay:

  SCR decay time constant in seconds.

## Value

A list with the `n x n` matrices `A` and `M`.
