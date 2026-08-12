# Exact cvxEDA Decomposition for a Single Channel (Greco QP)

Full-signal wrapper around
[`.cvxeda_qp_core`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cvxeda_qp_core.md):
it decimates the signal to a tractable rate for the dense quadratic
program (electrodermal activity is slow, so a few Hz suffices), solves
the exact cvxEDA QP, and interpolates the tonic / phasic / driver
components back to the native grid. If quadprog is unavailable, or the
solve fails, it degrades gracefully to the Wiener approximation
([`.cvxeda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cvxeda_decompose.md)).

## Usage

``` r
.cvxeda_qp(
  signal,
  sr,
  tau1 = 0.75,
  tau2 = 2,
  alpha = 0.01,
  gamma = 0.1,
  delta_knot = 10,
  qp_sr = 4,
  max_n = 2000
)
```

## Arguments

- signal:

  Numeric vector of EDA signal values.

- sr:

  Sampling rate in Hz.

- tau1:

  SCR rise time constant in seconds (default: 0.75).

- tau2:

  SCR decay time constant in seconds (default: 2.0).

- alpha:

  L1 sparsity penalty on the driver (default: 0.01).

- gamma:

  L2 smoothness weight on the tonic (default: 0.1).

- delta_knot:

  Tonic spline knot spacing in seconds (default: 10).

- qp_sr:

  Target rate (Hz) for the QP solve; the signal is decimated towards
  this rate (default: 4).

- max_n:

  Maximum number of samples in the QP; further decimation is applied
  above this (default: 2000).

## Value

A list with `tonic`, `phasic`, and `driver` on the native grid, plus
`driver_qp` / `t_qp` / `sr_qp` at the QP rate, the QP `objective` and
`obj_init`, `converged`, and the `engine` ("qp" or "wiener_fallback")
and `solver` used.

## References

Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
electrodermal activity processing." *IEEE Transactions on Biomedical
Engineering*, 63(4), 797-804.
[doi:10.1109/TBME.2015.2474131](https://doi.org/10.1109/TBME.2015.2474131)

## See also

[`.cvxeda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cvxeda_decompose.md)
for the fast Wiener approximation,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for the user-facing wrapper.
