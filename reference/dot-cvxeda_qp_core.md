# Exact cvxEDA Quadratic-Program Core (Single Channel, Fixed Rate)

Solves the Greco et al. (2016) cvxEDA convex problem for one
standardised signal at a single sampling rate. The signal is decomposed
as `y = M q + B l + C d + noise` with a sparse non-negative driver
`p = A q >= 0`, minimising
`0.5 * sum((M q + B l + C d - y)^2) + alpha * sum(A q) + 0.5 * gamma * sum(l^2)`
subject to `A q >= 0`. This is a quadratic program solved with an
active-set solver
([`quadprog::solve.QP`](https://rdrr.io/pkg/quadprog/man/solve.QP.html));
the returned optimum is exact.

## Usage

``` r
.cvxeda_qp_core(
  y,
  sr,
  tau1 = 0.75,
  tau2 = 2,
  alpha = 0.01,
  gamma = 0.1,
  delta_knot = 10
)
```

## Arguments

- y:

  Numeric EDA signal (single channel).

- sr:

  Sampling rate in Hz.

- tau1:

  SCR rise time constant in seconds (default: 0.75).

- tau2:

  SCR decay time constant in seconds (default: 2.0).

- alpha:

  L1 sparsity penalty on the driver (default: 0.01).

- gamma:

  L2 smoothness weight on the spline tonic (default: 0.1).

- delta_knot:

  Tonic spline knot spacing in seconds (default: 10).

## Value

A list with `tonic`, `phasic`, `driver` (each of length `length(y)`),
the final QP `objective`, the objective at the trivial (driver-free)
initialisation `obj_init`, and `converged`.

## References

Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
electrodermal activity processing." *IEEE Transactions on Biomedical
Engineering*, 63(4), 797-804.
[doi:10.1109/TBME.2015.2474131](https://doi.org/10.1109/TBME.2015.2474131)
