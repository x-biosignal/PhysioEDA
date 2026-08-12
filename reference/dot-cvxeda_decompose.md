# cvxEDA Decomposition for a Single Channel (Simplified)

Performs a simplified convex optimization-based EDA decomposition using
an iterative ADMM-like approach with FFT-based
convolution/deconvolution. Based on Greco et al. (2016), but avoids the
CVXR dependency by using an iterative proximal gradient method with
Wiener deconvolution.

## Usage

``` r
.cvxeda_decompose(
  signal,
  sr,
  tau1 = 0.75,
  tau2 = 2,
  alpha = 0.01,
  gamma = 0.1,
  max_iter = 50,
  tol = 1e-04
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

  L1 sparsity penalty on the driver signal (default: 0.01).

- gamma:

  Smoothness weight for tonic component (default: 0.1).

- max_iter:

  Maximum number of ADMM iterations (default: 50).

- tol:

  Convergence tolerance for the driver signal (default: 1e-4).

## Value

A list with components:

- tonic:

  Numeric vector of tonic (SCL) component.

- phasic:

  Numeric vector of phasic (SCR) component.

- driver:

  Numeric vector of the sparse driver signal.

- alpha:

  The alpha parameter used.

- gamma:

  The gamma parameter used.

- iterations:

  Number of iterations performed.

- converged:

  Logical indicating whether the algorithm converged.

## References

Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
electrodermal activity processing." *IEEE Transactions on Biomedical
Engineering*, 63(4), 797-804.
[doi:10.1109/TBME.2015.2474131](https://doi.org/10.1109/TBME.2015.2474131)

## See also

[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for the user-facing decomposition wrapper,
[`.cda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cda_decompose.md)
for the CDA alternative
