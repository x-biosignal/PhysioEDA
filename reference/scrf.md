# Skin-Conductance Response Function (SCRF)

First-class canonical skin-conductance response kernel used throughout
PhysioEDA as the single source of the SCR shape (simulation, CDA /
cvxEDA deconvolution, and GLM designs). The response is the
biexponential (Bateman) function
`h(t) = exp(-t / tau2) - exp(-t / tau1)` for `t >= 0` (0 otherwise),
optionally augmented with a Gaussian bump for the PsPM-style canonical
form (Bach et al., 2010). Both the kernel and its first/second
time-derivatives are available (the derivatives drive GLM
temporal/dispersion regressors).

## Usage

``` r
scrf(
  t,
  tau1 = 0.75,
  tau2 = 2,
  form = c("bateman", "pspm_canonical"),
  normalize = c("peak", "area", "none"),
  deriv = 0,
  bump_amp = 0.3,
  bump_center = NULL,
  bump_sd = 0.7
)
```

## Arguments

- t:

  Numeric vector of times in seconds. Values `< 0` return 0.

- tau1:

  SCR rise time constant in seconds (default: 0.75).

- tau2:

  SCR decay time constant in seconds (default: 2.0). Must differ from
  `tau1`; the canonical shape has `tau2 > tau1`.

- form:

  Kernel family: `"bateman"` (pure biexponential) or `"pspm_canonical"`
  (biexponential plus a Gaussian bump).

- normalize:

  Normalisation of the kernel: `"peak"` (unit maximum, the default),
  `"area"` (unit integral), or `"none"` (raw). Derivatives are scaled by
  the same constant.

- deriv:

  Derivative order to return: 0 (kernel), 1, or 2.

- bump_amp:

  Gaussian bump amplitude, relative to the raw Bateman peak, for
  `form = "pspm_canonical"` (default: 0.3; ignored otherwise).

- bump_center:

  Gaussian bump centre in seconds (default: the analytic Bateman peak
  time).

- bump_sd:

  Gaussian bump standard deviation in seconds (default: 0.7).

## Value

A numeric vector the same length as `t`: the requested derivative of the
(optionally normalised) response kernel.

## References

Bach, D.R., Flandin, G., Friston, K.J., & Dolan, R.J. (2010). "Modelling
event-related skin conductance responses." *International Journal of
Psychophysiology*, 75(3), 349-356.
[doi:10.1016/j.ijpsycho.2010.01.005](https://doi.org/10.1016/j.ijpsycho.2010.01.005)

## See also

[`scrfPeak`](https://x-biosignal.github.io/PhysioEDA/reference/scrfPeak.md)
for the analytic Bateman maximum,
[`edaSimulate`](https://x-biosignal.github.io/PhysioEDA/reference/edaSimulate.md)
and
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
which consume it.

## Examples

``` r
t <- seq(0, 15, by = 0.1)
k <- scrf(t)                       # unit-peak Bateman kernel
dk <- scrf(t, deriv = 1)           # its temporal derivative
ck <- scrf(t, form = "pspm_canonical")
```
