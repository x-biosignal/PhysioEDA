# CDA Decomposition for a Single EDA Channel

Performs Continuous Decomposition Analysis (CDA) on a single EDA channel
using the method described by Benedek & Kaernbach (2010). The signal is
deconvolved with a Bateman impulse response function to obtain a
sudomotor nerve activity (driver) signal, which is then reconvolved to
separate tonic and phasic components.

## Usage

``` r
.cda_decompose(signal, sr, tau1 = 0.75, tau2 = 2, optimize_tau = FALSE)
```

## Arguments

- signal:

  Numeric vector of EDA values.

- sr:

  Sampling rate in Hz.

- tau1:

  Rise time constant in seconds (default: 0.75).

- tau2:

  Decay time constant in seconds (default: 2.0).

- optimize_tau:

  Logical; if `TRUE`, optimise the decay constant `tau2` to minimise the
  negative energy of the deconvolved driver (a Ledalab-style criterion;
  the driver should be sparse and non-negative). `tau1` is held fixed.
  Default `FALSE`.

## Value

A list with components:

- tonic:

  Numeric vector of the tonic (SCL) component.

- phasic:

  Numeric vector of the phasic (SCR) component.

- driver:

  Numeric vector of the sudomotor nerve activity driver signal.

- tau1:

  Rise time constant used.

- tau2:

  Decay time constant used (optimised when `optimize_tau`).

## References

Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
electrodermal activity." *Journal of Neuroscience Methods*, 190(1),
80-91.
[doi:10.1016/j.jneumeth.2010.04.028](https://doi.org/10.1016/j.jneumeth.2010.04.028)

## See also

[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for the user-facing decomposition wrapper,
[`.dda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-dda_decompose.md)
for the discrete (DDA) alternative,
[`.cvxeda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cvxeda_decompose.md)
for the cvxEDA alternative
