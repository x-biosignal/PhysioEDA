# DDA Decomposition for a Single EDA Channel

Performs Discrete Decomposition Analysis (DDA) on a single EDA channel
(Benedek & Kaernbach, 2010). Unlike Continuous Decomposition Analysis
(CDA), which yields a continuous driver, DDA represents the phasic
activity as a small set of *discrete* SCRs: the continuous CDA driver is
reduced to a sparse impulse train (one impulse per significant SCR,
placed at its onset and weighted by the integrated driver of that SCR),
which is reconvolved with the Bateman response to obtain the phasic
component.

## Usage

``` r
.dda_decompose(
  signal,
  sr,
  tau1 = 0.75,
  tau2 = 2,
  optimize_tau = FALSE,
  amp_threshold = 0.01,
  min_interval_sec = 1
)
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

  Logical; passed to
  [`.cda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cda_decompose.md)
  to optimise `tau2` before discretisation (default: `FALSE`).

- amp_threshold:

  Minimum driver amplitude for a significant SCR (default: 0.01).

- min_interval_sec:

  Minimum separation in seconds between distinct SCRs (default: 1).

## Value

A list with components:

- tonic:

  Numeric vector of the tonic (SCL) component.

- phasic:

  Numeric vector of the phasic (SCR) component from the discrete driver.

- driver:

  Numeric vector of the discrete (impulse-train) driver.

- tau1:

  Rise time constant used.

- tau2:

  Decay time constant used.

- onsets:

  Onset times of the discrete SCRs (seconds).

- nSCR:

  Number of discrete SCRs.

## References

Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
electrodermal activity." *Journal of Neuroscience Methods*, 190(1),
80-91.
[doi:10.1016/j.jneumeth.2010.04.028](https://doi.org/10.1016/j.jneumeth.2010.04.028)

## See also

[`.cda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cda_decompose.md)
for the continuous (CDA) variant,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for the user-facing wrapper.
