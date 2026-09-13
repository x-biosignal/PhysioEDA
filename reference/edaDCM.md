# Dynamic Causal Model for Skin-Conductance Responses

Estimates trial-by-trial latent sudomotor input amplitudes that best
explain the skin-conductance signal, in the spirit of the PsPM dynamic
causal model for SCR (Bach et al., 2010, 2011). Each event contributes a
non-negative latent input whose skin-conductance response is the
canonical
[`scrf`](https://x-biosignal.github.io/PhysioEDA/reference/scrf.md)
kernel; the per-event amplitudes are estimated by non-negative least
squares against a low-frequency cosine drift basis, and the input
timings are optionally refined by a nonlinear local search (the
"dynamic" step).

## Usage

``` r
edaDCM(
  x,
  events,
  channel = 1,
  assay_name = NULL,
  tau1 = 0.75,
  tau2 = 2,
  drift = 6,
  refine_latency = TRUE,
  latency_window = 1.5
)
```

## Arguments

- x:

  A PhysioExperiment object with EDA data.

- events:

  A `data.frame` with numeric `onset` (seconds); a `condition` column is
  optional and carried through.

- channel:

  EDA channel index or label (default: 1).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- tau1, tau2:

  SCR rise/decay time constants in seconds.

- drift:

  Number of cosine drift terms (default: 6).

- refine_latency:

  Logical; refine each event's onset by a local search to improve the
  fit (default: `TRUE`).

- latency_window:

  Half-width in seconds of the latency search (default: 1.5).

## Value

A list of class `"eda_dcm"` with:

- estimates:

  a `data.frame` with `event`, `onset`, `onset_refined`, `amplitude`
  (per-event latent input), and `condition` if present.

- fitted:

  fitted signal; `residuals`; `r_squared`.

## Details

Unlike
[`edaGLM`](https://x-biosignal.github.io/PhysioEDA/reference/edaGLM.md),
which estimates one amplitude per condition, `edaDCM` estimates one
latent input amplitude per event.

## References

Bach, D.R., Daunizeau, J., Friston, K.J., & Dolan, R.J. (2010). "Dynamic
causal modelling of anticipatory skin conductance responses."
*Biological Psychology*, 85(1), 163-170.
[doi:10.1016/j.biopsycho.2010.06.007](https://doi.org/10.1016/j.biopsycho.2010.06.007)

## See also

[`edaGLM`](https://x-biosignal.github.io/PhysioEDA/reference/edaGLM.md)
for the condition-level GLM,
[`scrf`](https://x-biosignal.github.io/PhysioEDA/reference/scrf.md) for
the response kernel.

## Examples

``` r
ev <- data.frame(onset = c(10, 30, 50, 70))
pe <- edaSimulate(n_time = 1000, sr = 10, scr_count = 4, seed = 1)
fit <- edaDCM(pe, ev, refine_latency = FALSE)
fit$estimates
#>   event onset onset_refined     amplitude
#> 1     1    10            10  2.864800e-01
#> 2     2    30            30  0.000000e+00
#> 3     3    50            50 -4.356157e-18
#> 4     4    70            70  8.186079e-03
```
