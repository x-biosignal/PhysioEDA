# PsPM-style GLM for Skin-Conductance Responses

Fits a general linear model of the skin-conductance signal against
event-related regressors, each a stick function at the event onsets
convolved with the canonical skin-conductance response function
([`scrf`](https://x-biosignal.github.io/PhysioEDA/reference/scrf.md);
Bach & Friston, 2013). Events are grouped into conditions; each
condition yields an amplitude estimate (beta) with a standard error and
t statistic. A low-frequency cosine drift basis absorbs the tonic level.
Optionally, a temporal-derivative basis is added per condition to
accommodate response latency differences.

## Usage

``` r
edaGLM(
  x,
  events,
  conditions = NULL,
  channel = 1,
  assay_name = NULL,
  tau1 = 0.75,
  tau2 = 2,
  basis = c("canonical", "derivative"),
  drift = 6
)
```

## Arguments

- x:

  A PhysioExperiment object with EDA data.

- events:

  A `data.frame` with numeric `onset` (seconds) and a `condition` column
  grouping events.

- conditions:

  Optional character vector selecting/ordering the conditions (default:
  the unique conditions in `events`).

- channel:

  EDA channel index or label to model (default: 1).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- tau1, tau2:

  SCR rise/decay time constants in seconds.

- basis:

  `"canonical"` (SCRF only) or `"derivative"` (SCRF plus its temporal
  derivative per condition).

- drift:

  Number of cosine drift terms (default: 6).

## Value

A list of class `"eda_glm"` with:

- coefficients:

  a `data.frame` with `condition`, `beta`, `se`, `t`, `p` (canonical
  amplitude per condition; `NA` for a condition dropped due to
  collinearity).

- derivative:

  (basis = "derivative") the derivative-regressor betas.

- fitted:

  fitted signal; `residuals`; `sigma`; `df`; `r_squared`; `rank`;
  `n_regressors`; `collinear`.

## References

Bach, D.R., & Friston, K.J. (2013). "Model-based analysis of skin
conductance responses." *Psychophysiology*, 50(1), 15-22.
[doi:10.1111/j.1469-8986.2012.01483.x](https://doi.org/10.1111/j.1469-8986.2012.01483.x)

## See also

[`edaDCM`](https://x-biosignal.github.io/PhysioEDA/reference/edaDCM.md)
for trial-wise latent input estimation,
[`scrf`](https://x-biosignal.github.io/PhysioEDA/reference/scrf.md) for
the response kernel.

## Examples

``` r
ev <- data.frame(onset = c(10, 30, 50, 70), condition = c("A", "A", "B", "B"))
pe <- edaSimulate(n_time = 1000, sr = 10, scr_count = 4, seed = 1)
fit <- edaGLM(pe, ev)
fit$coefficients
#>   condition        beta         se         t            p
#> 1         A  0.04408197 0.02909868  1.514913 1.301132e-01
#> 2         B -0.13266133 0.02937474 -4.516170 7.053846e-06
```
