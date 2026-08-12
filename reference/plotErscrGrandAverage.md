# Grand-average event-related skin conductance response (ER-SCR)

Computes and plots the grand-average phasic skin-conductance waveform
time-locked to stimulus events, separately per event type (condition),
with a confidence-interval band. Each event's epoch over `window` is
baseline corrected against the pre-event `baseline` interval, epochs are
averaged within a condition, and the band is the standard error of the
mean scaled to `ci_level` (Boucsein 2012; Bach event-related SCR
averaging). The epoching reuses the event onsets from
[`PhysioCore::getEvents()`](https://rdrr.io/pkg/PhysioCore/man/getEvents.html)
on the same sample grid as
[`edaErscr()`](https://x-biosignal.github.io/PhysioEDA/reference/edaErscr.md).

## Usage

``` r
plotErscrGrandAverage(
  x,
  event_types = NULL,
  window = c(-1, 6),
  baseline = c(-1, 0),
  channel = 1,
  assay_name = NULL,
  ci = c("sem", "none"),
  ci_level = 0.95,
  colors = NULL,
  plot = TRUE,
  ...
)
```

## Arguments

- x:

  A PhysioExperiment object with a skin conductance assay and events.

- event_types:

  Character vector of event types to average and overlay; `NULL`
  (default) uses every event type present.

- window:

  Numeric `c(pre, post)` epoch limits in seconds relative to the event
  onset (default `c(-1, 6)`).

- baseline:

  Numeric `c(start, end)` pre-event baseline interval in seconds
  subtracted from each epoch, or `NULL` for no baseline correction
  (default `c(-1, 0)`).

- channel:

  Channel index to average (default 1).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- ci:

  Confidence band: `"sem"` (mean +/- z \* SEM, default) or `"none"`.

- ci_level:

  Confidence level for the SEM band (default 0.95).

- colors:

  Optional per-condition colours (recycled).

- plot:

  Draw the plot (default `TRUE`); set `FALSE` to only compute.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, an `erscr_grand_average` object with `grand_average` (a tidy
data frame of `type`, `time`, `mean`, `sd`, `sem`, `ci_lower`,
`ci_upper`, `n_trials`), the per-condition baseline-corrected `epochs`
matrices, and the `time` axis.

## References

Boucsein, W. (2012). Electrodermal Activity. 2nd ed. Springer.

## See also

[`edaErscr()`](https://x-biosignal.github.io/PhysioEDA/reference/edaErscr.md)

## Examples

``` r
x <- make_eda_with_scr()
ga <- plotErscrGrandAverage(x, plot = FALSE)
head(ga$grand_average)
#>                type time          mean          sd         sem     ci_lower
#> stimulus.1 stimulus -1.0  0.0276643418 0.035345176 0.017672588 -0.006973294
#> stimulus.2 stimulus -0.9  0.0188938019 0.033667749 0.016833875 -0.014099986
#> stimulus.3 stimulus -0.8  0.0138053383 0.023610669 0.011805335 -0.009332692
#> stimulus.4 stimulus -0.7  0.0078882259 0.013536197 0.006768099 -0.005377004
#> stimulus.5 stimulus -0.6 -0.0034103585 0.009070294 0.004535147 -0.012299083
#> stimulus.6 stimulus -0.5  0.0002261251 0.007947753 0.003973877 -0.007562530
#>               ci_upper n_trials
#> stimulus.1 0.062301978        4
#> stimulus.2 0.051887590        4
#> stimulus.3 0.036943369        4
#> stimulus.4 0.021153456        4
#> stimulus.5 0.005478366        4
#> stimulus.6 0.008014780        4
```
