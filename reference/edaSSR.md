# Sympathetic skin response (SSR) analysis

Extracts the sympathetic skin response evoked by each stimulus from a
skin conductance recording: the response onset latency, the
baseline-to-peak amplitude, the waveform morphology (monophasic /
biphasic) and a present/absent classification. A per-channel summary
reports the response rate, mean latency and amplitude, an
amplitude-habituation slope across the repeated stimuli, and the
clinical absent-SSR flag.

## Usage

``` r
edaSSR(
  x,
  stimulus_events,
  assay_name = NULL,
  response_window = c(0.5, 8),
  baseline_sec = 2,
  amplitude_min = 0.05,
  recovery_frac = 0.1,
  onset_frac = 0.05,
  morphology_frac = 0.3,
  events_unit = c("auto", "index", "seconds")
)
```

## Arguments

- x:

  A PhysioExperiment object with a skin conductance assay.

- stimulus_events:

  Stimulus times, as sample indices or times in seconds (see
  `events_unit`).

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- response_window:

  Numeric `c(start, end)` seconds after each stimulus in which the SSR
  is sought (default `c(0.5, 8)`).

- baseline_sec:

  Length in seconds of the pre-stimulus window used to estimate and
  remove the tonic trend (default 2).

- amplitude_min:

  Minimum detrended peak amplitude (uS) for a response to count as
  present (default 0.05).

- recovery_frac:

  Fraction of the peak amplitude the trace must fall back after the peak
  for the response to be accepted as phasic (rather than a monotone
  drift; default 0.1).

- onset_frac:

  Fraction of the SSR amplitude above baseline at which the response
  onset (and hence latency) is marked (default 0.05).

- morphology_frac:

  Fraction of the SSR amplitude an opposite-going deflection must reach
  for the response to be classed biphasic (default 0.3).

- events_unit:

  Unit of `stimulus_events`: `"auto"` (default), `"index"` or
  `"seconds"`.

## Value

An `eda_ssr` object with `responses` (one row per stimulus per channel:
`channel`, `stimulus`, `stim_time`, `present`, `latency`, `amplitude`,
`peak_time`, `morphology`) and `summary` (one row per channel:
`channel`, `n_stimuli`, `n_present`, `response_rate`, `mean_latency`,
`mean_amplitude`, `habituation_slope`, `ssr_absent`).

## Details

The tonic trend estimated from the pre-stimulus baseline is subtracted
from each response window, so slow skin-conductance drift is neither
counted as a response nor mistaken for a biphasic deflection. In the
detrended window the phasic peak is located; a response is accepted only
if it exceeds `amplitude_min` and recovers (rises then falls), and the
latency is measured from the stimulus to the onset of that phasic
deflection (Vetrugno et al. 2003).

## References

Vetrugno, R. et al. (2003). Sympathetic skin response: basic mechanisms
and clinical applications. Clin Auton Res 13(4):256-270.

## See also

[`edaSympatheticIndex()`](https://x-biosignal.github.io/PhysioEDA/reference/edaSympatheticIndex.md)

## Examples

``` r
# 3 stimuli at 10, 40 and 70 s, each followed by an SCR
pe <- edaSimulate(n_time = 900, sr = 10, scr_count = 0, seed = 1)
sig <- as.numeric(SummarizedExperiment::assay(pe))
sr <- 10
for (st in c(10, 40, 70)) {
  t <- (seq_along(sig) - 1) / sr - st
  sig <- sig + 0.5 * ifelse(t > 0, (exp(-t / 3) - exp(-t / 0.75)), 0)
}
SummarizedExperiment::assay(pe) <- matrix(sig, ncol = 1)
res <- edaSSR(pe, stimulus_events = c(10, 40, 70), events_unit = "seconds")
res$summary
#>   channel n_stimuli n_present response_rate mean_latency mean_amplitude
#> 1    EDA1         3         3             1          0.5       0.243441
#>   habituation_slope ssr_absent
#> 1      -0.001389532      FALSE
```
