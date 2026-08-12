# Extract EDA Features per Channel

Computes summary features from electrodermal activity data, including
skin conductance response (SCR) statistics and skin conductance level
(SCL) measures. Requires tonic/phasic decomposition to have been
performed (via
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)).

## Usage

``` r
edaFeatures(
  x,
  peaks = NULL,
  window = NULL,
  assay_name = NULL,
  ns_scr_amplitude_min = 0.05
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- peaks:

  Optional pre-computed peaks data.frame from
  [`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md).
  If NULL, peaks are computed automatically using the gradient method.

- window:

  Optional numeric vector of length 2 giving the time window in seconds
  as `c(start_sec, end_sec)` to restrict analysis.

- assay_name:

  Name of the input assay for peak detection. If NULL, uses
  `defaultAssay(x)`.

- ns_scr_amplitude_min:

  Minimum SCR amplitude (in uS) to count a peak as a non-specific SCR
  (default 0.05, the standard threshold).

## Value

A `data.frame` with one row per channel and the following columns:

- channel:

  Character channel label.

- scr_count:

  Integer number of detected SCR peaks.

- scr_rate_per_min:

  Numeric SCR count per minute.

- mean_amplitude:

  Numeric mean SCR peak amplitude in microsiemens (0 if no peaks
  detected).

- mean_scl:

  Numeric mean tonic skin conductance level.

- scl_sd:

  Numeric standard deviation of the tonic signal.

- auc_phasic:

  Numeric area under the phasic curve (positive values only, in `uS*s`),
  or `NA` if no phasic assay exists.

- ns_scr_count:

  Integer count of non-specific SCRs (amplitude \>=
  `ns_scr_amplitude_min`).

- ns_scr_freq:

  Numeric non-specific SCR frequency (per minute):
  `ns_scr_count / duration_min`.

- cda_nscr:

  Integer number of significant SCRs in the CDA/DDA driver, or `NA` if
  no driver assay exists.

- cda_ampsum:

  Numeric sum of the CDA/DDA driver SCR amplitudes (uS), or `NA`.

- cda_iscr:

  Numeric integrated SCR: time-integral of the phasic driver (`uS*s`),
  or `NA`.

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
electrodermal activity." *Journal of Neuroscience Methods*, 190(1),
80-91.
[doi:10.1016/j.jneumeth.2010.04.028](https://doi.org/10.1016/j.jneumeth.2010.04.028)

## See also

[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition (run first),
[`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md)
for SCR peak detection,
[`edaErscr`](https://x-biosignal.github.io/PhysioEDA/reference/edaErscr.md)
for event-related SCR analysis
