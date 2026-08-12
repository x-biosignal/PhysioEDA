# CDA/DDA Phasic-Driver Metrics

Computes the Ledalab-style phasic-driver metrics from a (non-negative)
sudomotor driver signal: the number of significant SCRs, the sum of
their amplitudes, and the integrated SCR (the time-integral of the
phasic driver). Significant SCRs are contiguous supra-threshold driver
excursions; excursions separated by less than `min_interval_sec` are
merged into one SCR.

## Usage

``` r
.cda_metrics(
  driver,
  sr,
  amp_threshold = 0.01,
  min_interval_sec = 1,
  baseline_sec = 4,
  edge_sec = 1
)
```

## Arguments

- driver:

  Numeric vector of the (non-negative) sudomotor driver. The slow
  tonic-driver baseline is removed (a running minimum over
  `baseline_sec`) to obtain the phasic driver before detecting SCRs.

- sr:

  Sampling rate in Hz.

- amp_threshold:

  Minimum phasic-driver amplitude for a significant SCR (default: 0.01).

- min_interval_sec:

  Minimum separation in seconds between distinct SCRs; closer excursions
  are merged (default: 1).

- baseline_sec:

  Window in seconds for the running-minimum tonic-driver baseline that
  is subtracted to form the phasic driver (default: 4).

- edge_sec:

  Seconds of phasic driver at each end to zero out before detection,
  suppressing deconvolution boundary artifacts (default: 1).

## Value

A list with:

- nSCR:

  Number of significant SCRs.

- AmpSum:

  Sum of the SCR peak amplitudes (uS).

- ISCR:

  Integrated SCR: time-integral of the phasic driver (`uS*s`).

- onset_sec:

  Onset times of the detected SCRs (seconds).

- peak_sec:

  Peak times of the detected SCRs (seconds).

- amplitudes:

  Peak amplitudes of the detected SCRs (uS).

- scr_area:

  Time-integral of the driver within each SCR (`uS*s`).

## References

Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
electrodermal activity." *Journal of Neuroscience Methods*, 190(1),
80-91.
[doi:10.1016/j.jneumeth.2010.04.028](https://doi.org/10.1016/j.jneumeth.2010.04.028)

## See also

[`.cda_decompose`](https://x-biosignal.github.io/PhysioEDA/reference/dot-cda_decompose.md),
[`edaFeatures`](https://x-biosignal.github.io/PhysioEDA/reference/edaFeatures.md)
