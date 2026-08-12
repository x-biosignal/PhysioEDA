# Assess EDA Signal Quality

Computes per-channel signal quality metrics for electrodermal activity
data, including basic statistics, flatline detection, artifact
estimation, and an overall quality score.

## Usage

``` r
edaQuality(x, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- assay_name:

  Name of the input assay. If NULL, uses `defaultAssay(x)`.

## Value

A data.frame with one row per channel and the following columns:

- channel:

  Channel label

- mean_sc:

  Mean skin conductance

- sd_sc:

  Standard deviation of skin conductance

- min_sc:

  Minimum skin conductance

- max_sc:

  Maximum skin conductance

- pct_negative:

  Percentage of samples \<= 0

- pct_flatline:

  Percentage of signal in flatline segments (runs of consecutive
  near-zero differences longer than 1 second)

- pct_artifact:

  Percentage of gradient-based artifact samples

- snr_db:

  Estimated signal-to-noise ratio in dB (capped at 60)

- quality_score:

  Overall quality score from 0 to 100

- quality_label:

  "good" (\>= 70), "acceptable" (\>= 40), or "poor"

## References

Kleckner, I.R., et al. (2018). "Simple, transparent, and flexible
automated quality assessment procedures for ambulatory electrodermal
activity data." *IEEE Transactions on Biomedical Engineering*, 65(7),
1460-1467.
[doi:10.1109/TBME.2017.2758643](https://doi.org/10.1109/TBME.2017.2758643)

## See also

[`edaArtifact`](https://x-biosignal.github.io/PhysioEDA/reference/edaArtifact.md)
for artifact detection and correction,
[`edaFilter`](https://x-biosignal.github.io/PhysioEDA/reference/edaFilter.md)
for frequency-domain filtering,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition
