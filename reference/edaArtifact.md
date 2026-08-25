# Detect and Correct EDA Artifacts

Identifies artifacts in electrodermal activity (EDA) signals using one
or more detection methods (threshold, gradient, flatline, and
epoch-based wavelet / machine-learning methods) and optionally corrects
them via interpolation or NA replacement. An accelerometer channel can
additionally gate out high-motion epochs (Taylor et al., 2015).

## Usage

``` r
edaArtifact(
  x,
  methods = c("threshold", "gradient", "flatline"),
  threshold_range = c(0.001, 60),
  gradient_max = NULL,
  flatline_sec = 5,
  epoch_sec = 5,
  wavelet_k = 5,
  ml_threshold = 0.5,
  acc_channel = NULL,
  acc_threshold = NULL,
  correct = "interpolate",
  assay_name = NULL,
  output_assay = "cleaned"
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- methods:

  Character vector of detection methods to apply. Any combination of
  `"threshold"`, `"gradient"`, `"flatline"`, `"wavelet"` (per-5 s-epoch
  stationary-wavelet fine-scale energy outliers), and `"ml"` (a shipped
  ridge-logistic classifier on wavelet/derivative epoch features).
  Default is the first three.

- threshold_range:

  Numeric vector of length 2 giving the acceptable range of EDA values
  in microsiemens (default: `c(0.001, 60)`). Values outside this range
  are flagged as artifacts.

- gradient_max:

  Maximum allowable absolute gradient in microsiemens per sample. If
  NULL (default), computed as `10 / samplingRate(x)` (equivalent to 10
  uS/sec).

- flatline_sec:

  Minimum duration in seconds of a constant-value segment to be flagged
  as a flatline artifact (default: 5).

- epoch_sec:

  Epoch length in seconds for the `"wavelet"` / `"ml"` methods and
  accelerometer gating (default: 5).

- wavelet_k:

  Robust-outlier multiplier for the `"wavelet"` method; an epoch is
  flagged when its peak level-1 detail coefficient exceeds
  `median + wavelet_k * MAD` across epochs (default: 5).

- ml_threshold:

  Probability threshold for the `"ml"` method's binary label (default:
  0.5).

- acc_channel:

  Optional channel index or label (or several) identifying accelerometer
  channel(s). When supplied, per-epoch motion energy (variance of the
  accelerometer magnitude) gates the EDA channels: epochs with energy
  above `acc_threshold` are flagged. The accelerometer channels are
  excluded from EDA artifact detection.

- acc_threshold:

  Motion-energy threshold for accelerometer gating. If `NULL` (default),
  computed robustly as `median + 3 * MAD` of the per-epoch motion
  energies.

- correct:

  Correction strategy: `"interpolate"` (linear interpolation across
  artifact regions), `"na"` (replace with NA), or `"none"` (detection
  only, no correction). Default is `"interpolate"`.

- assay_name:

  Name of the input assay. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name for the corrected output assay (default: `"cleaned"`). Only used
  when `correct != "none"`.

## Value

A modified
[`PhysioExperiment`](https://x-biosignal.github.io/PhysioCore//reference/PhysioExperiment.html)
with artifact information stored in `metadata(x)$eda_artifacts`, a list
containing:

- mask:

  Logical matrix (time x channels) where `TRUE` indicates an artifact
  sample.

- summary:

  A `data.frame` with columns `channel`, `method`, `n_artifacts`, and
  `pct` giving artifact counts per channel per detection method.

- epochs:

  (only when `"wavelet"`/`"ml"` or accelerometer gating is used) A
  `data.frame`, one row per EDA channel per epoch, with the epoch bounds
  and, as applicable, `wavelet_flag`, `ml_prob`, `ml_flag`,
  `motion_energy`, and `acc_flag`.

If `correct != "none"`, the corrected signal is stored in the
`output_assay`.

## References

Taylor, S., et al. (2015). "Automatic identification of artifacts in
electrodermal activity data." *IEEE EMBC*, 1934-1937.
[doi:10.1109/EMBC.2015.7318762](https://doi.org/10.1109/EMBC.2015.7318762)

Kleckner, I.R., et al. (2018). "Simple, transparent, and flexible
automated quality assessment procedures for ambulatory electrodermal
activity data." *IEEE Transactions on Biomedical Engineering*, 65(7),
1460-1467.
[doi:10.1109/TBME.2017.2758643](https://doi.org/10.1109/TBME.2017.2758643)

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

## See also

[`edaQuality`](https://x-biosignal.github.io/PhysioEDA/reference/edaQuality.md)
for signal quality assessment,
[`edaFilter`](https://x-biosignal.github.io/PhysioEDA/reference/edaFilter.md)
for frequency-domain filtering,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition
