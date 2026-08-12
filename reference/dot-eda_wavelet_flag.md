# Wavelet-rule Artifact Flag per Epoch

Flags an epoch as an artifact when its fine-scale wavelet energy
(`wav_max_d1`, the peak level-1 detail coefficient) is a robust outlier:
greater than `median + k * MAD` across epochs. Motion and contact
artifacts inject high-frequency energy that genuine (slow) skin
conductance does not.

## Usage

``` r
.eda_wavelet_flag(features, k = 5)
```

## Arguments

- features:

  Epoch features from
  [`.eda_epoch_features`](https://x-biosignal.github.io/PhysioEDA/reference/dot-eda_epoch_features.md).

- k:

  Robust-outlier multiplier (default: 5).

## Value

Logical vector, one per epoch (`TRUE` = artifact).
