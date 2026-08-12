# Per-epoch Artifact Probability from the Logistic Model

Applies the shipped standardised ridge-logistic model to the epoch
features.

## Usage

``` r
.eda_ml_predict(features, model)
```

## Arguments

- features:

  Epoch features from
  [`.eda_epoch_features`](https://x-biosignal.github.io/PhysioEDA/reference/dot-eda_epoch_features.md).

- model:

  The model `data.frame` from
  [`.load_eda_artifact_model`](https://x-biosignal.github.io/PhysioEDA/reference/dot-load_eda_artifact_model.md).

## Value

Numeric vector of artifact probabilities, one per epoch.
