# Transform EDA Data

Applies a mathematical transformation to EDA signal data. Common
transformations include log and square root (to reduce skewness of SCR
amplitudes), z-score standardization, and range normalization.
Transformation parameters are stored in metadata so that
[`edaUntransform`](https://x-biosignal.github.io/PhysioEDA/reference/edaUntransform.md)
can reverse the operation.

## Usage

``` r
edaTransform(
  x,
  method = c("log", "sqrt", "zscore", "range"),
  assay_name = NULL,
  output_assay = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- method:

  Transformation method: `"log"`, `"sqrt"`, `"zscore"`, or `"range"`.

- assay_name:

  Name of the input assay. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name for the output assay. If NULL, defaults to
  `"{method}_transformed"`.

## Value

A modified PhysioExperiment with the transformed assay added and
transformation parameters stored in `metadata(x)$eda_transform`.

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

## See also

[`edaUntransform`](https://x-biosignal.github.io/PhysioEDA/reference/edaUntransform.md)
for reversing the transformation,
[`edaFeatures`](https://x-biosignal.github.io/PhysioEDA/reference/edaFeatures.md)
for feature extraction,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition
