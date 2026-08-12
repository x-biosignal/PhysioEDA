# Reverse EDA Transformation

Reverses a transformation previously applied by
[`edaTransform`](https://x-biosignal.github.io/PhysioEDA/reference/edaTransform.md),
using parameters stored in `metadata(x)$eda_transform`.

## Usage

``` r
edaUntransform(x, assay_name = NULL, output_assay = NULL)
```

## Arguments

- x:

  A PhysioExperiment object that has been transformed via
  `edaTransform`.

- assay_name:

  Name of the transformed assay to reverse. If NULL, uses the
  output_assay recorded in the transform metadata.

- output_assay:

  Name for the untransformed output assay. If NULL, defaults to
  `"untransformed"`.

## Value

A modified PhysioExperiment with the untransformed assay added.

## References

Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
[doi:10.1007/978-1-4614-1126-0](https://doi.org/10.1007/978-1-4614-1126-0)

## See also

[`edaTransform`](https://x-biosignal.github.io/PhysioEDA/reference/edaTransform.md)
for applying the initial transformation,
[`edaDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
for tonic/phasic decomposition
