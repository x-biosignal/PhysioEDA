# Load the Shipped EDA Artifact Logistic Model

Reads the lightweight ridge-logistic artifact model shipped in
`inst/extdata/eda_artifact_model.csv` (columns `feature`, `mean`, `sd`,
`coef`); the result is cached.

## Usage

``` r
.load_eda_artifact_model()
```

## Value

A `data.frame` with the model coefficients and standardisation.
