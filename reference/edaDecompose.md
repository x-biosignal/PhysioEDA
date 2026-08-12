# Decompose EDA into Tonic and Phasic Components

Separates an electrodermal activity (EDA) signal into its slow-varying
tonic component (skin conductance level, SCL) and fast-varying phasic
component (skin conductance responses, SCR).

## Usage

``` r
edaDecompose(
  x,
  method = c("highpass", "median", "cda", "dda", "cvxeda", "cvxeda_fast"),
  cutoff = 0.05,
  window_sec = 4,
  tau1 = 0.75,
  tau2 = 2,
  alpha = 0.01,
  gamma = 0.1,
  delta_knot = 10,
  qp_sr = 4,
  optimize_tau = FALSE,
  assay_name = NULL,
  output_tonic = "tonic",
  output_phasic = "phasic"
)
```

## Arguments

- x:

  A PhysioExperiment object containing EDA data.

- method:

  Decomposition method: `"highpass"`, `"median"`, `"cda"`, `"dda"`,
  `"cvxeda"`, or `"cvxeda_fast"`. Default is `"highpass"`.

- cutoff:

  Cutoff frequency in Hz for the highpass method (default: 0.05).

- window_sec:

  Window length in seconds for the median method (default: 4).

- tau1:

  SCR rise time constant in seconds for CDA/cvxEDA (default: 0.75).

- tau2:

  SCR decay time constant in seconds for CDA/cvxEDA (default: 2.0).

- alpha:

  L1 sparsity penalty for cvxEDA (default: 0.01). Ignored by other
  methods.

- gamma:

  Smoothness weight for cvxEDA tonic component (default: 0.1). Ignored
  by other methods.

- delta_knot:

  Tonic spline knot spacing in seconds for the `"cvxeda"` QP (default:
  10).

- qp_sr:

  Target rate (Hz) at which the `"cvxeda"` QP is solved; the signal is
  decimated towards this rate for tractability (default: 4).

- optimize_tau:

  Logical; for `"cda"` / `"dda"`, optimise the decay constant `tau2` to
  minimise the driver's negative energy (default: `FALSE`).

- assay_name:

  Name of the input assay. If NULL, uses `defaultAssay(x)`.

- output_tonic:

  Name for the tonic output assay (default: "tonic").

- output_phasic:

  Name for the phasic output assay (default: "phasic").

## Value

A modified
[`PhysioExperiment`](https://rdrr.io/pkg/PhysioCore/man/PhysioExperiment.html)
with new assays:

- tonic:

  The slow-varying skin conductance level (SCL) component.

- phasic:

  The fast-varying skin conductance response (SCR) component.

- driver:

  (CDA and cvxEDA variants only) The sudomotor nerve activity driver
  signal.

Decomposition parameters are stored in `metadata(x)$eda_decompose`; for
`method = "cvxeda"` this also records the per-channel solver `engine`
("qp" or "wiener_fallback") and final QP `objective`.

## Details

Five methods are available:

- highpass:

  FFT-based highpass/lowpass separation at a cutoff frequency.

- median:

  Sliding median filter for tonic extraction.

- cda:

  Continuous Decomposition Analysis (Benedek & Kaernbach, 2010).
  Deconvolution with a Bateman impulse response, Gaussian smoothing, and
  a non-negativity constraint on the driver signal; optional `tau2`
  optimisation.

- dda:

  Discrete Decomposition Analysis (Benedek & Kaernbach, 2010): the
  continuous CDA driver is reduced to a sparse impulse train (one
  impulse per significant SCR) and reconvolved.

- cvxeda:

  Exact convex optimization decomposition (Greco et al., 2016), solved
  as a quadratic program: a cubic-knot spline tonic basis and a sparse,
  non-negative sudomotor driver through a biexponential system, with L1
  sparsity on the driver and L2 smoothness on the tonic. Requires the
  quadprog solver; degrades gracefully to the Wiener approximation if it
  is unavailable.

- cvxeda_fast:

  Fast Wiener/ADMM approximation of cvxEDA (no solver dependency); the
  original iterative approximation.

## References

Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
electrodermal activity." *Journal of Neuroscience Methods*, 190(1),
80-91.
[doi:10.1016/j.jneumeth.2010.04.028](https://doi.org/10.1016/j.jneumeth.2010.04.028)

Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
electrodermal activity processing." *IEEE Transactions on Biomedical
Engineering*, 63(4), 797-804.
[doi:10.1109/TBME.2015.2474131](https://doi.org/10.1109/TBME.2015.2474131)

## See also

[`edaPeaks`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md)
for SCR peak detection on the phasic signal,
[`edaFeatures`](https://x-biosignal.github.io/PhysioEDA/reference/edaFeatures.md)
for feature extraction,
[`plotDecompose`](https://x-biosignal.github.io/PhysioEDA/reference/plotDecompose.md)
for visualizing decomposition results
