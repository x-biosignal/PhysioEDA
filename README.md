# PhysioEDA <img src="man/figures/logo.png" align="right" height="139" alt="PhysioEDA logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/x-biosignal/PhysioEDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/x-biosignal/PhysioEDA/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/PhysioEDA)](https://CRAN.R-project.org/package=PhysioEDA)
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioEDA)](https://x-biosignal.r-universe.dev/PhysioEDA)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Electrodermal Activity Analysis for PhysioExperiment Objects**

PhysioEDA provides electrodermal activity (EDA / skin conductance) analysis functions for the PhysioExperiment ecosystem. With 19 exported functions, it delivers a complete EDA processing pipeline -- from preprocessing (filtering, downsampling, artifact handling) through tonic/phasic decomposition (highpass, median, CDA, cvxEDA), SCR peak detection, feature extraction, event-related SCR analysis, and visualization. All functions operate on PhysioExperiment objects and store results as new assays for seamless integration with downstream analysis.

## Installation

You can install PhysioEDA from [r-universe](https://x-biosignal.r-universe.dev):

```r
install.packages("PhysioEDA",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioEDA")
```

## Quick Start

```r
library(PhysioEDA)

# Generate a simulated EDA signal (10 minutes at 10 Hz, 5 embedded SCRs)
pe <- edaSimulate(n_time = 6000, sr = 10, scr_count = 5, seed = 42)

# Lowpass filter to remove high-frequency noise
pe <- edaFilter(pe, type = "lowpass", cutoff = 1.0)

# Decompose into tonic (SCL) and phasic (SCR) components
pe <- edaDecompose(pe, method = "highpass", cutoff = 0.05)

# Detect SCR peaks in the phasic signal
peaks <- edaPeaks(pe, amplitude_min = 0.01)
head(peaks)
#>   channel onset_sample onset_sec peak_sample peak_sec amplitude rise_time recovery_time
#> 1    EDA1          587      58.6         602     60.1     0.482      1.50          3.20
#> ...

# Extract comprehensive features per channel
features <- edaFeatures(pe, peaks)
features
#>   channel scr_count scr_rate_per_min mean_amplitude mean_scl scl_sd auc_phasic ns_scr_freq
#> 1    EDA1         5             0.50          0.453     5.02   0.21       4.32        0.50

# Visualize decomposition
plotDecompose(pe)
```

## Features

### Preprocessing

Frequency-domain filtering and downsampling with anti-aliasing:

- `edaFilter()` -- FFT-based lowpass, highpass, or bandpass filter with smooth frequency-domain transition
- `edaDownsample()` -- reduce sampling rate with automatic anti-aliasing lowpass filter and event preservation

### Artifact Handling

Detection and correction of motion artifacts and signal dropouts:

- `edaArtifact()` -- detect artifacts (amplitude spikes, rapid transients) and correct them via interpolation or replacement

### Tonic/Phasic Decomposition

Separate the slow-varying skin conductance level (SCL) from fast-varying skin conductance responses (SCR) using four methods:

- `edaDecompose()` with `method = "highpass"` -- FFT-based highpass/lowpass separation at a configurable cutoff frequency
- `edaDecompose()` with `method = "median"` -- sliding median filter for tonic extraction
- `edaDecompose()` with `method = "cda"` -- Continuous Decomposition Analysis (Benedek & Kaernbach, 2010) using Bateman impulse response deconvolution with non-negativity constraint
- `edaDecompose()` with `method = "cvxeda"` -- convex optimization approach (Greco et al., 2016) using iterative ADMM with Wiener deconvolution and L1 sparsity

### SCR Peak Detection

Detect individual skin conductance responses with configurable criteria:

- `edaPeaks()` -- detect SCR peaks using gradient-based (first-derivative zero-crossing) or amplitude threshold methods
- Returns onset, peak, amplitude, rise time, and 50% recovery time for each SCR
- Configurable minimum amplitude, rise time bounds, and automatic phasic assay selection

### Feature Extraction

Comprehensive per-channel summary metrics:

- `edaFeatures()` -- SCR count, SCR rate per minute, mean SCR amplitude, mean SCL, SCL standard deviation, phasic area under curve (AUC), non-specific SCR frequency

### Signal Quality Assessment

- `edaQuality()` -- assess EDA signal quality per channel

### Data Transformations

Normalize and transform EDA signals for statistical analysis:

- `edaTransform()` -- apply log, sqrt, range normalization, or z-score transformations
- `edaUntransform()` -- reverse transformations to recover original scale

### Event-Related Analysis

Time-locked SCR analysis for experimental paradigms:

- `edaErscr()` -- event-related SCR analysis with configurable onset and peak windows, returning per-event amplitude, latency, rise time, and recovery time

### Visualization

Publication-ready plotting functions using base R graphics:

- `plotEda()` -- plot raw or processed EDA signal
- `plotPeaks()` -- overlay detected SCR peaks on the signal
- `plotDecompose()` -- multi-panel plot showing original signal, tonic, and phasic components

### Simulated EDA Data

Synthetic EDA generators with known ground truth for testing and validation:

- `edaSimulate()` -- generate EDA with configurable SCR count, amplitude, SCL level, and noise using biexponential (Bateman function) impulse responses
- `make_eda()` -- create a basic simulated EDA PhysioExperiment
- `make_eda_with_scr()` -- create EDA with embedded SCR events at specified times

## Dependencies

- **R** (>= 4.2)
- **[PhysioCore](https://github.com/x-biosignal/PhysioCore)**
- **SummarizedExperiment**
- **S4Vectors**
- **stats**, **graphics**, **grDevices**

## PhysioExperiment Ecosystem

PhysioEDA is the electrodermal activity analysis module of the PhysioExperiment ecosystem, a suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|---------|-------------|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| [PhysioECG](https://github.com/x-biosignal/PhysioECG) | ECG analysis and HRV |
| **PhysioEDA** | Electrodermal activity analysis |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to browse all available packages.

## References

- Boucsein, W. (2012). *Electrodermal Activity*. 2nd ed. Springer.
- Benedek, M. & Kaernbach, C. (2010). "A continuous measure of phasic electrodermal activity." *Journal of Neuroscience Methods*, 190(1), 80--91.
- Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to electrodermal activity processing." *IEEE Transactions on Biomedical Engineering*, 63(4), 797--804.
- Bach, D.R., et al. (2010). "Modelling event-related skin conductance responses." *International Journal of Psychophysiology*, 75(3), 349--356.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

Yusuke Matsui
