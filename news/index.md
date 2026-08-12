# Changelog

## PhysioEDA 0.2.2

- [`edaErscr()`](https://x-biosignal.github.io/PhysioEDA/reference/edaErscr.md)
  no longer errors when a detected SCR peak falls at the very start of
  the analysis window: the backward onset search is now guarded (an
  unguarded `seq(peak - 1, start, by = -1)` had `from < to` and failed).
  The event-related vignette now runs clean. No change to results for
  peaks away from the window edge.

## PhysioEDA 0.2.1

- [`edaDecompose()`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
  now records a W3C-PROV provenance activity (shared `.recordProv`
  helper), so the decomposition step is visible to the reproducibility
  substrate run-tracing. No change to the tonic/phasic output.

## PhysioEDA 0.2.0

Initial release of PhysioEDA as a standalone package in the x-biosignal
ecosystem, split out from the PhysioExperiment monorepo. Provides a
complete electrodermal activity (EDA / skin conductance) analysis
pipeline for `PhysioExperiment` objects, from preprocessing through
tonic/phasic decomposition, SCR detection, feature extraction, and
visualization.

### New Features

- Preprocessing:
  - [`edaFilter()`](https://x-biosignal.github.io/PhysioEDA/reference/edaFilter.md)
    applies FFT-based lowpass, highpass, and bandpass filtering with a
    smooth frequency-domain transition.
  - [`edaDownsample()`](https://x-biosignal.github.io/PhysioEDA/reference/edaDownsample.md)
    decimates the signal after an anti-aliasing lowpass stage,
    preserving events and channel metadata.
- Artifact handling and quality:
  - [`edaArtifact()`](https://x-biosignal.github.io/PhysioEDA/reference/edaArtifact.md)
    detects artifacts via threshold, gradient, and flatline methods,
    with optional correction by linear interpolation or NA replacement
    and a per-channel, per-method artifact summary.
  - [`edaQuality()`](https://x-biosignal.github.io/PhysioEDA/reference/edaQuality.md)
    reports per-channel quality metrics (statistics, flatline and
    artifact percentages, Donoho-estimator SNR) and an overall quality
    score and label.
- Tonic/phasic decomposition:
  - [`edaDecompose()`](https://x-biosignal.github.io/PhysioEDA/reference/edaDecompose.md)
    separates skin conductance into tonic (SCL) and phasic (SCR)
    components via four methods: FFT `highpass`, sliding `median`,
    Continuous Decomposition Analysis (`cda`, Benedek & Kaernbach 2010),
    and convex optimization (`cvxeda`, Greco et al. 2016). CDA and
    cvxEDA also return the sudomotor nerve driver signal.
- SCR detection and features:
  - [`edaPeaks()`](https://x-biosignal.github.io/PhysioEDA/reference/edaPeaks.md)
    detects SCR peaks by gradient zero-crossing or amplitude threshold,
    reporting onset, peak, amplitude, rise time, and 50% recovery time
    per response.
  - [`edaFeatures()`](https://x-biosignal.github.io/PhysioEDA/reference/edaFeatures.md)
    summarizes per-channel SCR and SCL features including SCR count and
    rate, mean amplitude, tonic level, phasic AUC, and non-specific SCR
    count and frequency.
  - [`edaErscr()`](https://x-biosignal.github.io/PhysioEDA/reference/edaErscr.md)
    performs event-related SCR analysis, extracting amplitude, latency,
    rise time, and recovery time time-locked to experimental events.
- Transformations:
  - [`edaTransform()`](https://x-biosignal.github.io/PhysioEDA/reference/edaTransform.md)
    and
    [`edaUntransform()`](https://x-biosignal.github.io/PhysioEDA/reference/edaUntransform.md)
    apply and invert `log`, `sqrt`, `zscore`, and `range`
    transformations, storing parameters in metadata for exact reversal.
- Simulation and data helpers:
  - [`edaSimulate()`](https://x-biosignal.github.io/PhysioEDA/reference/edaSimulate.md)
    generates synthetic EDA with known tonic and phasic components using
    Bateman-shaped SCRs.
  - [`make_eda()`](https://x-biosignal.github.io/PhysioEDA/reference/make_eda.md)
    and
    [`make_eda_with_scr()`](https://x-biosignal.github.io/PhysioEDA/reference/make_eda_with_scr.md)
    provide ready-made simulated objects (the latter with stimulus
    events) for testing and examples.
- Visualization:
  - [`plotEda()`](https://x-biosignal.github.io/PhysioEDA/reference/plotEda.md),
    [`plotDecompose()`](https://x-biosignal.github.io/PhysioEDA/reference/plotDecompose.md),
    and
    [`plotPeaks()`](https://x-biosignal.github.io/PhysioEDA/reference/plotPeaks.md)
    render base-R plots of EDA time series, tonic/phasic/driver
    decomposition, and annotated SCR onsets and peaks.

### Bug Fixes

- Corrected the non-specific SCR (NS-SCR) frequency computation in
  [`edaFeatures()`](https://x-biosignal.github.io/PhysioEDA/reference/edaFeatures.md)
  to count spontaneous SCRs exceeding the standard amplitude threshold
  and report them per minute.
