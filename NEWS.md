# PhysioEDA 0.2.3

- `edaDecompose()` no longer errors when the input `PhysioExperiment` carries
  column names / `colData` (i.e. any object from a real reader such as `readWFDB`
  or `readEDF`): the new tonic/phasic assays are now stored with
  `assays(x, withDimnames = FALSE) <- ...`. Previously it stopped with
  "please use 'assay(x, withDimnames=FALSE) <-' ...". No change to the computed
  components.

# PhysioEDA 0.2.2

- `edaErscr()` no longer errors when a detected SCR peak falls at the very start of the analysis window: the backward onset search is now guarded (an unguarded `seq(peak - 1, start, by = -1)` had `from < to` and failed). The event-related vignette now runs clean. No change to results for peaks away from the window edge.

# PhysioEDA 0.2.1

- `edaDecompose()` now records a W3C-PROV provenance activity (shared `.recordProv` helper), so the decomposition step is visible to the reproducibility substrate run-tracing. No change to the tonic/phasic output.

# PhysioEDA 0.2.0

Initial release of PhysioEDA as a standalone package in the x-biosignal
ecosystem, split out from the PhysioExperiment monorepo. Provides a complete
electrodermal activity (EDA / skin conductance) analysis pipeline for
`PhysioExperiment` objects, from preprocessing through tonic/phasic
decomposition, SCR detection, feature extraction, and visualization.

## New Features

- Preprocessing:
  - `edaFilter()` applies FFT-based lowpass, highpass, and bandpass filtering
    with a smooth frequency-domain transition.
  - `edaDownsample()` decimates the signal after an anti-aliasing lowpass
    stage, preserving events and channel metadata.

- Artifact handling and quality:
  - `edaArtifact()` detects artifacts via threshold, gradient, and flatline
    methods, with optional correction by linear interpolation or NA
    replacement and a per-channel, per-method artifact summary.
  - `edaQuality()` reports per-channel quality metrics (statistics, flatline
    and artifact percentages, Donoho-estimator SNR) and an overall
    quality score and label.

- Tonic/phasic decomposition:
  - `edaDecompose()` separates skin conductance into tonic (SCL) and phasic
    (SCR) components via four methods: FFT `highpass`, sliding `median`,
    Continuous Decomposition Analysis (`cda`, Benedek & Kaernbach 2010), and
    convex optimization (`cvxeda`, Greco et al. 2016). CDA and cvxEDA also
    return the sudomotor nerve driver signal.

- SCR detection and features:
  - `edaPeaks()` detects SCR peaks by gradient zero-crossing or amplitude
    threshold, reporting onset, peak, amplitude, rise time, and 50% recovery
    time per response.
  - `edaFeatures()` summarizes per-channel SCR and SCL features including SCR
    count and rate, mean amplitude, tonic level, phasic AUC, and
    non-specific SCR count and frequency.
  - `edaErscr()` performs event-related SCR analysis, extracting
    amplitude, latency, rise time, and recovery time time-locked to
    experimental events.

- Transformations:
  - `edaTransform()` and `edaUntransform()` apply and invert `log`, `sqrt`,
    `zscore`, and `range` transformations, storing parameters in metadata for
    exact reversal.

- Simulation and data helpers:
  - `edaSimulate()` generates synthetic EDA with known tonic and phasic
    components using Bateman-shaped SCRs.
  - `make_eda()` and `make_eda_with_scr()` provide ready-made simulated
    objects (the latter with stimulus events) for testing and examples.

- Visualization:
  - `plotEda()`, `plotDecompose()`, and `plotPeaks()` render base-R plots of
    EDA time series, tonic/phasic/driver decomposition, and annotated SCR
    onsets and peaks.

## Bug Fixes

- Corrected the non-specific SCR (NS-SCR) frequency computation in
  `edaFeatures()` to count spontaneous SCRs exceeding the standard amplitude
  threshold and report them per minute.
