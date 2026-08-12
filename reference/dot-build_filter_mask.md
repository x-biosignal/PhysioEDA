# Build FFT filter mask

Constructs a smooth frequency-domain mask for lowpass, highpass, or
bandpass filtering using a power-law (raised cosine-like) transition.

## Usage

``` r
.build_filter_mask(freq_bins, type, cutoff, order)
```

## Arguments

- freq_bins:

  Numeric vector of frequencies for each FFT bin.

- type:

  Filter type.

- cutoff:

  Cutoff frequency or frequencies.

- order:

  Steepness parameter.

## Value

Numeric vector of mask values (0 to 1) for each FFT bin.
