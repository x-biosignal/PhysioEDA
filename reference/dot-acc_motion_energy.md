# Per-epoch Accelerometer Motion Energy

Motion energy per epoch: the variance of the accelerometer magnitude
within the epoch. A still sensor sees near-constant magnitude (low
variance); motion raises it.

## Usage

``` r
.acc_motion_energy(acc, sr, epoch_sec = 5)
```

## Arguments

- acc:

  Numeric matrix (time x axes) of accelerometer channels.

- sr:

  Sampling rate in Hz.

- epoch_sec:

  Epoch length in seconds (default: 5).

## Value

A `data.frame` of epoch bounds plus `motion_energy`.
