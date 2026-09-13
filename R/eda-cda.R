#' CDA Decomposition for a Single EDA Channel
#'
#' Performs Continuous Decomposition Analysis (CDA) on a single EDA channel
#' using the method described by Benedek & Kaernbach (2010). The signal is
#' deconvolved with a Bateman impulse response function to obtain a sudomotor
#' nerve activity (driver) signal, which is then reconvolved to separate tonic
#' and phasic components.
#'
#' @param signal Numeric vector of EDA values.
#' @param sr Sampling rate in Hz.
#' @param tau1 Rise time constant in seconds (default: 0.75).
#' @param tau2 Decay time constant in seconds (default: 2.0).
#' @param optimize_tau Logical; if \code{TRUE}, optimise the decay constant
#'   \code{tau2} to minimise the negative energy of the deconvolved driver
#'   (a Ledalab-style criterion; the driver should be sparse and non-negative).
#'   \code{tau1} is held fixed. Default \code{FALSE}.
#' @return A list with components:
#'   \describe{
#'     \item{tonic}{Numeric vector of the tonic (SCL) component.}
#'     \item{phasic}{Numeric vector of the phasic (SCR) component.}
#'     \item{driver}{Numeric vector of the sudomotor nerve activity driver signal.}
#'     \item{tau1}{Rise time constant used.}
#'     \item{tau2}{Decay time constant used (optimised when \code{optimize_tau}).}
#'   }
#' @references
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#' @seealso \code{\link{edaDecompose}} for the user-facing decomposition wrapper,
#'   \code{\link{.dda_decompose}} for the discrete (DDA) alternative,
#'   \code{\link{.cvxeda_decompose}} for the cvxEDA alternative
#' @keywords internal
.cda_decompose <- function(signal, sr, tau1 = 0.75, tau2 = 2.0,
                           optimize_tau = FALSE) {
  stopifnot(is.numeric(signal), length(signal) > 1)
  stopifnot(is.numeric(sr), sr > 0)
  stopifnot(is.numeric(tau1), tau1 > 0)
  stopifnot(is.numeric(tau2), tau2 > 0, tau2 > tau1)

  n <- length(signal)

  # Raw FFT deconvolution of the signal by a unit-peak Bateman kernel.
  S <- stats::fft(signal)
  deconv_raw <- function(t2) {
    hl <- min(n, as.integer(10 * t2 * sr))
    th <- seq(0, by = 1 / sr, length.out = hl)
    hk <- scrf(th, tau1 = tau1, tau2 = t2, form = "bateman", normalize = "none")
    hk <- hk / max(hk)
    Hk <- stats::fft(c(hk, rep(0, n - hl)))
    eps <- max(Mod(Hk)) * 1e-6
    list(driver = Re(stats::fft(S / (Hk + eps), inverse = TRUE)) / n,
         h = hk, h_len = hl, H = Hk)
  }

  # Optional tau2 optimisation: minimise the driver's negative energy.
  if (optimize_tau) {
    neg_energy <- function(t2) {
      dr <- deconv_raw(t2)$driver
      mean(pmax(-dr, 0)^2)
    }
    lo <- max(tau1 + 0.3, 1.0)
    opt <- stats::optimize(neg_energy, interval = c(lo, 8), tol = 1e-3)
    tau2 <- opt$minimum
  }

  # Step 1-2: Bateman kernel (single source: scrf) and deconvolution.
  dc <- deconv_raw(tau2)
  h <- dc$h
  h_len <- dc$h_len
  H <- dc$H
  driver_raw <- dc$driver

  # Step 3: Smooth driver with Gaussian kernel (sigma = 0.2s, width = 0.4s * sr)
  kernel_width <- max(1L, round(0.4 * sr))
  # Ensure odd kernel width
  if (kernel_width %% 2 == 0) {
    kernel_width <- kernel_width + 1L
  }
  half_k <- (kernel_width - 1L) %/% 2L
  sigma <- 0.2 * sr  # sigma in samples
  k_seq <- seq(-half_k, half_k)
  kernel <- stats::dnorm(k_seq, mean = 0, sd = sigma)
  kernel <- kernel / sum(kernel)

  driver <- as.numeric(stats::filter(driver_raw, kernel, sides = 2))
  # Fill NAs at edges with edge values
  na_idx <- which(is.na(driver))
  if (length(na_idx) > 0) {
    first_valid <- min(which(!is.na(driver)))
    last_valid <- max(which(!is.na(driver)))
    driver[na_idx[na_idx < first_valid]] <- driver[first_valid]
    driver[na_idx[na_idx > last_valid]] <- driver[last_valid]
  }

  # Step 4: Non-negativity
  driver <- pmax(driver, 0)

  # Step 5: Reconvolution: phasic = IFFT(FFT(driver) * H)
  D <- stats::fft(driver)
  phasic <- Re(stats::fft(D * H, inverse = TRUE)) / n

  # Step 6: Tonic = signal - phasic
  tonic <- signal - phasic

  list(
    tonic = tonic,
    phasic = phasic,
    driver = driver,
    tau1 = tau1,
    tau2 = tau2
  )
}


#' CDA/DDA Phasic-Driver Metrics
#'
#' Computes the Ledalab-style phasic-driver metrics from a (non-negative)
#' sudomotor driver signal: the number of significant SCRs, the sum of their
#' amplitudes, and the integrated SCR (the time-integral of the phasic driver).
#' Significant SCRs are contiguous supra-threshold driver excursions; excursions
#' separated by less than \code{min_interval_sec} are merged into one SCR.
#'
#' @param driver Numeric vector of the (non-negative) sudomotor driver. The
#'   slow tonic-driver baseline is removed (a running minimum over
#'   \code{baseline_sec}) to obtain the phasic driver before detecting SCRs.
#' @param sr Sampling rate in Hz.
#' @param amp_threshold Minimum phasic-driver amplitude for a significant SCR
#'   (default: 0.01).
#' @param min_interval_sec Minimum separation in seconds between distinct SCRs;
#'   closer excursions are merged (default: 1).
#' @param baseline_sec Window in seconds for the running-minimum tonic-driver
#'   baseline that is subtracted to form the phasic driver (default: 4).
#' @param edge_sec Seconds of phasic driver at each end to zero out before
#'   detection, suppressing deconvolution boundary artifacts (default: 1).
#' @return A list with:
#'   \describe{
#'     \item{nSCR}{Number of significant SCRs.}
#'     \item{AmpSum}{Sum of the SCR peak amplitudes (uS).}
#'     \item{ISCR}{Integrated SCR: time-integral of the phasic driver (`uS*s`).}
#'     \item{onset_sec}{Onset times of the detected SCRs (seconds).}
#'     \item{peak_sec}{Peak times of the detected SCRs (seconds).}
#'     \item{amplitudes}{Peak amplitudes of the detected SCRs (uS).}
#'     \item{scr_area}{Time-integral of the driver within each SCR (`uS*s`).}
#'   }
#' @references
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#' @seealso \code{\link{.cda_decompose}}, \code{\link{edaFeatures}}
#' @keywords internal
.cda_metrics <- function(driver, sr, amp_threshold = 0.01, min_interval_sec = 1,
                         baseline_sec = 4, edge_sec = 1) {
  drv <- pmax(as.numeric(driver), 0)
  n <- length(drv)

  # Remove the slow tonic-driver baseline (running minimum) to get the phasic
  # driver, so SCR detection is not swamped by the tonic level.
  w <- max(1L, as.integer(round(baseline_sec * sr)))
  if (w > 1L && n > w) {
    half <- w %/% 2L
    baseline <- vapply(seq_len(n), function(i)
      min(drv[max(1L, i - half):min(n, i + half)]), numeric(1))
    pd <- pmax(drv - baseline, 0)
  } else {
    pd <- drv
  }

  # Zero the phasic driver near the edges to suppress deconvolution boundary
  # artifacts (circular-convolution wrap-around) before detection and ISCR.
  em <- as.integer(round(edge_sec * sr))
  if (em > 0 && 2 * em < n) {
    pd[seq_len(em)] <- 0
    pd[(n - em + 1L):n] <- 0
  }

  # ISCR: time-integral of the phasic driver (uS*s).
  iscr <- sum(pd) / sr

  above <- pd > amp_threshold
  d <- diff(c(FALSE, above, FALSE))
  starts <- which(d == 1)
  ends <- which(d == -1) - 1L

  if (length(starts) == 0) {
    return(list(nSCR = 0L, AmpSum = 0, ISCR = iscr,
                onset_sec = numeric(0), peak_sec = numeric(0),
                amplitudes = numeric(0), scr_area = numeric(0)))
  }

  # Merge excursions closer than min_interval_sec (measured start-to-start).
  min_gap <- min_interval_sec * sr
  keep_start <- starts[1]
  merged_starts <- integer(0)
  merged_ends <- integer(0)
  cur_start <- starts[1]
  cur_end <- ends[1]
  for (i in seq_along(starts)[-1]) {
    if (starts[i] - cur_end <= min_gap) {
      cur_end <- ends[i]                    # merge
    } else {
      merged_starts <- c(merged_starts, cur_start)
      merged_ends <- c(merged_ends, cur_end)
      cur_start <- starts[i]
      cur_end <- ends[i]
    }
  }
  merged_starts <- c(merged_starts, cur_start)
  merged_ends <- c(merged_ends, cur_end)

  nscr <- length(merged_starts)
  amps <- vapply(seq_len(nscr),
                 function(i) max(pd[merged_starts[i]:merged_ends[i]]), numeric(1))
  peak_samp <- vapply(seq_len(nscr),
                      function(i) merged_starts[i] - 1L +
                        which.max(pd[merged_starts[i]:merged_ends[i]]), integer(1))
  scr_area <- vapply(seq_len(nscr),
                     function(i) sum(pd[merged_starts[i]:merged_ends[i]]) / sr,
                     numeric(1))

  list(
    nSCR = nscr,
    AmpSum = sum(amps),
    ISCR = iscr,
    onset_sec = (merged_starts - 1L) / sr,
    peak_sec = (peak_samp - 1L) / sr,
    amplitudes = amps,
    scr_area = scr_area
  )
}
