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
#' @return A list with components:
#'   \describe{
#'     \item{tonic}{Numeric vector of the tonic (SCL) component.}
#'     \item{phasic}{Numeric vector of the phasic (SCR) component.}
#'     \item{driver}{Numeric vector of the sudomotor nerve activity driver signal.}
#'     \item{tau1}{Rise time constant used.}
#'     \item{tau2}{Decay time constant used.}
#'   }
#' @references
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#' @seealso \code{\link{edaDecompose}} for the user-facing decomposition wrapper,
#'   \code{\link{.cvxeda_decompose}} for the cvxEDA alternative
#' @keywords internal
.cda_decompose <- function(signal, sr, tau1 = 0.75, tau2 = 2.0) {
  stopifnot(is.numeric(signal), length(signal) > 1)
  stopifnot(is.numeric(sr), sr > 0)
  stopifnot(is.numeric(tau1), tau1 > 0)
  stopifnot(is.numeric(tau2), tau2 > 0, tau2 > tau1)

  n <- length(signal)

  # Step 1: Build Bateman impulse response h(t) = exp(-t/tau2) - exp(-t/tau1)
  h_len <- min(n, as.integer(10 * tau2 * sr))
  t_h <- seq(0, by = 1 / sr, length.out = h_len)
  h <- exp(-t_h / tau2) - exp(-t_h / tau1)
  # Normalize so max(h) = 1
  h_max <- max(h)
  if (h_max > 0) {
    h <- h / h_max
  }

  # Step 2: Deconvolution via FFT
  # Zero-pad h to same length as signal
  h_padded <- c(h, rep(0, n - h_len))
  H <- stats::fft(h_padded)
  S <- stats::fft(signal)

  # Regularization: epsilon = max(|H|) * 1e-6
  epsilon <- max(Mod(H)) * 1e-6
  driver_raw <- Re(stats::fft(S / (H + epsilon), inverse = TRUE)) / n

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
