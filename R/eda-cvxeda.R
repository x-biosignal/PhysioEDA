#' cvxEDA Decomposition for a Single Channel (Simplified)
#'
#' Performs a simplified convex optimization-based EDA decomposition using
#' an iterative ADMM-like approach with FFT-based convolution/deconvolution.
#' Based on Greco et al. (2016), but avoids the CVXR dependency by using
#' an iterative proximal gradient method with Wiener deconvolution.
#'
#' @param signal Numeric vector of EDA signal values.
#' @param sr Sampling rate in Hz.
#' @param tau1 SCR rise time constant in seconds (default: 0.75).
#' @param tau2 SCR decay time constant in seconds (default: 2.0).
#' @param alpha L1 sparsity penalty on the driver signal (default: 0.01).
#' @param gamma Smoothness weight for tonic component (default: 0.1).
#' @param max_iter Maximum number of ADMM iterations (default: 50).
#' @param tol Convergence tolerance for the driver signal (default: 1e-4).
#' @return A list with components:
#'   \describe{
#'     \item{tonic}{Numeric vector of tonic (SCL) component.}
#'     \item{phasic}{Numeric vector of phasic (SCR) component.}
#'     \item{driver}{Numeric vector of the sparse driver signal.}
#'     \item{alpha}{The alpha parameter used.}
#'     \item{gamma}{The gamma parameter used.}
#'     \item{iterations}{Number of iterations performed.}
#'     \item{converged}{Logical indicating whether the algorithm converged.}
#'   }
#' @references
#' Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
#' electrodermal activity processing." \emph{IEEE Transactions on Biomedical
#' Engineering}, 63(4), 797-804. \doi{10.1109/TBME.2015.2474131}
#' @seealso \code{\link{edaDecompose}} for the user-facing decomposition wrapper,
#'   \code{\link{.cda_decompose}} for the CDA alternative
#' @keywords internal
.cvxeda_decompose <- function(signal, sr, tau1 = 0.75, tau2 = 2.0,
                               alpha = 0.01, gamma = 0.1,
                               max_iter = 50, tol = 1e-4) {
  n <- length(signal)
  stopifnot(n > 0, sr > 0, tau1 > 0, tau2 > 0, alpha >= 0, gamma >= 0)

  # --- Step 1: Build Bateman impulse response ---
  h_len <- min(n, as.integer(15 * sr))
  t_h <- seq(0, by = 1 / sr, length.out = h_len)
  h <- exp(-t_h / tau2) - exp(-t_h / tau1)
  h[h < 0] <- 0
  # Normalize to unit peak for interpretable driver amplitudes
  if (max(h) > 0) {
    h <- h / max(h)
  }

  # --- Step 2: Build FFT-based convolution/deconvolution operators ---
  h_padded <- c(h, rep(0, n - h_len))
  H_fft <- stats::fft(h_padded)
  H_conj <- Conj(H_fft)
  H_mag2 <- Re(H_fft * H_conj)  # |H(f)|^2

  # Wiener deconvolution filter: H*(f) / (|H(f)|^2 + lambda)
  # lambda controls regularization strength
  lambda <- 0.01 * max(H_mag2)
  W_fft <- H_conj / (H_mag2 + lambda)

  # --- Step 3: Initialize ---
  driver <- rep(0, n)
  tonic <- rep(mean(signal), n)

  # Precompute smoothing window for tonic
  smooth_window <- max(3L, as.integer(gamma * sr * 10))
  if (smooth_window %% 2 == 0) {
    smooth_window <- smooth_window + 1L
  }

  converged <- FALSE
  iterations <- max_iter

  # --- Step 4: Iterative ADMM-like updates ---
  for (iter in seq_len(max_iter)) {
    driver_prev <- driver

    # (a) Compute residual: signal minus current tonic estimate
    r <- signal - tonic

    # (b) Wiener deconvolution of residual to get driver estimate
    R_fft <- stats::fft(r)
    driver_raw <- Re(stats::fft(R_fft * W_fft, inverse = TRUE)) / n

    # (c) Proximal operator: non-negativity + L1 sparsity
    driver <- pmax(0, driver_raw - alpha)

    # (d) Reconvolve driver to get phasic component
    D_fft <- stats::fft(driver)
    phasic <- Re(stats::fft(D_fft * H_fft, inverse = TRUE)) / n

    # (e) Update tonic with smoothness constraint (running mean)
    tonic_raw <- signal - phasic
    tonic <- .running_mean(tonic_raw, smooth_window)

    # (f) Check convergence
    driver_max <- max(abs(driver))
    if (driver_max < 1e-10) {
      converged <- TRUE
      iterations <- iter
      break
    }
    rel_change <- max(abs(driver - driver_prev)) / (driver_max + 1e-10)
    if (rel_change < tol) {
      converged <- TRUE
      iterations <- iter
      break
    }
  }

  # --- Step 5: Final reconstruction ---
  # Reconvolve final driver to get phasic, derive tonic as residual
  D_fft <- stats::fft(driver)
  phasic <- Re(stats::fft(D_fft * H_fft, inverse = TRUE)) / n
  tonic <- signal - phasic

  list(
    tonic = tonic,
    phasic = phasic,
    driver = driver,
    alpha = alpha,
    gamma = gamma,
    iterations = iterations,
    converged = converged
  )
}


#' Running Mean Smoother
#'
#' Computes a symmetric running mean with edge handling via partial windows.
#'
#' @param x Numeric vector.
#' @param window Integer window size (must be odd).
#' @return Smoothed numeric vector of same length.
#' @keywords internal
.running_mean <- function(x, window) {
  n <- length(x)
  if (window >= n) {
    return(rep(mean(x), n))
  }
  half <- (window - 1L) %/% 2L
  result <- numeric(n)
  # Use cumsum for efficient running mean
  cs <- c(0, cumsum(x))
  for (i in seq_len(n)) {
    lo <- max(1L, i - half)
    hi <- min(n, i + half)
    result[i] <- (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
  }
  result
}
