#' DDA Decomposition for a Single EDA Channel
#'
#' Performs Discrete Decomposition Analysis (DDA) on a single EDA channel
#' (Benedek & Kaernbach, 2010). Unlike Continuous Decomposition Analysis (CDA),
#' which yields a continuous driver, DDA represents the phasic activity as a
#' small set of \emph{discrete} SCRs: the continuous CDA driver is reduced to a
#' sparse impulse train (one impulse per significant SCR, placed at its onset and
#' weighted by the integrated driver of that SCR), which is reconvolved with the
#' Bateman response to obtain the phasic component.
#'
#' @param signal Numeric vector of EDA values.
#' @param sr Sampling rate in Hz.
#' @param tau1 Rise time constant in seconds (default: 0.75).
#' @param tau2 Decay time constant in seconds (default: 2.0).
#' @param optimize_tau Logical; passed to \code{\link{.cda_decompose}} to
#'   optimise \code{tau2} before discretisation (default: \code{FALSE}).
#' @param amp_threshold Minimum driver amplitude for a significant SCR
#'   (default: 0.01).
#' @param min_interval_sec Minimum separation in seconds between distinct SCRs
#'   (default: 1).
#' @return A list with components:
#'   \describe{
#'     \item{tonic}{Numeric vector of the tonic (SCL) component.}
#'     \item{phasic}{Numeric vector of the phasic (SCR) component from the
#'       discrete driver.}
#'     \item{driver}{Numeric vector of the discrete (impulse-train) driver.}
#'     \item{tau1}{Rise time constant used.}
#'     \item{tau2}{Decay time constant used.}
#'     \item{onsets}{Onset times of the discrete SCRs (seconds).}
#'     \item{nSCR}{Number of discrete SCRs.}
#'   }
#' @references
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#' @seealso \code{\link{.cda_decompose}} for the continuous (CDA) variant,
#'   \code{\link{edaDecompose}} for the user-facing wrapper.
#' @keywords internal
.dda_decompose <- function(signal, sr, tau1 = 0.75, tau2 = 2.0,
                           optimize_tau = FALSE, amp_threshold = 0.01,
                           min_interval_sec = 1) {
  stopifnot(is.numeric(signal), length(signal) > 1, is.numeric(sr), sr > 0)
  n <- length(signal)

  # Continuous decomposition first, then discretise its driver.
  cda <- .cda_decompose(signal, sr, tau1 = tau1, tau2 = tau2,
                        optimize_tau = optimize_tau)
  tau2 <- cda$tau2
  m <- .cda_metrics(cda$driver, sr, amp_threshold = amp_threshold,
                    min_interval_sec = min_interval_sec)

  # Discrete driver: an impulse at each SCR onset, weighted so that its
  # reconvolution matches the integrated CDA driver of that SCR
  # (reconvolving an impulse of weight w with a unit-peak Bateman gives an SCR
  # of peak w; the CDA driver of the SCR integrates to scr_area, i.e. weight
  # scr_area * sr in sample units).
  disc <- numeric(n)
  if (m$nSCR > 0) {
    onset_samples <- as.integer(round(m$onset_sec * sr)) + 1L
    onset_samples <- pmax(1L, pmin(n, onset_samples))
    disc[onset_samples] <- disc[onset_samples] + m$scr_area * sr
  }

  # Reconvolve the discrete driver with the Bateman kernel (single source: scrf).
  h_len <- min(n, as.integer(10 * tau2 * sr))
  t_h <- seq(0, by = 1 / sr, length.out = h_len)
  h <- scrf(t_h, tau1 = tau1, tau2 = tau2, form = "bateman", normalize = "none")
  h <- h / max(h)
  H <- stats::fft(c(h, rep(0, n - h_len)))
  phasic <- Re(stats::fft(stats::fft(disc) * H, inverse = TRUE)) / n
  phasic <- pmax(phasic, 0)
  tonic <- signal - phasic

  list(
    tonic = tonic,
    phasic = phasic,
    driver = disc,
    tau1 = tau1,
    tau2 = tau2,
    onsets = m$onset_sec,
    nSCR = m$nSCR
  )
}
