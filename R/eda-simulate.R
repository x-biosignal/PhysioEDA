#' Simulate Synthetic EDA Signals
#'
#' Generates synthetic electrodermal activity (EDA) signals with known
#' tonic (SCL) and phasic (SCR) components for testing and demonstration.
#' SCRs are modeled as biexponential impulse responses (Bateman function).
#'
#' @param n_time Number of time points (default: 6000).
#' @param n_channels Number of EDA channels (default: 1).
#' @param sr Sampling rate in Hz (default: 10).
#' @param scr_count Number of SCRs to embed (default: 5).
#' @param scl_level Baseline skin conductance level in microsiemens (default: 5.0).
#' @param scr_amplitude Mean SCR amplitude in microsiemens (default: 0.5).
#' @param noise_sd Standard deviation of Gaussian noise (default: 0.01).
#' @param seed Random seed for reproducibility (default: NULL).
#' @return A \code{\link[PhysioCore]{PhysioExperiment}} object with a single
#'   \code{"raw"} assay containing the simulated EDA signal (time x channels
#'   matrix). Channel metadata has \code{type = "EDA"} and \code{unit = "uS"}.
#'   The sampling rate is set to \code{sr}. The ground-truth components used to
#'   build the signal are stored in \code{metadata(x)$eda_truth} as a list with:
#'   \describe{
#'     \item{tonic}{time x channels matrix of the true tonic (SCL) component.}
#'     \item{phasic}{time x channels matrix of the true phasic (SCR) component.}
#'     \item{onsets}{a \code{data.frame} of the true SCR onsets with columns
#'       \code{channel}, \code{sample}, \code{time_sec}, and \code{amplitude}.}
#'   }
#'   This lets \code{edaSimulate()} serve as a ground-truth oracle for testing
#'   decomposition (e.g. \code{\link{edaDecompose}} with \code{method="cvxeda"}).
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#'
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#' @seealso \code{\link{edaDecompose}} for tonic/phasic decomposition,
#'   \code{\link{edaPeaks}} for SCR detection,
#'   \code{\link{edaFilter}} for signal filtering
#' @export
edaSimulate <- function(n_time = 6000, n_channels = 1, sr = 10,
                        scr_count = 5, scl_level = 5.0,
                        scr_amplitude = 0.5, noise_sd = 0.01,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(n_time > 0, n_channels > 0, sr > 0, scr_count >= 0)

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  tonic_mat <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  phasic_mat <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  onset_ch <- integer(0)
  onset_sample <- integer(0)
  onset_amp <- numeric(0)

  for (ch in seq_len(n_channels)) {
    # Tonic component: slow drift around scl_level
    t_sec <- seq(0, (n_time - 1) / sr, length.out = n_time)
    tonic <- scl_level + 0.3 * sin(2 * pi * 0.01 * t_sec) +
      0.1 * sin(2 * pi * 0.005 * t_sec + runif(1, 0, 2 * pi))

    # Phasic component: sum of SCR impulse responses
    phasic <- numeric(n_time)
    if (scr_count > 0) {
      # Space SCRs evenly with jitter, avoiding edges
      margin <- as.integer(n_time * 0.1)
      usable <- n_time - 2 * margin
      positions <- margin + round(seq(1, usable, length.out = scr_count))
      jitter_range <- max(1L, as.integer(usable / scr_count / 4))
      positions <- positions + sample(-jitter_range:jitter_range,
                                      scr_count, replace = TRUE)
      positions <- pmax(1L, pmin(n_time, positions))

      for (p in positions) {
        amp <- scr_amplitude * (0.5 + runif(1))
        # Record the true SCR onset (sample index and peak amplitude)
        onset_ch <- c(onset_ch, ch)
        onset_sample <- c(onset_sample, as.integer(p))
        onset_amp <- c(onset_amp, amp)
        # Bateman function: biexponential SCR shape (single kernel source: scrf)
        tau_rise <- 0.75 + runif(1) * 0.5   # 0.75-1.25 sec
        tau_decay <- 3.0 + runif(1) * 2.0   # 3-5 sec
        scr_len <- min(n_time - p + 1, as.integer(15 * sr))
        t_scr <- seq(0, by = 1 / sr, length.out = scr_len)
        scr_shape <- amp * scrf(t_scr, tau1 = tau_rise, tau2 = tau_decay,
                                form = "bateman", normalize = "none")
        scr_shape[scr_shape < 0] <- 0
        # Normalize peak to amp
        if (max(scr_shape) > 0) {
          scr_shape <- scr_shape * amp / max(scr_shape)
        }
        idx_end <- min(p + length(scr_shape) - 1, n_time)
        len <- idx_end - p + 1
        phasic[p:idx_end] <- phasic[p:idx_end] + scr_shape[seq_len(len)]
      }
    }

    # Combine components
    signal <- tonic + phasic + rnorm(n_time, sd = noise_sd)
    # EDA must be positive
    signal <- pmax(signal, 0.001)
    data[, ch] <- signal
    tonic_mat[, ch] <- tonic
    phasic_mat[, ch] <- phasic
  }

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("EDA", seq_len(n_channels)),
      type = rep("EDA", n_channels),
      unit = rep("uS", n_channels)
    ),
    samplingRate = sr
  )

  # Store ground-truth components so edaSimulate() can act as an oracle for
  # decomposition tests (true onsets, tonic, and phasic).
  onsets <- data.frame(
    channel = onset_ch,
    sample = onset_sample,
    time_sec = (onset_sample - 1L) / sr,
    amplitude = onset_amp
  )
  onsets <- onsets[order(onsets$channel, onsets$sample), , drop = FALSE]
  rownames(onsets) <- NULL
  S4Vectors::metadata(pe)$eda_truth <- list(
    tonic = tonic_mat,
    phasic = phasic_mat,
    onsets = onsets
  )
  pe
}
