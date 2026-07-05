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
#'   The sampling rate is set to \code{sr}.
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
        # Bateman function: biexponential SCR shape
        tau_rise <- 0.75 + runif(1) * 0.5   # 0.75-1.25 sec
        tau_decay <- 3.0 + runif(1) * 2.0   # 3-5 sec
        scr_len <- min(n_time - p + 1, as.integer(15 * sr))
        t_scr <- seq(0, by = 1 / sr, length.out = scr_len)
        scr_shape <- amp * (exp(-t_scr / tau_decay) - exp(-t_scr / tau_rise))
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
  }

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("EDA", seq_len(n_channels)),
      type = rep("EDA", n_channels),
      unit = rep("uS", n_channels)
    ),
    samplingRate = sr
  )
}
