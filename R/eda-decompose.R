#' Decompose EDA into Tonic and Phasic Components
#'
#' Separates an electrodermal activity (EDA) signal into its slow-varying
#' tonic component (skin conductance level, SCL) and fast-varying phasic
#' component (skin conductance responses, SCR).
#'
#' Five methods are available:
#' \describe{
#'   \item{highpass}{FFT-based highpass/lowpass separation at a cutoff frequency.}
#'   \item{median}{Sliding median filter for tonic extraction.}
#'   \item{cda}{Continuous Decomposition Analysis (Benedek & Kaernbach, 2010).
#'     Deconvolution with a Bateman impulse response, Gaussian smoothing, and a
#'     non-negativity constraint on the driver signal; optional \code{tau2}
#'     optimisation.}
#'   \item{dda}{Discrete Decomposition Analysis (Benedek & Kaernbach, 2010): the
#'     continuous CDA driver is reduced to a sparse impulse train (one impulse
#'     per significant SCR) and reconvolved.}
#'   \item{cvxeda}{Exact convex optimization decomposition (Greco et al., 2016),
#'     solved as a quadratic program: a cubic-knot spline tonic basis and a
#'     sparse, non-negative sudomotor driver through a biexponential system, with
#'     L1 sparsity on the driver and L2 smoothness on the tonic. Requires the
#'     \pkg{quadprog} solver; degrades gracefully to the Wiener approximation if
#'     it is unavailable.}
#'   \item{cvxeda_fast}{Fast Wiener/ADMM approximation of cvxEDA (no solver
#'     dependency); the original iterative approximation.}
#' }
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param method Decomposition method: \code{"highpass"}, \code{"median"},
#'   \code{"cda"}, \code{"dda"}, \code{"cvxeda"}, or \code{"cvxeda_fast"}. Default
#'   is \code{"highpass"}.
#' @param cutoff Cutoff frequency in Hz for the highpass method (default: 0.05).
#' @param window_sec Window length in seconds for the median method (default: 4).
#' @param tau1 SCR rise time constant in seconds for CDA/cvxEDA (default: 0.75).
#' @param tau2 SCR decay time constant in seconds for CDA/cvxEDA (default: 2.0).
#' @param alpha L1 sparsity penalty for cvxEDA (default: 0.01). Ignored by
#'   other methods.
#' @param gamma Smoothness weight for cvxEDA tonic component (default: 0.1).
#'   Ignored by other methods.
#' @param delta_knot Tonic spline knot spacing in seconds for the \code{"cvxeda"}
#'   QP (default: 10).
#' @param qp_sr Target rate (Hz) at which the \code{"cvxeda"} QP is solved; the
#'   signal is decimated towards this rate for tractability (default: 4).
#' @param optimize_tau Logical; for \code{"cda"} / \code{"dda"}, optimise the
#'   decay constant \code{tau2} to minimise the driver's negative energy
#'   (default: \code{FALSE}).
#' @param assay_name Name of the input assay. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param output_tonic Name for the tonic output assay (default: "tonic").
#' @param output_phasic Name for the phasic output assay (default: "phasic").
#' @return A modified \code{\link[PhysioCore]{PhysioExperiment}} with new assays:
#'   \describe{
#'     \item{tonic}{The slow-varying skin conductance level (SCL) component.}
#'     \item{phasic}{The fast-varying skin conductance response (SCR) component.}
#'     \item{driver}{(CDA and cvxEDA variants only) The sudomotor nerve activity
#'       driver signal.}
#'   }
#'   Decomposition parameters are stored in \code{metadata(x)$eda_decompose};
#'   for \code{method = "cvxeda"} this also records the per-channel solver
#'   \code{engine} ("qp" or "wiener_fallback") and final QP \code{objective}.
#' @references
#' Benedek, M., & Kaernbach, C. (2010). "A continuous measure of phasic
#' electrodermal activity." \emph{Journal of Neuroscience Methods}, 190(1),
#' 80-91. \doi{10.1016/j.jneumeth.2010.04.028}
#'
#' Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
#' electrodermal activity processing." \emph{IEEE Transactions on Biomedical
#' Engineering}, 63(4), 797-804. \doi{10.1109/TBME.2015.2474131}
#' @seealso \code{\link{edaPeaks}} for SCR peak detection on the phasic signal,
#'   \code{\link{edaFeatures}} for feature extraction,
#'   \code{\link{plotDecompose}} for visualizing decomposition results
#' @export
edaDecompose <- function(x, method = c("highpass", "median", "cda", "dda",
                                       "cvxeda", "cvxeda_fast"),
                         cutoff = 0.05, window_sec = 4,
                         tau1 = 0.75, tau2 = 2.0,
                         alpha = 0.01, gamma = 0.1,
                         delta_knot = 10, qp_sr = 4,
                         optimize_tau = FALSE,
                         assay_name = NULL,
                         output_tonic = "tonic",
                         output_phasic = "phasic") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  tonic_mat <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  phasic_mat <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  driver_mat <- NULL

  if (method %in% c("cda", "dda", "cvxeda", "cvxeda_fast")) {
    driver_mat <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  }
  cvx_info <- NULL
  cda_info <- NULL

  if (method == "median") {
    # Sliding median filter for tonic extraction
    window_samples <- as.integer(sr * window_sec)
    # Ensure odd window size
    if (window_samples %% 2 == 0) {
      window_samples <- window_samples + 1L
    }
    half_win <- (window_samples - 1L) %/% 2L

    for (ch in seq_len(n_channels)) {
      signal <- data[, ch]
      tonic_ch <- numeric(n_time)
      for (i in seq_len(n_time)) {
        lo <- max(1L, i - half_win)
        hi <- min(n_time, i + half_win)
        tonic_ch[i] <- median(signal[lo:hi])
      }
      tonic_mat[, ch] <- tonic_ch
      phasic_mat[, ch] <- signal - tonic_ch
    }

  } else if (method == "cda") {
    # Continuous Decomposition Analysis (Benedek & Kaernbach, 2010)
    taus <- matrix(NA_real_, nrow = n_channels, ncol = 2,
                   dimnames = list(NULL, c("tau1", "tau2")))
    metrics <- vector("list", n_channels)
    for (ch in seq_len(n_channels)) {
      result <- .cda_decompose(data[, ch], sr, tau1 = tau1, tau2 = tau2,
                               optimize_tau = optimize_tau)
      tonic_mat[, ch] <- result$tonic
      phasic_mat[, ch] <- result$phasic
      driver_mat[, ch] <- result$driver
      taus[ch, ] <- c(result$tau1, result$tau2)
      metrics[[ch]] <- .cda_metrics(result$driver, sr)
    }
    cda_info <- list(method = "cda", tau = taus, metrics = metrics)

  } else if (method == "dda") {
    # Discrete Decomposition Analysis (Benedek & Kaernbach, 2010)
    taus <- matrix(NA_real_, nrow = n_channels, ncol = 2,
                   dimnames = list(NULL, c("tau1", "tau2")))
    metrics <- vector("list", n_channels)
    for (ch in seq_len(n_channels)) {
      result <- .dda_decompose(data[, ch], sr, tau1 = tau1, tau2 = tau2,
                               optimize_tau = optimize_tau)
      tonic_mat[, ch] <- result$tonic
      phasic_mat[, ch] <- result$phasic
      driver_mat[, ch] <- result$driver
      taus[ch, ] <- c(result$tau1, result$tau2)
      metrics[[ch]] <- .cda_metrics(result$driver, sr)
    }
    cda_info <- list(method = "dda", tau = taus, metrics = metrics)

  } else if (method == "cvxeda") {
    # Exact cvxEDA quadratic program (Greco et al., 2016)
    engines <- character(n_channels)
    objectives <- rep(NA_real_, n_channels)
    for (ch in seq_len(n_channels)) {
      result <- .cvxeda_qp(data[, ch], sr, tau1 = tau1, tau2 = tau2,
                           alpha = alpha, gamma = gamma,
                           delta_knot = delta_knot, qp_sr = qp_sr)
      tonic_mat[, ch] <- result$tonic
      phasic_mat[, ch] <- result$phasic
      driver_mat[, ch] <- result$driver
      engines[ch] <- result$engine
      objectives[ch] <- result$objective
    }
    cvx_info <- list(engine = engines, objective = objectives)

  } else if (method == "cvxeda_fast") {
    # Fast Wiener/ADMM approximation of cvxEDA
    for (ch in seq_len(n_channels)) {
      result <- .cvxeda_decompose(data[, ch], sr, tau1 = tau1, tau2 = tau2,
                                  alpha = alpha, gamma = gamma)
      tonic_mat[, ch] <- result$tonic
      phasic_mat[, ch] <- result$phasic
      driver_mat[, ch] <- result$driver
    }

  } else {
    # FFT-based highpass/lowpass decomposition
    freqs <- seq(0, sr / 2, length.out = floor(n_time / 2) + 1)

    for (ch in seq_len(n_channels)) {
      signal <- data[, ch]
      ft <- fft(signal)

      # Build frequency vector for all FFT bins
      freq_bins <- (seq_len(n_time) - 1) * sr / n_time
      # Mirror frequencies above Nyquist
      freq_bins[freq_bins > sr / 2] <- sr - freq_bins[freq_bins > sr / 2]

      # Smooth highpass mask: quadratic transition from 0 to cutoff
      mask <- numeric(n_time)
      for (i in seq_len(n_time)) {
        f <- freq_bins[i]
        if (f >= cutoff) {
          mask[i] <- 1.0
        } else {
          # Smooth quadratic transition
          mask[i] <- (f / cutoff)^2
        }
      }

      phasic_ft <- ft * mask
      tonic_ft <- ft * (1 - mask)

      phasic_mat[, ch] <- Re(fft(phasic_ft, inverse = TRUE)) / n_time
      tonic_mat[, ch] <- Re(fft(tonic_ft, inverse = TRUE)) / n_time
    }
  }

  # Store results as new assays
  current_assays <- SummarizedExperiment::assays(x)
  current_assays[[output_tonic]] <- tonic_mat
  current_assays[[output_phasic]] <- phasic_mat
  if (!is.null(driver_mat)) {
    current_assays[["driver"]] <- driver_mat
  }
  SummarizedExperiment::assays(x, withDimnames = FALSE) <- current_assays

  # Store decomposition parameters in metadata
  meta <- S4Vectors::metadata(x)
  meta$eda_decompose <- list(
    method = method,
    cutoff = cutoff,
    window_sec = window_sec,
    tau1 = tau1,
    tau2 = tau2,
    alpha = alpha,
    gamma = gamma,
    delta_knot = delta_knot,
    qp_sr = qp_sr,
    optimize_tau = optimize_tau,
    cvxeda = cvx_info,
    cda = cda_info
  )
  S4Vectors::metadata(x) <- meta

  .recordProv(x, input_assay = assay_name, output_assay = output_phasic,
              .package = "PhysioEDA")
}
