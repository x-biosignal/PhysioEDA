#' Dynamic Causal Model for Skin-Conductance Responses
#'
#' Estimates trial-by-trial latent sudomotor input amplitudes that best explain
#' the skin-conductance signal, in the spirit of the PsPM dynamic causal model
#' for SCR (Bach et al., 2010, 2011). Each event contributes a non-negative
#' latent input whose skin-conductance response is the canonical
#' \code{\link{scrf}} kernel; the per-event amplitudes are estimated by
#' non-negative least squares against a low-frequency cosine drift basis, and
#' the input timings are optionally refined by a nonlinear local search
#' (the "dynamic" step).
#'
#' Unlike \code{\link{edaGLM}}, which estimates one amplitude per condition,
#' \code{edaDCM} estimates one latent input amplitude per event.
#'
#' @param x A PhysioExperiment object with EDA data.
#' @param events A \code{data.frame} with numeric \code{onset} (seconds); a
#'   \code{condition} column is optional and carried through.
#' @param channel EDA channel index or label (default: 1).
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @param tau1,tau2 SCR rise/decay time constants in seconds.
#' @param drift Number of cosine drift terms (default: 6).
#' @param refine_latency Logical; refine each event's onset by a local search to
#'   improve the fit (default: \code{TRUE}).
#' @param latency_window Half-width in seconds of the latency search
#'   (default: 1.5).
#' @return A list of class \code{"eda_dcm"} with:
#'   \describe{
#'     \item{estimates}{a \code{data.frame} with \code{event}, \code{onset},
#'       \code{onset_refined}, \code{amplitude} (per-event latent input), and
#'       \code{condition} if present.}
#'     \item{fitted}{fitted signal; \code{residuals}; \code{r_squared}.}
#'   }
#' @references
#' Bach, D.R., Daunizeau, J., Friston, K.J., & Dolan, R.J. (2010). "Dynamic
#' causal modelling of anticipatory skin conductance responses." \emph{Biological
#' Psychology}, 85(1), 163-170. \doi{10.1016/j.biopsycho.2010.06.007}
#' @seealso \code{\link{edaGLM}} for the condition-level GLM,
#'   \code{\link{scrf}} for the response kernel.
#' @examples
#' ev <- data.frame(onset = c(10, 30, 50, 70))
#' pe <- edaSimulate(n_time = 1000, sr = 10, scr_count = 4, seed = 1)
#' fit <- edaDCM(pe, ev, refine_latency = FALSE)
#' fit$estimates
#' @export
edaDCM <- function(x, events, channel = 1, assay_name = NULL,
                   tau1 = 0.75, tau2 = 2.0, drift = 6,
                   refine_latency = TRUE, latency_window = 1.5) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.data.frame(events), "onset" %in% names(events),
            is.numeric(events$onset), nrow(events) >= 1)
  if (is.null(assay_name)) assay_name <- defaultAssay(x)

  ch <- .eda_resolve_channel(x, channel)
  y <- SummarizedExperiment::assay(x, assay_name)[, ch]
  n <- length(y)
  sr <- samplingRate(x)
  ne <- nrow(events)
  kern <- .scrf_kernel(sr, tau1, tau2, deriv = 0)
  Xdr <- .drift_basis(n, drift)
  nn <- ncol(Xdr)

  has_nnls <- requireNamespace("quadprog", quietly = TRUE)

  # Fit per-event amplitudes (non-negative) plus drift for given onsets.
  fit_amps <- function(ons) {
    Xev <- vapply(ons, function(o) .scrf_regressor(o, n, sr, kern), numeric(n))
    Xf <- cbind(Xev, Xdr)
    if (has_nnls) {
      D <- crossprod(Xf)
      D <- D + diag(1e-8 * mean(diag(D)) + 1e-10, ncol(D))
      d <- crossprod(Xf, y)
      Amat <- rbind(diag(ne), matrix(0, nn, ne))    # event amps >= 0
      sol <- tryCatch(
        quadprog::solve.QP(D, d, Amat, rep(0, ne), meq = 0)$solution,
        error = function(e) NULL)
      if (!is.null(sol)) {
        return(list(amps = sol[seq_len(ne)], nuis = sol[ne + seq_len(nn)]))
      }
    }
    # Fallback: OLS then clamp negatives to zero.
    b <- stats::.lm.fit(Xf, y)$coefficients
    list(amps = pmax(b[seq_len(ne)], 0), nuis = b[ne + seq_len(nn)])
  }

  ons <- events$onset
  fitres <- fit_amps(ons)
  amps <- fitres$amps
  nuis <- fitres$nuis

  if (refine_latency && ne >= 1) {
    shifts <- seq(-latency_window, latency_window, by = 1 / sr * max(1, round(sr / 4)))
    for (iter in seq_len(3)) {
      base <- as.numeric(Xdr %*% nuis)
      full <- base + rowSums(vapply(seq_len(ne),
        function(k) amps[k] * .scrf_regressor(ons[k], n, sr, kern), numeric(n)))
      for (k in seq_len(ne)) {
        contrib_k <- amps[k] * .scrf_regressor(ons[k], n, sr, kern)
        pr <- y - full + contrib_k                  # partial residual for event k
        best <- ons[k]; bestv <- Inf
        for (sh in shifts) {
          o2 <- ons[k] + sh
          v <- sum((pr - amps[k] * .scrf_regressor(o2, n, sr, kern))^2)
          if (v < bestv) { bestv <- v; best <- o2 }
        }
        full <- full - contrib_k + amps[k] * .scrf_regressor(best, n, sr, kern)
        ons[k] <- best
      }
      fitres <- fit_amps(ons)
      amps <- fitres$amps
      nuis <- fitres$nuis
    }
  }

  Xev <- vapply(ons, function(o) .scrf_regressor(o, n, sr, kern), numeric(n))
  fitted <- as.numeric(cbind(Xev, Xdr) %*% c(amps, nuis))
  resid <- y - fitted
  r2 <- 1 - sum(resid^2) / sum((y - mean(y))^2)

  est <- data.frame(
    event = seq_len(ne),
    onset = events$onset,
    onset_refined = ons,
    amplitude = amps,
    row.names = NULL, stringsAsFactors = FALSE
  )
  if ("condition" %in% names(events)) est$condition <- as.character(events$condition)

  structure(list(
    estimates = est,
    fitted = fitted,
    residuals = resid,
    r_squared = r2,
    refined = refine_latency
  ), class = "eda_dcm")
}
