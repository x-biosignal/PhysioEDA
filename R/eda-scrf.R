#' Analytic Bateman peak time
#'
#' Closed-form time of the maximum of the biexponential (Bateman) response
#' \code{h(t) = exp(-t / tau2) - exp(-t / tau1)}, obtained from the stationarity
#' condition \code{dh/dt = 0}.
#'
#' @param tau1 Rise time constant in seconds.
#' @param tau2 Decay time constant in seconds (\code{tau2 > tau1}).
#' @return The peak time in seconds.
#' @keywords internal
.scrf_peak_time <- function(tau1, tau2) {
  log(tau2 / tau1) / (1 / tau1 - 1 / tau2)
}

#' Raw biexponential (Bateman) kernel and its analytic time-derivatives
#'
#' @param x Numeric time vector (seconds).
#' @param tau1 Rise time constant in seconds.
#' @param tau2 Decay time constant in seconds.
#' @param deriv Derivative order: 0 (kernel), 1, or 2.
#' @return Numeric vector of the (\code{deriv}-th derivative of the) raw kernel.
#'   Values are not zeroed for \code{x < 0}; the caller does that.
#' @keywords internal
.scrf_bateman <- function(x, tau1, tau2, deriv = 0) {
  e1 <- exp(-x / tau1)
  e2 <- exp(-x / tau2)
  if (deriv == 0) {
    e2 - e1
  } else if (deriv == 1) {
    e1 / tau1 - e2 / tau2
  } else {
    e2 / tau2^2 - e1 / tau1^2
  }
}

#' Skin-Conductance Response Function (SCRF)
#'
#' First-class canonical skin-conductance response kernel used throughout
#' \pkg{PhysioEDA} as the single source of the SCR shape (simulation,
#' CDA / cvxEDA deconvolution, and GLM designs). The response is the
#' biexponential (Bateman) function
#' \code{h(t) = exp(-t / tau2) - exp(-t / tau1)} for \code{t >= 0} (0 otherwise),
#' optionally augmented with a Gaussian bump for the PsPM-style canonical form
#' (Bach et al., 2010). Both the kernel and its first/second time-derivatives are
#' available (the derivatives drive GLM temporal/dispersion regressors).
#'
#' @param t Numeric vector of times in seconds. Values \code{< 0} return 0.
#' @param tau1 SCR rise time constant in seconds (default: 0.75).
#' @param tau2 SCR decay time constant in seconds (default: 2.0). Must differ
#'   from \code{tau1}; the canonical shape has \code{tau2 > tau1}.
#' @param form Kernel family: \code{"bateman"} (pure biexponential) or
#'   \code{"pspm_canonical"} (biexponential plus a Gaussian bump).
#' @param normalize Normalisation of the kernel: \code{"peak"} (unit maximum,
#'   the default), \code{"area"} (unit integral), or \code{"none"} (raw).
#'   Derivatives are scaled by the same constant.
#' @param deriv Derivative order to return: 0 (kernel), 1, or 2.
#' @param bump_amp Gaussian bump amplitude, relative to the raw Bateman peak,
#'   for \code{form = "pspm_canonical"} (default: 0.3; ignored otherwise).
#' @param bump_center Gaussian bump centre in seconds (default: the analytic
#'   Bateman peak time).
#' @param bump_sd Gaussian bump standard deviation in seconds (default: 0.7).
#' @return A numeric vector the same length as \code{t}: the requested derivative
#'   of the (optionally normalised) response kernel.
#' @references
#' Bach, D.R., Flandin, G., Friston, K.J., & Dolan, R.J. (2010). "Modelling
#' event-related skin conductance responses." \emph{International Journal of
#' Psychophysiology}, 75(3), 349-356. \doi{10.1016/j.ijpsycho.2010.01.005}
#' @seealso \code{\link{scrfPeak}} for the analytic Bateman maximum,
#'   \code{\link{edaSimulate}} and \code{\link{edaDecompose}} which consume it.
#' @examples
#' t <- seq(0, 15, by = 0.1)
#' k <- scrf(t)                       # unit-peak Bateman kernel
#' dk <- scrf(t, deriv = 1)           # its temporal derivative
#' ck <- scrf(t, form = "pspm_canonical")
#' @export
scrf <- function(t, tau1 = 0.75, tau2 = 2.0,
                 form = c("bateman", "pspm_canonical"),
                 normalize = c("peak", "area", "none"),
                 deriv = 0,
                 bump_amp = 0.3, bump_center = NULL, bump_sd = 0.7) {
  form <- match.arg(form)
  normalize <- match.arg(normalize)
  stopifnot(is.numeric(t), length(tau1) == 1L, length(tau2) == 1L,
            tau1 > 0, tau2 > 0, tau1 != tau2, deriv %in% c(0, 1, 2),
            bump_sd > 0)

  tpk <- .scrf_peak_time(tau1, tau2)
  hpk <- .scrf_bateman(tpk, tau1, tau2, 0)          # analytic raw peak value
  if (is.null(bump_center)) bump_center <- tpk

  # Unit-amplitude Gaussian bump and its analytic derivatives.
  gauss <- function(x, d) {
    z <- x - bump_center
    g0 <- exp(-z^2 / (2 * bump_sd^2))
    if (d == 0) g0
    else if (d == 1) -(z / bump_sd^2) * g0
    else (z^2 / bump_sd^4 - 1 / bump_sd^2) * g0
  }

  # Raw kernel (and derivatives): Bateman, plus a scaled bump for the PsPM form.
  raw <- function(x, d) {
    v <- .scrf_bateman(x, tau1, tau2, d)
    if (form == "pspm_canonical") v <- v + bump_amp * hpk * gauss(x, d)
    v
  }

  # Normalisation constant (analytic for the Bateman form; numeric for PsPM).
  cst <- switch(
    normalize,
    none = 1,
    peak = if (form == "bateman") 1 / hpk
           else 1 / stats::optimize(function(x) raw(x, 0),
                                    interval = c(0, 10 * tau2),
                                    maximum = TRUE)$objective,
    area = if (form == "bateman") 1 / (tau2 - tau1)
           else 1 / stats::integrate(function(x) raw(x, 0), 0, 20 * tau2,
                                     rel.tol = 1e-8)$value)

  out <- cst * raw(t, deriv)
  out[t < 0] <- 0
  out
}

#' Analytic peak of the Bateman SCRF
#'
#' Returns the closed-form time and amplitude of the maximum of the raw
#' biexponential (Bateman) response used by \code{\link{scrf}}.
#'
#' @param tau1 SCR rise time constant in seconds (default: 0.75).
#' @param tau2 SCR decay time constant in seconds (default: 2.0).
#' @return A list with \code{time} (seconds) and \code{amplitude} (the raw
#'   Bateman value at the peak).
#' @seealso \code{\link{scrf}}
#' @examples
#' scrfPeak()
#' @export
scrfPeak <- function(tau1 = 0.75, tau2 = 2.0) {
  stopifnot(tau1 > 0, tau2 > 0, tau1 != tau2)
  tp <- .scrf_peak_time(tau1, tau2)
  list(time = tp, amplitude = .scrf_bateman(tp, tau1, tau2, 0))
}
