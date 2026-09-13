#' Sampled SCRF Basis Kernel
#'
#' Samples the skin-conductance response function (\code{\link{scrf}}) on the
#' signal grid for use as a GLM/DCM basis regressor kernel.
#'
#' @param sr Sampling rate in Hz.
#' @param tau1 SCR rise time constant in seconds.
#' @param tau2 SCR decay time constant in seconds.
#' @param deriv Derivative order (0 canonical, 1 temporal derivative).
#' @param len_sec Kernel length in seconds (default: \code{10 * tau2}).
#' @return Numeric kernel vector, scaled to unit peak absolute value.
#' @keywords internal
.scrf_kernel <- function(sr, tau1 = 0.75, tau2 = 2.0, deriv = 0,
                         len_sec = NULL) {
  if (is.null(len_sec)) len_sec <- 10 * tau2
  t <- seq(0, len_sec, by = 1 / sr)
  k <- scrf(t, tau1 = tau1, tau2 = tau2, form = "bateman",
            normalize = "none", deriv = deriv)
  mx <- max(abs(k))
  if (mx > 0) k <- k / mx
  k
}

#' Build an Event Regressor
#'
#' Places the basis \code{kern} at each onset and sums (a stick function
#' convolved with the SCRF).
#'
#' @param onsets Numeric event onsets in seconds.
#' @param n Number of samples.
#' @param sr Sampling rate in Hz.
#' @param kern Basis kernel from \code{\link{.scrf_kernel}}.
#' @return Numeric regressor vector of length \code{n}.
#' @keywords internal
.scrf_regressor <- function(onsets, n, sr, kern) {
  reg <- numeric(n)
  kl0 <- length(kern)
  for (o in onsets) {
    s <- as.integer(round(o * sr)) + 1L
    if (s < 1L || s > n) next
    kl <- min(kl0, n - s + 1L)
    reg[s:(s + kl - 1L)] <- reg[s:(s + kl - 1L)] + kern[seq_len(kl)]
  }
  reg
}

#' Low-frequency Drift (High-pass) Nuisance Basis
#'
#' A discrete-cosine drift basis (intercept plus low-frequency cosines) that
#' absorbs the slow tonic level, equivalent to a high-pass filter in the GLM.
#'
#' @param n Number of samples.
#' @param k Number of cosine terms (default: 6).
#' @return An \code{n x (k + 1)} matrix.
#' @keywords internal
.drift_basis <- function(n, k = 6) {
  tt <- seq(0, 1, length.out = n)
  B <- matrix(1, n, 1)
  if (k > 0) {
    for (j in seq_len(k)) B <- cbind(B, cos(pi * j * tt))
  }
  colnames(B) <- c("drift0", if (k > 0) paste0("drift", seq_len(k)))
  B
}

#' PsPM-style GLM for Skin-Conductance Responses
#'
#' Fits a general linear model of the skin-conductance signal against
#' event-related regressors, each a stick function at the event onsets convolved
#' with the canonical skin-conductance response function (\code{\link{scrf}};
#' Bach & Friston, 2013). Events are grouped into conditions; each condition
#' yields an amplitude estimate (beta) with a standard error and t statistic.
#' A low-frequency cosine drift basis absorbs the tonic level. Optionally, a
#' temporal-derivative basis is added per condition to accommodate response
#' latency differences.
#'
#' @param x A PhysioExperiment object with EDA data.
#' @param events A \code{data.frame} with numeric \code{onset} (seconds) and a
#'   \code{condition} column grouping events.
#' @param conditions Optional character vector selecting/ordering the conditions
#'   (default: the unique conditions in \code{events}).
#' @param channel EDA channel index or label to model (default: 1).
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @param tau1,tau2 SCR rise/decay time constants in seconds.
#' @param basis \code{"canonical"} (SCRF only) or \code{"derivative"} (SCRF plus
#'   its temporal derivative per condition).
#' @param drift Number of cosine drift terms (default: 6).
#' @return A list of class \code{"eda_glm"} with:
#'   \describe{
#'     \item{coefficients}{a \code{data.frame} with \code{condition},
#'       \code{beta}, \code{se}, \code{t}, \code{p} (canonical amplitude per
#'       condition; \code{NA} for a condition dropped due to collinearity).}
#'     \item{derivative}{(basis = "derivative") the derivative-regressor betas.}
#'     \item{fitted}{fitted signal; \code{residuals}; \code{sigma}; \code{df};
#'       \code{r_squared}; \code{rank}; \code{n_regressors}; \code{collinear}.}
#'   }
#' @references
#' Bach, D.R., & Friston, K.J. (2013). "Model-based analysis of skin conductance
#' responses." \emph{Psychophysiology}, 50(1), 15-22. \doi{10.1111/j.1469-8986.2012.01483.x}
#' @seealso \code{\link{edaDCM}} for trial-wise latent input estimation,
#'   \code{\link{scrf}} for the response kernel.
#' @examples
#' ev <- data.frame(onset = c(10, 30, 50, 70), condition = c("A", "A", "B", "B"))
#' pe <- edaSimulate(n_time = 1000, sr = 10, scr_count = 4, seed = 1)
#' fit <- edaGLM(pe, ev)
#' fit$coefficients
#' @export
edaGLM <- function(x, events, conditions = NULL, channel = 1,
                   assay_name = NULL, tau1 = 0.75, tau2 = 2.0,
                   basis = c("canonical", "derivative"), drift = 6) {
  stopifnot(inherits(x, "PhysioExperiment"))
  basis <- match.arg(basis)
  stopifnot(is.data.frame(events), all(c("onset", "condition") %in% names(events)),
            is.numeric(events$onset))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)

  ch <- .eda_resolve_channel(x, channel)
  y <- SummarizedExperiment::assay(x, assay_name)[, ch]
  n <- length(y)
  sr <- samplingRate(x)

  if (is.null(conditions)) conditions <- unique(as.character(events$condition))
  kern <- .scrf_kernel(sr, tau1, tau2, deriv = 0)
  kern_d <- if (basis == "derivative") .scrf_kernel(sr, tau1, tau2, deriv = 1) else NULL

  # Condition regressors (+ derivative), named.
  Xc <- vapply(conditions, function(c)
    .scrf_regressor(events$onset[as.character(events$condition) == c], n, sr, kern),
    numeric(n))
  colnames(Xc) <- conditions
  Xd <- NULL
  if (basis == "derivative") {
    Xd <- vapply(conditions, function(c)
      .scrf_regressor(events$onset[as.character(events$condition) == c], n, sr, kern_d),
      numeric(n))
    colnames(Xd) <- paste0(conditions, "_d")
  }
  Xdr <- .drift_basis(n, drift)

  design <- cbind(Xc, Xd, Xdr)
  rank <- qr(design)$rank
  collinear <- rank < ncol(design)
  if (collinear) {
    warning(sprintf(paste0("edaGLM: design matrix is rank-deficient (rank %d of ",
                           "%d columns); collinear regressors are dropped."),
                    rank, ncol(design)), call. = FALSE)
  }

  dm <- as.data.frame(design)
  dm$.y <- y
  fit <- stats::lm(.y ~ . - 1, data = dm)
  sm <- stats::summary.lm(fit)
  co <- sm$coefficients

  # per-condition canonical amplitude table (NA if dropped)
  get_row <- function(nm) if (nm %in% rownames(co)) co[nm, ] else rep(NA_real_, 4)
  coef_df <- data.frame(
    condition = conditions,
    beta = vapply(conditions, function(c) get_row(c)[1], numeric(1)),
    se = vapply(conditions, function(c) get_row(c)[2], numeric(1)),
    t = vapply(conditions, function(c) get_row(c)[3], numeric(1)),
    p = vapply(conditions, function(c) get_row(c)[4], numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )

  deriv_df <- NULL
  if (basis == "derivative") {
    deriv_df <- data.frame(
      condition = conditions,
      beta = vapply(paste0(conditions, "_d"), function(c) get_row(c)[1], numeric(1)),
      se = vapply(paste0(conditions, "_d"), function(c) get_row(c)[2], numeric(1)),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }

  structure(list(
    coefficients = coef_df,
    derivative = deriv_df,
    fitted = stats::fitted(fit),
    residuals = stats::residuals(fit),
    sigma = sm$sigma,
    df = fit$df.residual,
    r_squared = sm$r.squared,
    rank = rank,
    n_regressors = ncol(design),
    collinear = collinear,
    basis = basis,
    design = design
  ), class = "eda_glm")
}

#' Resolve a channel index from an index or label
#' @keywords internal
.eda_resolve_channel <- function(x, channel) {
  if (is.character(channel)) {
    labs <- as.character(SummarizedExperiment::colData(x)$label)
    idx <- match(channel, labs)
    if (is.na(idx)) stop("channel label not found: ", channel, call. = FALSE)
    idx
  } else {
    idx <- as.integer(channel)
    if (idx < 1 || idx > ncol(SummarizedExperiment::assay(x)))
      stop("channel index out of range", call. = FALSE)
    idx
  }
}
