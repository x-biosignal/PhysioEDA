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

  # --- Step 1: Build Bateman impulse response (single kernel source: scrf) ---
  h_len <- min(n, as.integer(15 * sr))
  t_h <- seq(0, by = 1 / sr, length.out = h_len)
  h <- scrf(t_h, tau1 = tau1, tau2 = tau2, form = "bateman", normalize = "none")
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


#' Bateman ARMA Matrices for cvxEDA
#'
#' Builds the banded auto-regressive (\code{A}) and moving-average (\code{M})
#' matrices of the discretised biexponential (Bateman) system used by the exact
#' cvxEDA formulation (Greco et al., 2016). The latent variable \code{q} relates
#' to the SMNA driver by \code{p = A \%*\% q} and to the phasic component by
#' \code{r = M \%*\% q}.
#'
#' @param n Integer number of samples.
#' @param dt Sampling interval in seconds (\code{1 / sr}).
#' @param tau_rise SCR rise time constant in seconds.
#' @param tau_decay SCR decay time constant in seconds.
#' @return A list with the \code{n x n} matrices \code{A} and \code{M}.
#' @keywords internal
.bateman_arma <- function(n, dt, tau_rise, tau_decay) {
  a1 <- 1 / min(tau_rise, tau_decay)   # fast rate
  a0 <- 1 / max(tau_rise, tau_decay)   # slow rate
  if (abs(a1 - a0) < 1e-8) a1 <- a1 + 1e-6
  ar <- c((a1 * dt + 2) * (a0 * dt + 2),
          2 * a1 * a0 * dt^2 - 8,
          (a1 * dt - 2) * (a0 * dt - 2)) / ((a1 - a0) * dt^2)
  ma <- c(1, 2, 1)
  A <- matrix(0, n, n)
  M <- matrix(0, n, n)
  for (i in 3:n) {
    A[i, i] <- ar[1]; A[i, i - 1] <- ar[2]; A[i, i - 2] <- ar[3]
    M[i, i] <- ma[1]; M[i, i - 1] <- ma[2]; M[i, i - 2] <- ma[3]
  }
  list(A = A, M = M)
}


#' Triangular (order-1) B-spline Tonic Basis for cvxEDA
#'
#' Builds the triangular spline regressor matrix \code{B} whose columns are
#' overlapping tent functions spaced every \code{delta_knot} seconds, following
#' Greco et al. (2016). The tonic (SCL) component is modelled as
#' \code{B \%*\% l + C \%*\% d} where \code{C} is a linear trend.
#'
#' @param n Integer number of samples.
#' @param dt Sampling interval in seconds.
#' @param delta_knot Knot spacing in seconds (default: 10).
#' @return An \code{n x nB} matrix of spline regressors.
#' @keywords internal
.spline_tonic_basis <- function(n, dt, delta_knot = 10) {
  dk <- max(2L, as.integer(round(delta_knot / dt)))
  spl <- c(seq_len(dk - 1), rev(seq_len(dk))) / dk   # triangular, length 2dk-1
  L <- length(spl)
  half <- L %/% 2L
  knot_starts <- seq(0L, n - 1L, by = dk)
  nB <- length(knot_starts)
  B <- matrix(0, n, nB)
  offsets <- seq(-half, L - 1L - half)
  for (j in seq_len(nB)) {
    idx <- knot_starts[j] + offsets + 1L
    ok <- idx >= 1 & idx <= n
    B[idx[ok], j] <- spl[ok]
  }
  B
}


#' Exact cvxEDA Quadratic-Program Core (Single Channel, Fixed Rate)
#'
#' Solves the Greco et al. (2016) cvxEDA convex problem for one standardised
#' signal at a single sampling rate. The signal is decomposed as
#' \code{y = M q + B l + C d + noise} with a sparse non-negative driver
#' \code{p = A q >= 0}, minimising
#' \code{0.5 * sum((M q + B l + C d - y)^2) + alpha * sum(A q) + 0.5 * gamma * sum(l^2)}
#' subject to \code{A q >= 0}. This is a quadratic program solved with an
#' active-set solver (\code{quadprog::solve.QP}); the returned optimum is exact.
#'
#' @param y Numeric EDA signal (single channel).
#' @param sr Sampling rate in Hz.
#' @param tau1 SCR rise time constant in seconds (default: 0.75).
#' @param tau2 SCR decay time constant in seconds (default: 2.0).
#' @param alpha L1 sparsity penalty on the driver (default: 0.01).
#' @param gamma L2 smoothness weight on the spline tonic (default: 0.1).
#' @param delta_knot Tonic spline knot spacing in seconds (default: 10).
#' @return A list with \code{tonic}, \code{phasic}, \code{driver} (each of length
#'   \code{length(y)}), the final QP \code{objective}, the objective at the
#'   trivial (driver-free) initialisation \code{obj_init}, and \code{converged}.
#' @references
#' Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
#' electrodermal activity processing." \emph{IEEE Transactions on Biomedical
#' Engineering}, 63(4), 797-804. \doi{10.1109/TBME.2015.2474131}
#' @keywords internal
.cvxeda_qp_core <- function(y, sr, tau1 = 0.75, tau2 = 2.0,
                            alpha = 0.01, gamma = 0.1, delta_knot = 10) {
  n <- length(y)
  dt <- 1 / sr
  # Standardise so alpha/gamma act on a scale-invariant signal.
  mu <- mean(y)
  sdv <- stats::sd(y)
  if (sdv < 1e-12) sdv <- 1
  ys <- (y - mu) / sdv

  ba <- .bateman_arma(n, dt, tau1, tau2)
  A <- ba$A
  M <- ba$M
  B <- .spline_tonic_basis(n, dt, delta_knot)
  C <- cbind(1, seq_len(n) / n)
  nb <- ncol(B)
  nc <- ncol(C)
  p <- n + nb + nc

  W <- cbind(M, B, C)
  D <- crossprod(W)
  gl <- (n + 1):(n + nb)
  diag(D)[gl] <- diag(D)[gl] + gamma
  # Ridge to guarantee a positive-definite Dmat for the Cholesky in solve.QP.
  D <- D + diag(1e-8 * mean(diag(D)) + 1e-10, p)
  dvec <- as.numeric(crossprod(W, ys))
  dvec[1:n] <- dvec[1:n] - alpha * as.numeric(crossprod(A, rep(1, n)))
  # Constraint A q >= 0 (one row of A per constraint).
  Amat <- rbind(t(A), matrix(0, nb + nc, n))
  bvec <- rep(0, n)

  sol <- quadprog::solve.QP(D, dvec, Amat, bvec, meq = 0)
  x <- sol$solution
  q <- x[1:n]
  l <- x[gl]
  d <- x[(n + nb + 1):p]

  phasic_s <- as.numeric(M %*% q)
  tonic_s <- as.numeric(B %*% l + C %*% d)
  driver_s <- pmax(0, as.numeric(A %*% q))

  obj_fn <- function(qq, ll, dd) {
    rec <- as.numeric(M %*% qq + B %*% ll + C %*% dd)
    0.5 * sum((rec - ys)^2) + alpha * sum(pmax(0, as.numeric(A %*% qq))) +
      0.5 * gamma * sum(ll^2)
  }
  # Objective at the trivial driver-free start (tonic = LS fit of y on [B, C]).
  BC <- cbind(B, C)
  init <- tryCatch(as.numeric(qr.coef(qr(BC), ys)),
                   error = function(e) rep(0, nb + nc))
  init[is.na(init)] <- 0
  obj_init <- obj_fn(rep(0, n), init[1:nb], init[(nb + 1):(nb + nc)])
  obj_fin <- obj_fn(q, l, d)

  list(
    tonic = tonic_s * sdv + mu,
    phasic = phasic_s * sdv,
    driver = driver_s,
    objective = obj_fin,
    obj_init = obj_init,
    converged = is.finite(obj_fin) && obj_fin <= obj_init + 1e-6
  )
}


#' Exact cvxEDA Decomposition for a Single Channel (Greco QP)
#'
#' Full-signal wrapper around \code{\link{.cvxeda_qp_core}}: it decimates the
#' signal to a tractable rate for the dense quadratic program (electrodermal
#' activity is slow, so a few Hz suffices), solves the exact cvxEDA QP, and
#' interpolates the tonic / phasic / driver components back to the native grid.
#' If \pkg{quadprog} is unavailable, or the solve fails, it degrades gracefully
#' to the Wiener approximation (\code{\link{.cvxeda_decompose}}).
#'
#' @param signal Numeric vector of EDA signal values.
#' @param sr Sampling rate in Hz.
#' @param tau1 SCR rise time constant in seconds (default: 0.75).
#' @param tau2 SCR decay time constant in seconds (default: 2.0).
#' @param alpha L1 sparsity penalty on the driver (default: 0.01).
#' @param gamma L2 smoothness weight on the tonic (default: 0.1).
#' @param delta_knot Tonic spline knot spacing in seconds (default: 10).
#' @param qp_sr Target rate (Hz) for the QP solve; the signal is decimated
#'   towards this rate (default: 4).
#' @param max_n Maximum number of samples in the QP; further decimation is
#'   applied above this (default: 2000).
#' @return A list with \code{tonic}, \code{phasic}, and \code{driver} on the
#'   native grid, plus \code{driver_qp} / \code{t_qp} / \code{sr_qp} at the QP
#'   rate, the QP \code{objective} and \code{obj_init}, \code{converged}, and the
#'   \code{engine} ("qp" or "wiener_fallback") and \code{solver} used.
#' @references
#' Greco, A., et al. (2016). "cvxEDA: A convex optimization approach to
#' electrodermal activity processing." \emph{IEEE Transactions on Biomedical
#' Engineering}, 63(4), 797-804. \doi{10.1109/TBME.2015.2474131}
#' @seealso \code{\link{.cvxeda_decompose}} for the fast Wiener approximation,
#'   \code{\link{edaDecompose}} for the user-facing wrapper.
#' @keywords internal
.cvxeda_qp <- function(signal, sr, tau1 = 0.75, tau2 = 2.0,
                       alpha = 0.01, gamma = 0.1, delta_knot = 10,
                       qp_sr = 4, max_n = 2000) {
  n0 <- length(signal)
  stopifnot(n0 > 0, sr > 0, tau1 > 0, tau2 > 0, alpha >= 0, gamma >= 0)

  fallback <- function(msg) {
    warning(msg, call. = FALSE)
    res <- .cvxeda_decompose(signal, sr, tau1 = tau1, tau2 = tau2,
                             alpha = alpha, gamma = gamma)
    res$engine <- "wiener_fallback"
    res$solver <- "none"
    res$objective <- NA_real_
    res
  }

  if (!requireNamespace("quadprog", quietly = TRUE)) {
    return(fallback(paste0("cvxEDA QP requires the 'quadprog' package; ",
                           "falling back to the Wiener approximation ",
                           "(equivalent to method = 'cvxeda_fast').")))
  }

  # Decimate to keep the dense QP tractable.
  dec <- max(1L, as.integer(round(sr / qp_sr)))
  if (n0 %/% dec > max_n) dec <- as.integer(ceiling(n0 / max_n))
  if (dec > 1L) {
    grp <- rep(seq_len(ceiling(n0 / dec)), each = dec)[seq_len(n0)]
    y_qp <- as.numeric(tapply(signal, grp, mean))
    t_qp <- as.numeric(tapply((seq_len(n0) - 1L) / sr, grp, mean))
    sr_qp <- sr / dec
  } else {
    y_qp <- signal
    t_qp <- (seq_len(n0) - 1L) / sr
    sr_qp <- sr
  }

  core <- tryCatch(
    .cvxeda_qp_core(y_qp, sr_qp, tau1 = tau1, tau2 = tau2,
                    alpha = alpha, gamma = gamma, delta_knot = delta_knot),
    error = function(e) NULL)
  if (is.null(core)) {
    return(fallback(paste0("cvxEDA QP solve failed; falling back to the ",
                           "Wiener approximation.")))
  }

  t_full <- (seq_len(n0) - 1L) / sr
  up <- function(v) stats::approx(t_qp, v, xout = t_full, rule = 2)$y
  list(
    tonic = up(core$tonic),
    phasic = up(core$phasic),
    driver = pmax(0, up(core$driver)),
    driver_qp = core$driver,
    t_qp = t_qp,
    sr_qp = sr_qp,
    alpha = alpha,
    gamma = gamma,
    objective = core$objective,
    obj_init = core$obj_init,
    converged = core$converged,
    engine = "qp",
    solver = "quadprog"
  )
}
