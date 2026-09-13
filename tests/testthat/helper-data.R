#' Create test EDA PhysioExperiment
#' @param n_time Number of time points
#' @param n_channels Number of EDA channels
#' @param sr Sampling rate in Hz
#' @return PhysioExperiment with simulated EDA data
make_eda <- function(n_time = 6000, n_channels = 1, sr = 10) {
  edaSimulate(n_time = n_time, n_channels = n_channels, sr = sr,
              scr_count = 5, scl_level = 5.0, scr_amplitude = 0.5,
              noise_sd = 0.01, seed = 42)
}

#' Create test EDA PhysioExperiment with stimulus events
#' @param n_time Number of time points
#' @param n_channels Number of EDA channels
#' @param sr Sampling rate in Hz
#' @param n_events Number of stimulus events
#' @param event_interval Interval between events in seconds
#' @return PhysioExperiment with simulated EDA data and events
make_eda_event <- function(n_time = 6000, n_channels = 1, sr = 10,
                           n_events = 4, event_interval = 10) {
  x <- edaSimulate(n_time = n_time, n_channels = n_channels, sr = sr,
                   scr_count = n_events, scl_level = 5.0,
                   scr_amplitude = 0.8, noise_sd = 0.01, seed = 42)

  # Add stimulus events at regular intervals
  duration_sec <- (n_time - 1) / sr
  margin <- duration_sec * 0.1
  event_onsets <- seq(margin, duration_sec - margin,
                      length.out = n_events)

  x <- setEvents(x, PhysioEvents(
    onset = event_onsets,
    duration = rep(0.5, n_events),
    type = rep("stimulus", n_events),
    value = paste0("S", seq_len(n_events))
  ))

  x
}

#' Create a clean EDA fixture with SCRs whose decay matches the CDA/DDA tau2
#'
#' Builds a single-channel EDA signal from a smooth tonic plus SCRs at known
#' onsets, each a unit-peak Bateman response (tau1 = 0.75, tau2 = 2.0) so there
#' is no tau mismatch with the default CDA/DDA decomposition. Used to test SCR
#' recovery and CDA/DDA agreement on clean data.
#' @param onsets SCR onset times in seconds.
#' @param amps SCR peak amplitudes in microsiemens.
#' @param sr Sampling rate in Hz.
#' @param dur Duration in seconds.
#' @param noise_sd Gaussian noise SD.
#' @param seed RNG seed.
#' @return A list with the PhysioExperiment \code{pe} and the true \code{onsets}.
make_eda_scr_clean <- function(onsets = c(10, 28, 46, 64),
                               amps = c(0.8, 0.6, 0.9, 0.7),
                               sr = 10, dur = 80, noise_sd = 0.002, seed = 1) {
  set.seed(seed)
  n <- as.integer(sr * dur)
  t <- (seq_len(n) - 1) / sr
  phasic <- numeric(n)
  for (i in seq_along(onsets)) {
    phasic <- phasic + amps[i] *
      scrf(t - onsets[i], tau1 = 0.75, tau2 = 2.0, form = "bateman",
           normalize = "peak")
  }
  sig <- 5 + 0.2 * sin(2 * pi * 0.01 * t) + phasic + rnorm(n, 0, noise_sd)
  sig <- pmax(sig, 0.001)
  pe <- PhysioExperiment(
    assays = list(raw = matrix(sig, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "EDA1", type = "EDA", unit = "uS"),
    samplingRate = sr
  )
  list(pe = pe, onsets = onsets, amps = amps)
}

#' Create an EDA artifact-detection fixture with known artifact epochs
#'
#' Builds a clean SCL + SCR signal and injects spike / step / motion-burst
#' artifacts into a set of 5 s epochs. Optionally adds a 3-axis accelerometer
#' whose motion energy rises exactly in the injected epochs. Returns the
#' PhysioExperiment plus the ground-truth epoch labels.
#' @param seed RNG seed.
#' @param sr Sampling rate in Hz.
#' @param dur Duration in seconds.
#' @param epoch_sec Epoch length in seconds.
#' @param n_art Number of artifact epochs (0 = clean).
#' @param with_acc Add a 3-axis accelerometer channel.
#' @return A list with \code{pe}, \code{label} (per-epoch logical), \code{n_ep}.
make_eda_artifact_fixture <- function(seed = 1, sr = 8, dur = 240,
                                      epoch_sec = 5, n_art = 8, with_acc = FALSE) {
  set.seed(seed)
  n <- sr * dur
  t <- (0:(n - 1)) / sr
  sig <- 5 + 0.3 * sin(2 * pi * 0.008 * t) + rnorm(n, 0, 0.01)
  for (o in seq(15, dur - 15, length.out = 10)) {
    tt <- t - o; tt[tt < 0] <- NA
    sh <- 0.5 * (exp(-tt / 3) - exp(-tt / 0.75)); sh[is.na(sh)] <- 0; sh[sh < 0] <- 0
    sig <- sig + sh
  }
  L <- as.integer(round(epoch_sec * sr))
  n_ep <- floor(n / L)
  label <- rep(FALSE, n_ep)
  acc <- matrix(rnorm(n * 3, 0, 0.01), ncol = 3)
  acc[, 3] <- acc[, 3] + 1                          # gravity on z-axis
  if (n_art > 0) {
    ae <- sort(sample(2:(n_ep - 1), n_art))
    label[ae] <- TRUE
    for (e in ae) {
      rng <- ((e - 1) * L + 1):(e * L)
      kind <- sample(c("spike", "step", "motion"), 1)
      if (kind == "spike") {
        p <- rng[sample(length(rng), 1)]
        sig[p] <- sig[p] + sample(c(-1, 1), 1) * runif(1, 2, 5)
      } else if (kind == "step") {
        sig[rng] <- sig[rng] + sample(c(-1, 1), 1) * runif(1, 1, 3)
      } else {
        sig[rng] <- sig[rng] + rnorm(length(rng), 0, 0.6) + 0.8 * sin(2 * pi * 2 * t[rng])
      }
      acc[rng, ] <- acc[rng, ] + rnorm(length(rng) * 3, 0, 0.5)
    }
  }
  if (with_acc) {
    assays <- list(raw = cbind(sig, acc))
    cd <- S4Vectors::DataFrame(
      label = c("EDA1", "ACCx", "ACCy", "ACCz"),
      type = c("EDA", "ACC", "ACC", "ACC"),
      unit = c("uS", "g", "g", "g"))
  } else {
    assays <- list(raw = matrix(sig, ncol = 1))
    cd <- S4Vectors::DataFrame(label = "EDA1", type = "EDA", unit = "uS")
  }
  pe <- PhysioExperiment(assays = assays, colData = cd, samplingRate = sr)
  list(pe = pe, label = label, n_ep = n_ep)
}

#' GLM fixture: SCRs with known per-condition amplitudes
#'
#' Builds an EDA signal from a smooth tonic plus condition-grouped SCRs
#' (unit-peak Bateman responses scaled by the per-condition amplitude), for
#' testing \code{\link{edaGLM}} amplitude recovery.
#' @return A list with \code{pe}, \code{events} (onset/condition), \code{conds}
#'   (condition/amp table).
make_eda_glm_fixture <- function(seed = 1, sr = 10, dur = 200, noise = 0.02,
                                 amps = c(A = 1.0, B = 0.5, C = 0.25)) {
  set.seed(seed)
  n <- sr * dur
  t <- (0:(n - 1)) / sr
  conds <- data.frame(condition = names(amps), amp = as.numeric(amps),
                      stringsAsFactors = FALSE)
  events <- do.call(rbind, lapply(seq_len(nrow(conds)), function(i) {
    on <- seq(10 + (i - 1) * 4, dur - 10, by = 18)[1:8]
    data.frame(onset = on, condition = conds$condition[i], stringsAsFactors = FALSE)
  }))
  amp_of <- conds$amp[match(events$condition, conds$condition)]
  sig <- 5 + 0.3 * sin(2 * pi * 0.005 * t)
  for (k in seq_len(nrow(events))) {
    sig <- sig + amp_of[k] *
      scrf(t - events$onset[k], tau1 = 0.75, tau2 = 2.0, normalize = "peak")
  }
  sig <- sig + rnorm(n, 0, noise)
  pe <- PhysioExperiment(
    assays = list(raw = matrix(pmax(sig, 0.001), ncol = 1)),
    colData = S4Vectors::DataFrame(label = "EDA1", type = "EDA", unit = "uS"),
    samplingRate = sr)
  list(pe = pe, events = events, conds = conds)
}

#' DCM fixture: SCRs with known per-event graded amplitudes
#'
#' Builds an EDA signal with well-separated SCRs of continuously varying
#' amplitude, for testing \code{\link{edaDCM}} per-event amplitude recovery.
#' @return A list with \code{pe}, \code{events} (onset), \code{amps} (true
#'   per-event amplitudes).
make_eda_dcm_fixture <- function(seed = 1, sr = 10, dur = 300, noise = 0.02,
                                 gap = 14) {
  set.seed(seed)
  n <- sr * dur
  t <- (0:(n - 1)) / sr
  onsets <- seq(12, dur - 15, by = gap)
  ne <- length(onsets)
  amps <- runif(ne, 0.3, 1.5)
  sig <- 5 + 0.3 * sin(2 * pi * 0.005 * t)
  for (k in seq_len(ne)) {
    sig <- sig + amps[k] *
      scrf(t - onsets[k], tau1 = 0.75, tau2 = 2.0, normalize = "peak")
  }
  sig <- sig + rnorm(n, 0, noise)
  pe <- PhysioExperiment(
    assays = list(raw = matrix(pmax(sig, 0.001), ncol = 1)),
    colData = S4Vectors::DataFrame(label = "EDA1", type = "EDA", unit = "uS"),
    samplingRate = sr)
  list(pe = pe, events = data.frame(onset = onsets), amps = amps)
}
