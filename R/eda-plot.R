#' Plot EDA Time Series
#'
#' Plots electrodermal activity signals as time series using base R graphics.
#' Multiple channels are displayed in stacked panels. Optionally overlays
#' vertical event markers.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param channels Integer or character vector specifying which channels to
#'   plot. If NULL (default), all channels are plotted.
#' @param time_range Numeric vector of length 2 giving start and end times
#'   in seconds for zooming (e.g., \code{c(10, 30)}). If NULL, the full
#'   signal is plotted.
#' @param assay_name Name of the assay to plot. If NULL, uses
#'   \code{defaultAssay(x)}.
#' @param show_events Logical; if TRUE (default) and events exist, draws
#'   vertical dashed lines at event onsets.
#' @param main Character string for the plot title (default: "EDA Signal").
#' @param col Vector of colors for the channels. If NULL, uses default palette.
#' @param ... Additional arguments passed to \code{plot()}.
#' @return Invisible NULL. Called for its side effect of producing a plot.
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{plotDecompose}} for decomposition visualization,
#'   \code{\link{plotPeaks}} for SCR peak visualization,
#'   \code{\link{edaDecompose}} for tonic/phasic decomposition
#' @export
plotEda <- function(x, channels = NULL, time_range = NULL,
                    assay_name = NULL, show_events = TRUE,
                    main = "EDA Signal", col = NULL, ...) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  # Channel labels
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% colnames(cd)) {
    ch_labels <- as.character(cd$label)
  } else {
    ch_labels <- paste0("Ch", seq_len(n_channels))
  }

  # Resolve channel selection
  if (!is.null(channels)) {
    if (is.character(channels)) {
      ch_idx <- match(channels, ch_labels)
      ch_idx <- ch_idx[!is.na(ch_idx)]
    } else {
      ch_idx <- as.integer(channels)
      ch_idx <- ch_idx[ch_idx >= 1L & ch_idx <= n_channels]
    }
    if (length(ch_idx) == 0) {
      stop("No valid channels selected", call. = FALSE)
    }
  } else {
    ch_idx <- seq_len(n_channels)
  }

  # Time vector
  time_sec <- (seq_len(n_time) - 1) / sr

  # Apply time range
  if (!is.null(time_range)) {
    stopifnot(is.numeric(time_range), length(time_range) == 2)
    t_idx <- which(time_sec >= time_range[1] & time_sec <= time_range[2])
    if (length(t_idx) == 0) {
      stop("time_range results in no data points", call. = FALSE)
    }
  } else {
    t_idx <- seq_len(n_time)
  }

  time_sub <- time_sec[t_idx]

  # Colors
  n_plot <- length(ch_idx)
  if (is.null(col)) {
    if (n_plot <= 8) {
      col <- seq_len(n_plot)
    } else {
      col <- rep_len(seq_len(8), n_plot)
    }
  } else {
    col <- rep_len(col, n_plot)
  }

  # Get events if needed
  event_onsets <- NULL
  if (show_events) {
    events <- getEvents(x)
    ev_df <- events@events
    if (nrow(ev_df) > 0) {
      event_onsets <- as.numeric(ev_df$onset)
      # Filter to time range
      if (!is.null(time_range)) {
        event_onsets <- event_onsets[event_onsets >= time_range[1] &
                                      event_onsets <= time_range[2]]
      }
      if (length(event_onsets) == 0) event_onsets <- NULL
    }
  }

  # Plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (n_plot > 1) {
    par(mfrow = c(n_plot, 1), mar = c(4, 4, 2, 1))
  }

  for (i in seq_along(ch_idx)) {
    ch <- ch_idx[i]
    sig <- data[t_idx, ch]
    plot_title <- if (n_plot == 1) main else paste0(main, " - ", ch_labels[ch])

    plot(time_sub, sig, type = "l", col = col[i],
         xlab = "Time (s)", ylab = expression(mu * S),
         main = plot_title, ...)

    if (!is.null(event_onsets)) {
      abline(v = event_onsets, lty = 2, col = "gray50")
    }
  }

  invisible(NULL)
}


#' Plot EDA Decomposition
#'
#' Displays the tonic/phasic decomposition of an EDA signal in a multi-panel
#' layout. Requires that \code{edaDecompose()} has been run first so that
#' "tonic" and "phasic" assays exist.
#'
#' @param x A PhysioExperiment object with tonic and phasic assays.
#' @param channel Integer or character specifying which channel to plot
#'   (default: 1).
#' @param time_range Numeric vector of length 2 giving start and end times
#'   in seconds. If NULL, the full signal is plotted.
#' @param main Character string for the overall plot title
#'   (default: "EDA Decomposition").
#' @param ... Additional arguments passed to \code{plot()}.
#' @return Invisible NULL. Called for its side effect of producing a plot.
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{plotEda}} for basic EDA time series plotting,
#'   \code{\link{edaDecompose}} for performing the decomposition,
#'   \code{\link{plotPeaks}} for SCR peak visualization
#' @export
plotDecompose <- function(x, channel = 1, time_range = NULL,
                          main = "EDA Decomposition", ...) {
  stopifnot(inherits(x, "PhysioExperiment"))

  anames <- SummarizedExperiment::assayNames(x)
  if (!("tonic" %in% anames) || !("phasic" %in% anames)) {
    stop("Assays 'tonic' and 'phasic' are required. Run edaDecompose() first.",
         call. = FALSE)
  }

  sr <- samplingRate(x)
  n_time <- nrow(SummarizedExperiment::assay(x, "tonic"))

  # Resolve channel
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% colnames(cd)) {
    ch_labels <- as.character(cd$label)
  } else {
    ch_labels <- paste0("Ch", seq_len(ncol(SummarizedExperiment::assay(x, "tonic"))))
  }

  if (is.character(channel)) {
    ch_idx <- match(channel, ch_labels)
    if (is.na(ch_idx)) stop("Channel not found: ", channel, call. = FALSE)
  } else {
    ch_idx <- as.integer(channel)
  }

  # Get data
  raw_assay <- defaultAssay(x)
  original <- SummarizedExperiment::assay(x, raw_assay)[, ch_idx]
  tonic <- SummarizedExperiment::assay(x, "tonic")[, ch_idx]
  phasic <- SummarizedExperiment::assay(x, "phasic")[, ch_idx]
  has_driver <- "driver" %in% anames
  if (has_driver) {
    driver <- SummarizedExperiment::assay(x, "driver")[, ch_idx]
  }

  # Time vector
  time_sec <- (seq_len(n_time) - 1) / sr

  # Apply time range
  if (!is.null(time_range)) {
    stopifnot(is.numeric(time_range), length(time_range) == 2)
    t_idx <- which(time_sec >= time_range[1] & time_sec <= time_range[2])
    if (length(t_idx) == 0) {
      stop("time_range results in no data points", call. = FALSE)
    }
  } else {
    t_idx <- seq_len(n_time)
  }

  time_sub <- time_sec[t_idx]
  original <- original[t_idx]
  tonic <- tonic[t_idx]
  phasic <- phasic[t_idx]
  if (has_driver) driver <- driver[t_idx]

  # Layout
  n_panels <- if (has_driver) 4L else 3L
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(n_panels, 1), mar = c(4, 4, 2, 1))

  ch_label <- ch_labels[ch_idx]

  # Panel 1: Original + tonic overlay
  plot(time_sub, original, type = "l", col = "black",
       xlab = "Time (s)", ylab = expression(mu * S),
       main = paste0(main, " - ", ch_label, " (Original + Tonic)"), ...)
  lines(time_sub, tonic, col = "blue", lwd = 2)
  legend("topright", legend = c("Original", "Tonic"),
         col = c("black", "blue"), lty = 1, lwd = c(1, 2), cex = 0.8)

  # Panel 2: Tonic component
  plot(time_sub, tonic, type = "l", col = "blue",
       xlab = "Time (s)", ylab = expression(mu * S),
       main = "Tonic (SCL)", ...)

  # Panel 3: Phasic component
  plot(time_sub, phasic, type = "l", col = "red",
       xlab = "Time (s)", ylab = expression(mu * S),
       main = "Phasic (SCR)", ...)

  # Panel 4: Driver (if exists)
  if (has_driver) {
    plot(time_sub, driver, type = "l", col = "darkgreen",
         xlab = "Time (s)", ylab = "Driver",
         main = "Sudomotor Nerve Activity (Driver)", ...)
  }

  invisible(NULL)
}


#' Plot EDA Signal with SCR Peaks
#'
#' Plots an EDA signal with detected SCR peaks, onsets, and amplitude lines
#' annotated. If no peaks data.frame is supplied, peaks are computed
#' automatically via \code{\link{edaPeaks}}.
#'
#' @param x A PhysioExperiment object containing EDA data.
#' @param peaks A data.frame as returned by \code{\link{edaPeaks}}, or NULL
#'   to compute peaks automatically.
#' @param channel Integer or character specifying which channel to plot
#'   (default: 1).
#' @param time_range Numeric vector of length 2 giving start and end times
#'   in seconds. If NULL, the full signal is plotted.
#' @param assay_name Name of the assay to plot. If NULL, uses "phasic" if
#'   available, otherwise \code{defaultAssay(x)}.
#' @param main Character string for the plot title (default: "SCR Peaks").
#' @param ... Additional arguments passed to \code{plot()}.
#' @return Invisible NULL. Called for its side effect of producing a plot.
#' @references
#' Boucsein, W. (2012). \emph{Electrodermal Activity}. 2nd ed. Springer.
#' \doi{10.1007/978-1-4614-1126-0}
#' @seealso \code{\link{plotEda}} for basic EDA time series plotting,
#'   \code{\link{plotDecompose}} for decomposition visualization,
#'   \code{\link{edaPeaks}} for SCR peak detection
#' @export
plotPeaks <- function(x, peaks = NULL, channel = 1,
                      time_range = NULL, assay_name = NULL,
                      main = "SCR Peaks", ...) {
  stopifnot(inherits(x, "PhysioExperiment"))

  # Resolve assay
  if (is.null(assay_name)) {
    anames <- SummarizedExperiment::assayNames(x)
    if ("phasic" %in% anames) {
      assay_name <- "phasic"
    } else {
      assay_name <- defaultAssay(x)
    }
  }

  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)

  # Channel labels
  cd <- SummarizedExperiment::colData(x)
  if ("label" %in% colnames(cd)) {
    ch_labels <- as.character(cd$label)
  } else {
    ch_labels <- paste0("Ch", seq_len(ncol(data)))
  }

  # Resolve channel
  if (is.character(channel)) {
    ch_idx <- match(channel, ch_labels)
    if (is.na(ch_idx)) stop("Channel not found: ", channel, call. = FALSE)
  } else {
    ch_idx <- as.integer(channel)
  }

  ch_label <- ch_labels[ch_idx]

  # Compute peaks if not provided
  if (is.null(peaks)) {
    peaks <- edaPeaks(x, assay_name = assay_name)
  }

  # Filter peaks to current channel
  ch_peaks <- peaks[peaks$channel == ch_label, , drop = FALSE]

  # Signal for the selected channel
  sig <- data[, ch_idx]

  # Time vector
  time_sec <- (seq_len(n_time) - 1) / sr

  # Apply time range
  if (!is.null(time_range)) {
    stopifnot(is.numeric(time_range), length(time_range) == 2)
    t_idx <- which(time_sec >= time_range[1] & time_sec <= time_range[2])
    if (length(t_idx) == 0) {
      stop("time_range results in no data points", call. = FALSE)
    }
    # Filter peaks to time range
    if (nrow(ch_peaks) > 0) {
      ch_peaks <- ch_peaks[ch_peaks$peak_sec >= time_range[1] &
                             ch_peaks$peak_sec <= time_range[2], , drop = FALSE]
    }
  } else {
    t_idx <- seq_len(n_time)
  }

  time_sub <- time_sec[t_idx]
  sig_sub <- sig[t_idx]

  # Plot signal
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  plot(time_sub, sig_sub, type = "l", col = "black",
       xlab = "Time (s)", ylab = expression(mu * S),
       main = paste0(main, " - ", ch_label), ...)

  # Annotate peaks
  n_peaks <- nrow(ch_peaks)
  if (n_peaks > 0) {
    # Onset markers (triangles)
    onset_vals <- sig[ch_peaks$onset_sample]
    points(ch_peaks$onset_sec, onset_vals, pch = 17, col = "blue", cex = 1.2)

    # Peak markers (filled circles)
    peak_vals <- sig[ch_peaks$peak_sample]
    points(ch_peaks$peak_sec, peak_vals, pch = 19, col = "red", cex = 1.2)

    # Amplitude lines
    segments(ch_peaks$peak_sec, onset_vals,
             ch_peaks$peak_sec, peak_vals,
             col = "green", lty = 2, lwd = 1.5)

    legend("topright",
           legend = c(paste0("Onsets (", n_peaks, ")"),
                      paste0("Peaks (", n_peaks, ")"),
                      "Amplitude"),
           pch = c(17, 19, NA),
           lty = c(NA, NA, 2),
           col = c("blue", "red", "green"),
           cex = 0.8)
  } else {
    legend("topright", legend = "No peaks detected", cex = 0.8)
  }

  invisible(NULL)
}
