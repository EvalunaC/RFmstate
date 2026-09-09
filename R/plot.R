#' Plot Aalen-Johansen Estimates
#'
#' Visualizes calendar-time state occupation, product-integral point estimates,
#' Nelson--Aalen cumulative cause-specific hazards, or hazard increments from
#' the Aalen--Johansen estimator.
#'
#' @param x An \code{aj_estimate} object.
#' @param type Character, one of \code{"state_occupation"} (default),
#'   \code{"stacked_transition_prob"}, \code{"cumulative_hazard"},
#'   \code{"hazard_increment"}. The legacy \code{"transition_intensity"}
#'   alias warns because increments are not intensities.
#' @param states Character vector of states to plot (default: all).
#'   For cumulative-hazard and hazard-increment plots, transitions are filtered
#'   by destination state.
#' @param ci One nonmissing logical value. This is a deprecated compatibility
#'   argument: \code{TRUE} warns because no validated AJ covariance or
#'   confidence band is returned.
#' @param col One or more colors for selected states/transitions. Values are
#'   recycled when necessary; \code{NULL} uses the package palette.
#' @param main Title (default: auto-generated).
#' @param xlab,ylab Axis labels.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#'
#' @return The input \code{x} object, returned invisibly. Called for its
#'   side effect of producing a plot.
#'
#' @section Limitations:
#' AJ plots are covariate-free calendar-time point estimates from the common
#' baseline. Confidence bands are unavailable. Hazard increments are discrete
#' Nelson--Aalen increments, not smoothed transition intensities. The legacy
#' \code{"transition_intensity"} type is deprecated and draws the same
#' hazard-increment plot with a warning. The historical
#' \code{"stacked_transition_prob"} label displays the selected state-occupation
#' components from the recorded initial state.
#'
#' @examples
#' ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
#' dat <- data.frame(id = 1:6, x = 1:6, time_B = 1:6,
#'                   censor = NA_real_)
#' long <- prepare_data(dat, "id", ms, list(B = "time_B"), "censor", "x")
#' aj <- aalen_johansen(long)
#' plot(aj, states = "B", xlab = "Study day", ylab = "Probability")
#' plot(aj, type = "hazard_increment", states = "B")
#'
#' @export
plot.aj_estimate <- function(x, type = c("state_occupation",
                                          "stacked_transition_prob",
                                          "cumulative_hazard",
                                          "hazard_increment",
                                          "transition_intensity"),
                             states = NULL, ci = FALSE, col = NULL,
                             main = NULL, xlab = "Time", ylab = NULL, ...) {
  type <- match.arg(type)

  states <- .validate_plot_states(states, x$structure$state_names)
  .validate_scalar_logical(ci, "ci")
  if (ci) {
    warning("AJ confidence intervals are unavailable: no validated covariance ",
            "is returned.", call. = FALSE)
  }
  ns <- length(states)
  if (is.null(col)) col <- .default_palette(ns)
  if (!length(col)) stop("'col' must contain at least one color.")
  col <- rep_len(col, ns)

  switch(type,
    "state_occupation" = {
      if (is.null(main)) main <- "State Occupation Probabilities (AJ)"
      .plot_state_occ(x$time, x$state_occ, states,
                      x$structure$state_names, col, main, xlab,
                      ylab %||% "State-occupation probability", ...)
    },
    "stacked_transition_prob" = {
      if (is.null(main)) main <- "Stacked Transition Probabilities (AJ)"
      .plot_trans_prob_aj(x, states, col, main, xlab,
                          ylab %||% "State-occupation probability", ...)
    },
    "cumulative_hazard" = {
      if (is.null(main)) main <- "Cumulative Hazards (Nelson-Aalen)"
      .plot_cum_hazard_aj(x, states, col, main, xlab,
                          ylab %||% "Cumulative cause-specific hazard", ...)
    },
    "hazard_increment" = {
      if (is.null(main)) main <- "Nelson-Aalen Hazard Increments"
      .plot_hazard_increment_aj(x, states, col, main, xlab,
                                ylab %||% "Nelson-Aalen hazard increment", ...)
    },
    "transition_intensity" = {
      warning("type = 'transition_intensity' is deprecated; use ",
              "type = 'hazard_increment'.", call. = FALSE)
      if (is.null(main)) main <- "Nelson-Aalen Hazard Increments"
      .plot_hazard_increment_aj(x, states, col, main, xlab,
                                ylab %||% "Nelson-Aalen hazard increment", ...)
    }
  )
  invisible(x)
}

#' Plot RF Multistate Predictions
#'
#' Visualizes predicted state occupation probabilities and transition
#' probabilities for individual patients.
#'
#' @param x An \code{rfmstate_pred} object.
#' @param type Character, one of \code{"state_occupation"} (default),
#'   \code{"transition_prob"}.
#' @param subject Integer, which subject to plot (default 1). Use 0 for
#'   mean across all subjects.
#' @param col One or more state colors. Values are recycled when necessary;
#'   \code{NULL} uses the package palette.
#' @param main Title.
#' @param states Occupied states to display.
#' @param xlab,ylab Axis labels.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#'
#' @return The input \code{x} object, returned invisibly. Called for its
#'   side effect of producing a plot.
#'
#' @section Limitations:
#' Curves are entry-conditioned on fresh entry into the fitted
#' \code{start_state}; they are not general Markov \eqn{P(s,t)} curves and have
#' no confidence bands. \code{subject = 0} is an arithmetic mean of profile
#' predictions, not a population-standardized estimator. The historical
#' \code{"transition_prob"} plot type displays the occupied-state components
#' conditional on the selected fresh-entry state.
#'
#' @examples
#' \donttest{
#' ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
#' dat <- data.frame(id = 1:60, x = seq(-1, 1, length.out = 60),
#'                   time_B = 1:60, censor = NA_real_)
#' long <- prepare_data(dat, "id", ms, list(B = "time_B"), "censor", "x")
#' fit <- rfmstate(long, num.trees = 20, min.node.size = 3,
#'                 min_events = 3, sparse_warning = Inf, seed = 42)
#' pred <- predict(fit, data.frame(x = 0), times = c(0, 5, 10),
#'                 target_grid_points = 64, max_grid_points = 1025)
#' plot(pred, states = c("A", "B"), xlab = "Elapsed day")
#' plot(pred, type = "transition_prob", states = "B")
#' }
#'
#' @export
plot.rfmstate_pred <- function(x, type = c("state_occupation",
                                            "transition_prob"),
                               subject = 1L, col = NULL,
                               main = NULL, states = NULL,
                               xlab = "Elapsed duration", ylab = NULL, ...) {
  type <- match.arg(type)
  state_names <- x$structure$state_names
  states <- .validate_plot_states(states, state_names)
  ns <- length(state_names)
  subject <- .validate_plot_subject(subject, x$n_subjects)

  if (is.null(col)) col <- .default_palette(ns)
  if (!length(col)) stop("'col' must contain at least one color.")
  col <- rep_len(col, ns)

  if (subject == 0) {
    # Mean across subjects
    occ <- apply(x$state_occ, c(2, 3), mean)
    subj_label <- "Mean"
  } else {
    if (subject > x$n_subjects) {
      stop("Subject ", subject, " exceeds number of subjects (",
           x$n_subjects, ").")
    }
    occ <- x$state_occ[subject, , ]
    subj_label <- paste("Subject", subject)
  }

  switch(type,
    "state_occupation" = {
      if (is.null(main)) {
        main <- paste("Predicted State Occupation -", subj_label)
      }
      .plot_state_occ_pred(x$time, occ, state_names, states, col, main,
                           xlab, ylab %||% "State-occupation probability", ...)
    },
    "transition_prob" = {
      if (is.null(main)) {
        main <- paste("Predicted Transition Probabilities -", subj_label)
      }
      if (subject == 0) {
        P_mean <- apply(x$P, c(2, 3, 4), mean)
      } else {
        P_mean <- x$P[subject, , , , drop = FALSE]
      }
      # Preserve the selected-start dimension even when it has length one.
      # Both a single profile and a profile mean are represented here as
      # start state x occupied state x time.
      P_mean <- array(
        as.numeric(P_mean),
        dim = dim(x$P)[-1L],
        dimnames = dimnames(x$P)[-1L]
      )
      .plot_trans_prob_pred(x$time, P_mean, state_names, states,
                            x$start_state, col, main, xlab,
                            ylab %||% "Entry-conditioned state probability", ...)
    }
  )
  invisible(x)
}

#' Plot Feature Importance
#'
#' Visualizes per-transition feature importance as a grouped barplot or
#' heatmap.
#'
#' @param x An \code{rfmstate_importance} object.
#' @param type Character, one of \code{"barplot"} (default),
#'   \code{"heatmap"}.
#' @param col Bar colors or heatmap gradient-anchor colors. \code{NULL} uses
#'   the package defaults.
#' @param main Title.
#' @param xlab,ylab Axis labels.
#' @param ... Additional graphical arguments. For a heatmap they are passed to
#'   \code{image}; for a barplot they are passed to \code{barplot}.
#'
#' @return The input \code{x} object, returned invisibly. Called for its
#'   side effect of producing a plot.
#'
#' @section Limitations:
#' Values are transition-specific ranger importance scores. Negative values may
#' arise from Monte Carlo noise, sparse events, correlated predictors, or
#' irrelevant variables and are not protective or causal effects. Scales may
#' differ across transitions.
#'
#' @examples
#' imp <- structure(list(
#'   importance_matrix = matrix(
#'     c(0.2, -0.1), nrow = 2,
#'     dimnames = list(c("x", "z"), "A->B")
#'   )
#' ), class = "rfmstate_importance")
#' plot(imp, xlab = "OOB loss increase", ylab = "Predictor")
#'
#' @export
plot.rfmstate_importance <- function(x,
                                     type = c("barplot", "heatmap"),
                                     col = NULL, main = NULL,
                                     xlab = NULL, ylab = NULL, ...) {
  type <- match.arg(type)

  switch(type,
    "barplot" = {
      if (is.null(main)) main <- "Feature Importance by Transition"
      .plot_importance_bar(
        x, col, main, xlab %||% "Permutation importance",
        ylab %||% "Variable", ...
      )
    },
    "heatmap" = {
      if (is.null(main)) main <- "Feature Importance Heatmap"
      .plot_importance_heat(
        x, col, main, xlab %||% "Transition", ylab %||% "Variable", ...
      )
    }
  )
  invisible(x)
}

#' Plot Diagnostics
#'
#' Visualizes genuine edge OOB concordance or patient-level cross-validated
#' full-state Brier scores.
#'
#' @param x An \code{rfmstate_diag} object.
#' @param type Character, one of \code{"brier"} (default),
#'   \code{"concordance"}.
#' @param col One or more plot colors; recycled when necessary.
#' @param main Title.
#' @param xlab,ylab Axis labels.
#' @param ... Additional arguments.
#'
#' @return The input \code{x} object, returned invisibly. Called for its
#'   side effect of producing a plot.
#'
#' @section Limitations:
#' Concordance plots contain separate edge-level ranger OOB statistics, not
#' assembled full-state validation. Brier plots require an object produced by
#' \code{diagnose(..., method = "cv")} and inherit its censoring/support
#' limitations. No bias--variance plot is available.
#'
#' @examples
#' d <- structure(list(
#'   concordance = data.frame(transition = "A->B", c_index = 0.7),
#'   brier = NULL
#' ), class = "rfmstate_diag")
#' plot(d, type = "concordance", xlab = "Edge", ylab = "OOB C-index")
#'
#' @export
plot.rfmstate_diag <- function(x, type = c("brier", "concordance"),
                               col = NULL, main = NULL,
                               xlab = NULL, ylab = NULL, ...) {
  type <- match.arg(type)

  switch(type,
    "brier" = {
      if (is.null(main)) main <- "Time-Dependent Brier Score"
      .plot_brier(x, col, main, xlab %||% "Time",
                  ylab %||% "Cross-validated Brier score", ...)
    },
    "concordance" = {
      if (is.null(main)) main <- "Concordance Index by Transition"
      .plot_concordance(
        x, col, main, xlab %||% "Transition", ylab %||% "OOB C-index", ...
      )
    }
  )
  invisible(x)
}

# ---- Internal plotting functions ----

#' Default color palette
#' @noRd
.default_palette <- function(n) {
  if (n <= 8) {
    cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
              "#66a61e", "#e6ab02", "#a6761d", "#666666")
    cols[seq_len(n)]
  } else {
    grDevices::rainbow(n, s = 0.7, v = 0.8)
  }
}

#' Validate a public plot state filter
#' @noRd
.validate_plot_states <- function(states, state_names) {
  if (is.null(states)) return(state_names)
  if (!is.character(states) || !length(states) || anyNA(states) ||
      any(!nzchar(states)) || anyDuplicated(states)) {
    stop("'states' must be a nonempty character vector of unique state names.")
  }
  invalid <- setdiff(states, state_names)
  if (length(invalid)) {
    stop("Unknown state(s): ", paste(invalid, collapse = ", "), ".")
  }
  states
}

#' Validate the public prediction-plot subject selector
#' @noRd
.validate_plot_subject <- function(subject, n_subjects) {
  if (!is.numeric(subject) || length(subject) != 1L || !is.finite(subject) ||
      subject < 0 || subject != as.integer(subject)) {
    stop("'subject' must be zero or one positive integer profile index.")
  }
  subject <- as.integer(subject)
  if (subject > n_subjects) {
    stop("Subject ", subject, " exceeds number of subjects (", n_subjects, ").")
  }
  subject
}

#' Filter transitions by selected destination states
#' @noRd
.filter_plot_transitions <- function(trans_list, states) {
  trans_list[trans_list$to %in% states, , drop = FALSE]
}

#' Plot state occupation from AJ
#' @noRd
.plot_state_occ <- function(times, state_occ, states, all_states,
                            col, main, xlab, ylab, ...) {
  state_idx <- match(states, all_states)

  plot(NULL, xlim = range(times), ylim = c(0, 1),
       xlab = xlab, ylab = ylab,
       main = main, ...)

  for (i in seq_along(state_idx)) {
    idx <- state_idx[i]
    lines(times, state_occ[, idx], col = col[i], lwd = 2)

  }

  legend("topright", legend = states, col = col[seq_along(states)],
         lwd = 2, bty = "n", cex = 0.8)
}

#' Plot transition probabilities from AJ
#' @noRd
.plot_trans_prob_aj <- function(aj, states, col, main, xlab, ylab, ...) {
  state_names <- aj$structure$state_names
  n_times <- length(aj$time)
  ns <- length(states)

  # Extract the recorded initial-state row for each destination state.
  initial_idx <- match(aj$initial_state, state_names)
  prob_mat <- matrix(0, nrow = n_times, ncol = ns)
  for (j in seq_len(ns)) {
    j_idx <- match(states[j], state_names)
    prob_mat[, j] <- vapply(aj$trans_prob, function(P) P[initial_idx, j_idx],
                            numeric(1))
  }

  # Cumulative sums for stacking (bottom to top)
  cum_lower <- matrix(0, nrow = n_times, ncol = ns)
  cum_upper <- matrix(0, nrow = n_times, ncol = ns)
  cum_upper[, 1] <- prob_mat[, 1]
  if (ns > 1) {
    for (j in 2:ns) {
      cum_lower[, j] <- cum_upper[, j - 1]
      cum_upper[, j] <- cum_lower[, j] + prob_mat[, j]
    }
  }

  plot(NULL, xlim = range(aj$time), ylim = c(0, 1),
       xlab = xlab, ylab = ylab,
       main = main, ...)

  # Draw stacked polygons from top to bottom so borders layer correctly
  for (j in rev(seq_len(ns))) {
    polygon(c(aj$time, rev(aj$time)),
            c(cum_upper[, j], rev(cum_lower[, j])),
            col = adjustcolor(col[j], alpha.f = 0.4),
            border = col[j], lwd = 1.5)
  }

  legend("topright", legend = states,
         fill = adjustcolor(col[seq_len(ns)], alpha.f = 0.4),
         border = col[seq_len(ns)], bty = "n", cex = 0.8)
}

#' Plot cumulative hazards from AJ
#' @noRd
.plot_cum_hazard_aj <- function(aj, states, col, main, xlab, ylab, ...) {
  trans_list <- .filter_plot_transitions(aj$structure$trans_list, states)
  n_trans <- nrow(trans_list)
  if (!n_trans) stop("No transitions lead to the selected state(s).")
  tcol <- rep_len(col, n_trans)

  # First pass: compute max cumulative hazard for y-axis
  max_ch <- 0
  ch_data <- vector("list", n_trans)
  for (tr in seq_len(n_trans)) {
    from <- trans_list$from[tr]
    to <- trans_list$to[tr]
    fi <- match(from, aj$structure$state_names)
    ti <- match(to, aj$structure$state_names)
    ch_vals <- vapply(aj$cum_hazard, function(M) M[fi, ti], numeric(1))
    ch_data[[tr]] <- ch_vals
    max_ch <- max(max_ch, max(ch_vals, na.rm = TRUE))
  }

  if (max_ch <= 0) max_ch <- 1

  plot(NULL, xlim = range(aj$time), ylim = c(0, max_ch * 1.1),
       xlab = xlab, ylab = ylab,
       main = main, ...)

  for (tr in seq_len(n_trans)) {
    lines(aj$time, ch_data[[tr]], col = tcol[tr], lwd = 2)
  }

  labels <- paste(trans_list$from, "->", trans_list$to)
  legend("topleft", legend = labels, col = tcol[seq_len(n_trans)],
         lwd = 2, bty = "n", cex = 0.7)
}

#' Plot Nelson-Aalen hazard increments from AJ
#' @noRd
.plot_hazard_increment_aj <- function(aj, states, col, main, xlab, ylab, ...) {
  state_names <- aj$structure$state_names
  trans_list <- aj$structure$trans_list

  # Keep transitions whose destination is in selected states
  trans_sub <- .filter_plot_transitions(trans_list, states)
  n_sub <- nrow(trans_sub)

  if (n_sub == 0) {
    stop("No transitions lead to the selected state(s).")
  }

  # Extract hazard increments for each transition
  int_data <- vector("list", n_sub)
  max_int <- 0
  for (tr in seq_len(n_sub)) {
    fi <- match(trans_sub$from[tr], state_names)
    ti <- match(trans_sub$to[tr], state_names)
    vals <- vapply(aj$hazard_inc, function(M) M[fi, ti], numeric(1))
    int_data[[tr]] <- vals
    max_int <- max(max_int, max(vals, na.rm = TRUE))
  }

  if (max_int <= 0) max_int <- 0.1

  plot(NULL, xlim = range(aj$time), ylim = c(0, max_int * 1.1),
       xlab = xlab, ylab = ylab,
       main = main, ...)

  # Point shape by origin state, color by destination state
  origin_states <- unique(trans_sub$from)
  n_origins <- length(origin_states)
  pch_base <- c(16, 17, 15, 18, 8, 3, 4, 1, 2, 0)
  pch_set <- rep_len(pch_base, n_origins)
  pch_map <- integer(n_sub)
  for (tr in seq_len(n_sub)) {
    pch_map[tr] <- pch_set[match(trans_sub$from[tr], origin_states)]
  }

  for (tr in seq_len(n_sub)) {
    col_idx <- match(trans_sub$to[tr], states)
    nz <- int_data[[tr]] > 0
    if (any(nz)) {
      points(aj$time[nz], int_data[[tr]][nz], col = col[col_idx],
             pch = pch_map[tr], cex = 0.8)
    }
  }

  labels <- paste(trans_sub$from, "->", trans_sub$to)
  legend_cols <- col[match(trans_sub$to, states)]
  legend("topright", legend = labels, col = legend_cols,
         pch = pch_map, pt.cex = 0.8, bty = "n", cex = 0.7)
}

#' Plot predicted state occupation
#' @noRd
.plot_state_occ_pred <- function(times, occ, state_names, states, col, main,
                                 xlab, ylab, ...) {
  state_idx <- match(states, state_names)

  plot(NULL, xlim = range(times), ylim = c(0, 1),
       xlab = xlab, ylab = ylab,
       main = main, ...)

  for (i in seq_along(state_idx)) {
    lines(times, occ[state_idx[i], ], col = col[state_idx[i]], lwd = 2)
  }

  legend("topright", legend = states, col = col[state_idx],
         lwd = 2, bty = "n", cex = 0.8)
}

#' Plot predicted transition probabilities
#' @noRd
.plot_trans_prob_pred <- function(times, P, state_names, states, start_state,
                                  col, main, xlab, ylab, ...) {
  state_idx <- match(states, state_names)
  start_idx <- match(start_state, dimnames(P)[[1L]])
  if (is.na(start_idx)) {
    stop("The prediction array does not contain the recorded start state.")
  }

  plot(NULL, xlim = range(times), ylim = c(0, 1),
       xlab = xlab, ylab = ylab,
       main = main, ...)

  for (j in seq_along(state_idx)) {
    probs <- P[start_idx, state_idx[j], ]
    lines(times, probs, col = col[state_idx[j]], lwd = 2)
  }

  legend("topright", legend = paste(start_state, "->", states),
         col = col[state_idx], lwd = 2, bty = "n", cex = 0.8)
}

#' Plot importance barplot
#' @noRd
.plot_importance_bar <- function(imp, col, main, xlab, ylab, ...) {
  mat <- imp$importance_matrix
  mat[is.na(mat)] <- 0
  n_vars <- nrow(mat)

  if (is.null(col)) col <- .default_palette(n_vars)
  if (!length(col)) stop("'col' must contain at least one color.")
  col <- rep_len(col, n_vars)

  # mat: rows = features, cols = transitions
  # barplot(mat, beside=TRUE) groups by columns (transitions),
  # bars within each group = rows (features)
  op <- par(mar = c(4, 10, 3, 8), xpd = TRUE)
  on.exit(par(op))

  barplot(mat, beside = TRUE, col = col[seq_len(n_vars)],
          main = main, xlab = xlab, ylab = ylab, horiz = TRUE,
          las = 1, cex.names = 0.7, ...)

  legend("topright", inset = c(-0.2, 0),
         legend = rownames(mat), fill = col[seq_len(n_vars)],
         bty = "n", cex = 0.65)
}

#' Plot importance heatmap
#' @noRd
.plot_importance_heat <- function(imp, col, main, xlab, ylab, ...) {
  mat <- imp$importance_matrix
  mat[is.na(mat)] <- 0

  nr <- nrow(mat)
  nc <- ncol(mat)

  # Color gradient
  n_cols <- 100
  if (is.null(col)) {
    col <- c("white", "#fee0d2", "#fc9272", "#de2d26")
  }
  if (!length(col)) stop("'col' must contain at least one color.")
  col_pal <- if (length(col) == 1L) {
    rep(col, n_cols)
  } else {
    grDevices::colorRampPalette(col)(n_cols)
  }

  # Normalize
  mat_norm <- (mat - min(mat)) / (max(mat) - min(mat) + 1e-10)

  op <- par(mar = c(8, 6, 3, 2))
  on.exit(par(op))

  image(seq_len(nc), seq_len(nr), t(mat_norm),
        col = col_pal, axes = FALSE,
        xlab = xlab, ylab = ylab, main = main, ...)

  axis(1, at = seq_len(nc), labels = colnames(mat),
       las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nr), labels = rownames(mat),
       las = 1, cex.axis = 0.8)

  # Add text values
  for (i in seq_len(nc)) {
    for (j in seq_len(nr)) {
      text(i, j, round(mat[j, i], 3), cex = 0.6)
    }
  }
  box()
}

#' Plot Brier scores
#' @noRd
.plot_brier <- function(diag, col, main, xlab, ylab, ...) {
  brier <- diag$brier
  if (is.null(brier) || !nrow(brier)) {
    stop("No patient-level cross-validated Brier score is present; rerun ",
         "diagnose(..., method = 'cv').")
  }
  if (is.null(col)) col <- "#2c7fb8"
  if (!length(col)) stop("'col' must contain at least one color.")
  plot(brier$time, brier$brier, type = "l", lwd = 2, col = col[1L],
       xlab = xlab, ylab = ylab, main = main, ...)
}

#' Plot concordance indices
#' @noRd
.plot_concordance <- function(diag, col, main, xlab, ylab, ...) {
  cdf <- diag$concordance
  n <- nrow(cdf)
  if (is.null(col)) col <- .default_palette(n)
  if (!length(col)) stop("'col' must contain at least one color.")
  col <- rep_len(col, n)

  op <- par(mar = c(8, 4, 3, 2))
  on.exit(par(op))

  bp <- barplot(cdf$c_index, names.arg = cdf$transition,
                col = col[seq_len(n)], main = main,
                xlab = xlab, ylab = ylab, ylim = c(0, 1),
                las = 2, cex.names = 0.7, ...)
  abline(h = 0.5, lty = 2, col = "gray50")
  text(bp, cdf$c_index + 0.03, round(cdf$c_index, 3), cex = 0.8)
}

#' Plot Transition Diagram
#'
#' Draws a state transition diagram with event counts annotated on edges.
#' Uses a layered layout that adapts to any number of states and
#' automatically routes arrows around intermediate state boxes using
#' Bezier curves when needed.
#'
#' @param structure An \code{mstate_structure} object.
#' @param msdata Optional \code{msdata} object to annotate with counts.
#' @param col Node colors. Default uses the standard palette.
#' @param main Title.
#' @param xlab,ylab Optional axis labels; defaults are blank because graph
#'   coordinates have no statistical scale.
#' @param ... Ignored.
#'
#' @return No return value, called for its side effect of producing a plot.
#'
#' @section Limitations:
#' The layout is descriptive and preserves the validated DAG. It does not
#' represent cycles/recurrent visits, transition intensity, uncertainty, or
#' edge direction over calendar time. Counts are shown only when a compatible
#' \code{msdata} object is supplied.
#'
#' @examples
#' ms <- clinical_states()
#' plot_transition_diagram(ms, xlab = "Graph layer")
#'
#' @export
plot_transition_diagram <- function(structure, msdata = NULL,
                                    col = NULL, main = "Transition Diagram",
                                    xlab = "", ylab = "", ...) {
  if (!inherits(structure, "mstate_structure")) {
    stop("'structure' must be an mstate_structure object.")
  }
  if (!is.null(msdata) && !inherits(msdata, "msdata")) {
    stop("'msdata' must be NULL or an msdata object.")
  }
  state_names <- structure$state_names
  ns <- structure$n_states
  trans_list <- structure$trans_list

  if (is.null(col)) col <- .default_palette(ns)
  if (!length(col)) stop("'col' must contain at least one color.")
  col <- rep_len(col, ns)

  # Compute layered layout
  layout <- .layout_states(state_names, structure)
  pos_x <- layout$x
  pos_y <- layout$y
  state_layer <- layout$layer

  # Box dimensions scaled to longest state name
  max_chars <- max(nchar(state_names))
  box_w <- max(0.4, max_chars * 0.055 + 0.15)
  box_h <- 0.25

  # Plot area with padding
  x_pad <- box_w + 1.0
  y_pad <- max(box_h + 1.0, 1.5)

  op <- par(mar = c(1, 1, 3, 1))
  on.exit(par(op))

  plot(NULL,
       xlim = c(min(pos_x) - x_pad, max(pos_x) + x_pad),
       ylim = c(min(pos_y) - y_pad, max(pos_y) + y_pad),
       xlab = xlab, ylab = ylab, main = main, axes = FALSE, asp = 1)

  # Color palette for transitions
  n_trans <- nrow(trans_list)
  trans_cols <- .transition_palette(n_trans)

  # Draw arrows first (below boxes)
  for (tr in seq_len(n_trans)) {
    from_name <- trans_list$from[tr]
    to_name <- trans_list$to[tr]
    fi <- match(from_name, state_names)
    ti <- match(to_name, state_names)

    x0 <- pos_x[fi]; y0 <- pos_y[fi]
    x1 <- pos_x[ti]; y1 <- pos_y[ti]

    arrow_col <- trans_cols[tr]

    # Event count label
    n_ev <- NULL
    if (!is.null(msdata)) {
      n_ev <- sum(msdata$status == 1 & msdata$from == from_name &
                    msdata$to == to_name, na.rm = TRUE)
    }

    # Find states in strictly intermediate layers
    from_l <- state_layer[fi]
    to_l <- state_layer[ti]
    if (abs(to_l - from_l) > 1) {
      min_l <- min(from_l, to_l) + 1L
      max_l <- max(from_l, to_l) - 1L
      between_idx <- which(state_layer >= min_l & state_layer <= max_l)
      between_idx <- setdiff(between_idx, c(fi, ti))
    } else {
      between_idx <- integer(0)
    }

    # Check if any box blocks the straight path
    other_idx <- setdiff(seq_len(ns), c(fi, ti))
    blocked <- any(vapply(other_idx, function(bi) {
      .line_near_box(x0, y0, x1, y1, pos_x[bi], pos_y[bi], box_w, box_h)
    }, logical(1)))

    if (blocked && length(between_idx) > 0) {
      # Curve around intermediate boxes
      between_y <- pos_y[between_idx]
      y_mid <- (y0 + y1) / 2

      if (y0 >= y_mid) {
        curve_dir <- 1  # above
        max_y <- max(between_y + box_h)
        top_ref <- max(y0, y1)
        clearance <- max(max_y - top_ref, 0) + box_h + 0.4
      } else {
        curve_dir <- -1  # below
        min_y <- min(between_y - box_h)
        bot_ref <- min(y0, y1)
        clearance <- max(bot_ref - min_y, 0) + box_h + 0.4
      }

      offset <- curve_dir * clearance * 2
      .draw_curved_arrow(x0, y0, x1, y1, box_w, box_h, offset, n_ev,
                         arrow_col)

    } else if (blocked) {
      # Blocked by adjacent-layer box: small curve
      block_sides <- vapply(other_idx, function(bi) {
        if (.line_near_box(x0, y0, x1, y1,
                           pos_x[bi], pos_y[bi], box_w, box_h)) {
          .perpendicular_side(x0, y0, x1, y1, pos_x[bi], pos_y[bi])
        } else {
          0
        }
      }, numeric(1))
      side_sum <- sum(sign(block_sides))
      curve_dir <- if (side_sum >= 0) -1 else 1
      offset <- curve_dir * box_h * 4
      .draw_curved_arrow(x0, y0, x1, y1, box_w, box_h, offset, n_ev,
                         arrow_col)

    } else {
      .draw_straight_arrow(x0, y0, x1, y1, box_w, box_h, n_ev, arrow_col)
    }
  }

  # Legend mapping colors to transitions
  trans_labels <- paste(trans_list$from, "->", trans_list$to)
  legend("bottomright", legend = trans_labels, col = trans_cols,
         lwd = 2, bty = "n", cex = 0.55)

  # Draw state boxes on top of arrows
  for (i in seq_len(ns)) {
    is_absorb <- state_names[i] %in% structure$absorbing
    border_col <- if (is_absorb) "red3" else "black"
    border_lwd <- if (is_absorb) 2.5 else 1.5

    rect(pos_x[i] - box_w, pos_y[i] - box_h,
         pos_x[i] + box_w, pos_y[i] + box_h,
         col = adjustcolor(col[i], alpha.f = 0.3),
         border = border_col, lwd = border_lwd)
    text(pos_x[i], pos_y[i], state_names[i], cex = 0.7, font = 2)
  }
  invisible(NULL)
}

#' Layered layout for state diagram
#' @noRd
.layout_states <- function(state_names, structure) {
  ns <- length(state_names)

  # Longest-path layer assignment for transient states
  layer <- integer(ns)
  names(layer) <- state_names

  for (iter in seq_len(ns)) {
    updated <- FALSE
    for (from in names(structure$transitions)) {
      for (to in structure$transitions[[from]]) {
        if (!(to %in% structure$absorbing)) {
          if (layer[from] + 1L > layer[to]) {
            layer[to] <- layer[from] + 1L
            updated <- TRUE
          }
        }
      }
    }
    if (!updated) break
  }

  # Absorbing states at max_transient + 1
  max_t <- if (length(structure$transient) > 0) {
    max(layer[structure$transient])
  } else {
    0L
  }
  for (a in structure$absorbing) {
    layer[a] <- max_t + 1L
  }

  # Position states: x from layer, y spread within layer
  x_spacing <- 2.0
  y_spacing <- 1.5

  x <- numeric(ns)
  y <- numeric(ns)
  names(x) <- state_names
  names(y) <- state_names

  for (l in sort(unique(layer))) {
    states_in_l <- state_names[layer[state_names] == l]
    n_l <- length(states_in_l)
    x[states_in_l] <- l * x_spacing
    if (n_l == 1L) {
      y[states_in_l] <- 0
    } else {
      offsets <- seq(-(n_l - 1) / 2, (n_l - 1) / 2, length.out = n_l)
      y[states_in_l] <- offsets * y_spacing
    }
  }

  data.frame(x = x, y = y, layer = layer)
}

#' Compute point on box edge toward a target
#' @noRd
.box_edge_point <- function(cx, cy, bw, bh, tx, ty) {
  dx <- tx - cx
  dy <- ty - cy
  if (abs(dx) < 1e-10 && abs(dy) < 1e-10) return(c(cx + bw, cy))
  if (abs(dx) < 1e-10) return(c(cx, cy + sign(dy) * bh))
  if (abs(dy) < 1e-10) return(c(cx + sign(dx) * bw, cy))
  t_x <- bw / abs(dx)
  t_y <- bh / abs(dy)
  t <- min(t_x, t_y)
  c(cx + dx * t, cy + dy * t)
}

#' Check if line segment passes near a box
#' @noRd
.line_near_box <- function(x0, y0, x1, y1, bx, by, bw, bh) {
  dx <- x1 - x0; dy <- y1 - y0
  len_sq <- dx * dx + dy * dy
  if (len_sq < 1e-10) return(FALSE)
  t <- ((bx - x0) * dx + (by - y0) * dy) / len_sq
  t <- max(0.05, min(0.95, t))
  cx <- x0 + t * dx
  cy <- y0 + t * dy
  abs(cx - bx) < (bw + 0.15) && abs(cy - by) < (bh + 0.15)
}

#' Which side of a directed line a point falls on
#' @noRd
.perpendicular_side <- function(x0, y0, x1, y1, px, py) {
  (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0)
}

#' Color palette for transition arrows
#' @noRd
.transition_palette <- function(n) {
  base <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
            "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
            "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494",
            "#b3b3b3")
  if (n <= length(base)) {
    base[seq_len(n)]
  } else {
    grDevices::rainbow(n, s = 0.75, v = 0.8)
  }
}

#' Draw a straight arrow between two box centers with edge clipping
#' @noRd
.draw_straight_arrow <- function(x0, y0, x1, y1, bw, bh, n_ev,
                                 col = "gray30") {
  start <- .box_edge_point(x0, y0, bw, bh, x1, y1)
  end   <- .box_edge_point(x1, y1, bw, bh, x0, y0)
  arrows(start[1], start[2], end[1], end[2],
         length = 0.1, lwd = 1.5, col = col)

  if (!is.null(n_ev)) {
    lx <- start[1] + (end[1] - start[1]) / 3
    ly <- start[2] + (end[2] - start[2]) / 3
    dx <- end[1] - start[1]; dy <- end[2] - start[2]
    len <- sqrt(dx * dx + dy * dy)
    if (len > 0) {
      lx <- lx + (-dy / len) * 0.18
      ly <- ly + ( dx / len) * 0.18
    }
    text(lx, ly, n_ev, cex = 0.65, col = col, font = 2)
  }
}

#' Draw a curved (Bezier) arrow between two box centers
#' @param offset Signed perpendicular offset for the control point.
#'   Positive = left of arrow direction (usually up for left-to-right).
#' @noRd
.draw_curved_arrow <- function(x0, y0, x1, y1, bw, bh, offset, n_ev,
                               col = "gray30") {
  dx <- x1 - x0; dy <- y1 - y0
  len <- sqrt(dx * dx + dy * dy)
  if (len < 1e-10) return()

  # Perpendicular unit vector (left of direction)
  px <- -dy / len; py <- dx / len

  # Control point at midpoint + perpendicular offset
  mx <- (x0 + x1) / 2 + px * offset
  my <- (y0 + y1) / 2 + py * offset

  # Start/end at box edges toward control point
  start <- .box_edge_point(x0, y0, bw, bh, mx, my)
  end   <- .box_edge_point(x1, y1, bw, bh, mx, my)

  # Quadratic Bezier curve
  t_vals <- seq(0, 1, length.out = 60)
  cx <- (1 - t_vals)^2 * start[1] +
    2 * (1 - t_vals) * t_vals * mx +
    t_vals^2 * end[1]
  cy <- (1 - t_vals)^2 * start[2] +
    2 * (1 - t_vals) * t_vals * my +
    t_vals^2 * end[2]

  lines(cx, cy, lwd = 1.5, col = col)

  # Arrowhead at the end
  n <- length(cx)
  arrows(cx[n - 3], cy[n - 3], cx[n], cy[n],
         length = 0.1, lwd = 1.5, col = col)

  # Label at 1/3 of curve with perpendicular offset
  if (!is.null(n_ev)) {
    tl <- 1 / 3
    lx <- (1 - tl)^2 * start[1] + 2 * (1 - tl) * tl * mx + tl^2 * end[1]
    ly <- (1 - tl)^2 * start[2] + 2 * (1 - tl) * tl * my + tl^2 * end[2]
    # Tangent at t for perpendicular label offset
    tx <- 2 * (1 - tl) * (mx - start[1]) + 2 * tl * (end[1] - mx)
    ty <- 2 * (1 - tl) * (my - start[2]) + 2 * tl * (end[2] - my)
    tlen <- sqrt(tx * tx + ty * ty)
    if (tlen > 0) {
      lx <- lx + (-ty / tlen) * 0.18
      ly <- ly + ( tx / tlen) * 0.18
    }
    text(lx, ly, n_ev, cex = 0.65, col = col, font = 2)
  }
}
