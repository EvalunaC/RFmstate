#' Entry-Conditioned Semi-Markov State Probabilities
#'
#' Combines one clock-reset cumulative cause-specific hazard curve per allowed
#' edge using a regular-grid entry-mass/sojourn convolution. The result is not
#' a general Markov \eqn{P(s,t)} matrix: each row is conditional on fresh entry
#' into its named starting state at elapsed duration zero.
#'
#' @param cum_hazards Named list of data frames, one for every allowed edge.
#'   Names must be \code{"from->to"}; each frame contains increasing
#'   \code{time} and nondecreasing \code{hazard} (or \code{cum_hazard}).
#' @param structure An \code{mstate_structure} object.
#' @param s Legacy starting-time argument. Only zero is supported.
#' @param times Finite nonnegative elapsed durations. If \code{NULL}, the union
#'   of input curve times is used.
#' @param start_state State in which fresh entry occurs. The default is the
#'   common initial state stored in \code{structure}.
#' @param grid_step Optional initial internal grid step.
#' @param target_grid_points Initial number of regular-grid intervals when
#'   \code{grid_step} is not supplied.
#' @param max_grid_points Maximum number of internal grid points.
#' @param grid_tol Maximum probability change allowed on grid refinement.
#'   The default \code{5e-4} reflects discontinuous ranger hazard curves;
#'   analytic tests retain their separate \code{1e-3} truth-error gate.
#' @param max_grid_refinements Maximum number of step halvings.
#' @param check_grid Whether to require grid-refinement convergence.
#' @param prob_tol Probability invariant tolerance.
#' @param extrapolate Either \code{"error"} (default) or explicit
#'   \code{"flat"} cumulative-hazard extension.
#'
#' @return A \code{trans_prob} object. \code{entry_prob} and its temporary
#'   alias \code{P} are arrays ordered starting state by occupied state by
#'   elapsed time. \code{state_occ} is the selected starting-state slice.
#'
#' @details Cumulative hazards are evaluated as right-continuous step
#' functions. Within each regular interval, total exit probability is allocated
#' to causes in proportion to their cumulative-hazard increments. Entry mass is
#' assigned to the interval's right endpoint. Grid refinement controls these
#' approximations. Probabilities are never clipped or row-normalized.
#'
#' @examples
#' ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
#' tt <- seq(0, 2, length.out = 2001)
#' hazards <- list("A->B" = data.frame(time = tt, hazard = 0.4 * tt))
#' pr <- compute_trans_prob(hazards, ms, times = c(0, 1, 2),
#'                          target_grid_points = 512)
#' pr$state_occ
#'
#' @export
compute_trans_prob <- function(cum_hazards, structure, s = 0, times = NULL,
                               start_state = NULL, grid_step = NULL,
                               target_grid_points = 1024L,
                               max_grid_points = 131073L,
                               grid_tol = 5e-4,
                               max_grid_refinements = 7L,
                               check_grid = TRUE, prob_tol = 1e-8,
                               extrapolate = c("error", "flat")) {
  if (!inherits(structure, "mstate_structure")) {
    stop("'structure' must be an 'mstate_structure' object.")
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s != 0) {
    stop("Only s = 0 fresh-entry prediction is supported; delayed-entry and ",
         "ongoing-sojourn prediction are unavailable.")
  }
  if (is.null(start_state)) start_state <- structure$initial_state
  if (!is.character(start_state) || length(start_state) != 1L ||
      is.na(start_state) || !(start_state %in% structure$state_names)) {
    stop("'start_state' must be one state in the structure.")
  }
  extrapolate <- match.arg(extrapolate)
  curves <- .validate_hazard_curves(cum_hazards, structure)

  if (is.null(times)) {
    times <- sort(unique(c(0, unlist(lapply(curves, `[[`, "time")))))
  }
  if (!is.numeric(times) || !length(times) || anyNA(times) ||
      any(!is.finite(times)) || any(times < 0)) {
    stop("'times' must contain finite, nonnegative elapsed durations.")
  }
  times <- sort(unique(as.numeric(times)))

  support <- .curve_support_from_start(curves, structure, start_state)
  tau <- max(times)
  if (is.finite(support) && tau > support) {
    if (extrapolate == "error") {
      stop("Requested horizon ", format(tau), " exceeds conservative support ",
           format(support), " from starting state '", start_state,
           "'. Use extrapolate = 'flat' only for an explicit sensitivity analysis.")
    }
    warning("Using flat cumulative-hazard extension beyond observed support; ",
            "this imposes zero additional exit hazard and is unsuitable for ",
            "primary validation.", call. = FALSE)
  }

  solved <- .solve_semi_markov(
    curves = curves,
    structure = structure,
    requested_times = times,
    grid_step = grid_step,
    target_grid_points = target_grid_points,
    max_grid_points = max_grid_points,
    grid_tol = grid_tol,
    max_grid_refinements = max_grid_refinements,
    check_grid = check_grid,
    prob_tol = prob_tol
  )
  selected <- solved$prob[start_state, , , drop = FALSE]
  state_occ <- t(selected[1L, , , drop = TRUE])
  if (length(times) == 1L) {
    state_occ <- matrix(selected[1L, , 1L], nrow = 1L,
                        dimnames = list(NULL, structure$state_names))
  }

  structure(
    list(
      time = times,
      entry_prob = solved$prob,
      P = solved$prob,
      state_occ = state_occ,
      structure = structure,
      s = 0,
      start_state = start_state,
      time_scale = "clock-reset",
      conditioning = "fresh entry at elapsed duration zero",
      support_horizon = support,
      extrapolation = if (tau > support) "flat cumulative hazard" else "none",
      grid = solved$grid,
      grid_step = solved$grid_step,
      grid_points = length(solved$grid),
      grid_converged = solved$converged,
      grid_error = solved$error,
      grid_refinements = solved$refinements,
      hazard_roundoff_corrections = solved$hazard_corrections,
      prob_tol = prob_tol
    ),
    class = "trans_prob"
  )
}

#' @noRd
.validate_hazard_curves <- function(cum_hazards, structure) {
  if (!is.list(cum_hazards) || !length(cum_hazards) ||
      is.null(names(cum_hazards)) || anyNA(names(cum_hazards)) ||
      any(!nzchar(names(cum_hazards))) || anyDuplicated(names(cum_hazards))) {
    stop("'cum_hazards' must be a uniquely named list with names 'from->to'.")
  }
  expected <- paste0(structure$trans_list$from, "->", structure$trans_list$to)
  missing_edges <- setdiff(expected, names(cum_hazards))
  unknown_edges <- setdiff(names(cum_hazards), expected)
  if (length(missing_edges) || length(unknown_edges)) {
    stop("Cumulative-hazard edges do not match the declared graph. Missing: ",
         if (length(missing_edges)) paste(missing_edges, collapse = ", ") else "none",
         "; unknown: ",
         if (length(unknown_edges)) paste(unknown_edges, collapse = ", ") else "none",
         ".")
  }
  curves <- cum_hazards[expected]
  for (edge in expected) {
    curve <- curves[[edge]]
    if (!is.data.frame(curve) || !("time" %in% names(curve))) {
      stop("Hazard curve '", edge, "' must be a data frame with a 'time' column.")
    }
    value_name <- if ("cum_hazard" %in% names(curve)) {
      "cum_hazard"
    } else if ("hazard" %in% names(curve)) {
      "hazard"
    } else {
      stop("Hazard curve '", edge,
           "' must contain 'hazard' or 'cum_hazard'.")
    }
    tt <- curve$time
    aa <- curve[[value_name]]
    if (!is.numeric(tt) || !is.numeric(aa) || length(tt) != length(aa) ||
        !length(tt) || anyNA(tt) || anyNA(aa) || any(!is.finite(tt)) ||
        any(!is.finite(aa)) || any(tt < 0) || any(aa < 0)) {
      stop("Hazard curve '", edge,
           "' must have finite nonnegative times and cumulative hazards.")
    }
    if (is.unsorted(tt, strictly = TRUE)) {
      stop("Hazard curve '", edge, "' times must be strictly increasing.")
    }
    delta <- diff(aa)
    tol <- pmax(1e-10, 1e-8 * pmax(1, head(aa, -1L), tail(aa, -1L)))
    if (any(delta < -tol)) {
      k <- which(delta < -tol)[1L]
      stop("Cumulative hazard decreases materially for edge '", edge,
           "' between times ", tt[k], " and ", tt[k + 1L], ".")
    }
    if (length(delta) && any(delta < 0)) {
      for (k in which(delta < 0)) aa[k + 1L] <- aa[k]
    }
    curves[[edge]] <- data.frame(time = tt, hazard = aa)
  }
  curves
}

#' @noRd
.reachable_states <- function(structure, start_state) {
  reached <- start_state
  frontier <- start_state
  while (length(frontier)) {
    children <- unique(structure$trans_list$to[
      structure$trans_list$from %in% frontier
    ])
    new <- setdiff(children, reached)
    reached <- c(reached, new)
    frontier <- new
  }
  reached
}

#' @noRd
.curve_support_from_start <- function(curves, structure, start_state) {
  if (start_state %in% structure$absorbing) return(Inf)
  reachable <- intersect(.reachable_states(structure, start_state),
                         structure$transient)
  origin_support <- vapply(reachable, function(state) {
    edges <- paste0(state, "->", structure$transitions[[state]])
    min(vapply(curves[edges], function(curve) max(curve$time), numeric(1)))
  }, numeric(1))
  min(origin_support)
}

#' @noRd
.step_cumulative_hazard <- function(curve, grid) {
  idx <- findInterval(grid, curve$time)
  out <- numeric(length(grid))
  keep <- idx > 0L
  out[keep] <- curve$hazard[idx[keep]]
  out
}

#' @noRd
.linear_convolution_fft <- function(x, y, n_out) {
  full_length <- length(x) + length(y) - 1L
  n_fft <- 2L^ceiling(log2(full_length))
  xx <- c(x, numeric(n_fft - length(x)))
  yy <- c(y, numeric(n_fft - length(y)))
  Re(stats::fft(stats::fft(xx) * stats::fft(yy), inverse = TRUE) / n_fft)[
    seq_len(n_out)
  ]
}

#' Direct reference convolution retained for regression tests
#' @noRd
.linear_convolution_direct <- function(x, y, n_out) {
  out <- numeric(n_out)
  for (k in seq_len(n_out)) {
    out[k] <- sum(x[seq_len(k)] * y[k:1L])
  }
  out
}

#' @noRd
.semi_markov_once <- function(curves, structure, grid, prob_tol = 1e-8,
                              hazard_tol_abs = 1e-10,
                              hazard_tol_rel = 1e-8,
                              engine = c("fft", "direct")) {
  engine <- match.arg(engine)
  conv <- if (engine == "fft") .linear_convolution_fft else
    .linear_convolution_direct
  state_names <- structure$state_names
  ns <- length(state_names)
  trans_list <- structure$trans_list
  n_edges <- nrow(trans_list)
  n_grid <- length(grid)
  n_intervals <- n_grid - 1L
  edge_names <- paste0(trans_list$from, "->", trans_list$to)

  survival <- stats::setNames(vector("list", ns), state_names)
  dF <- stats::setNames(vector("list", n_edges), edge_names)
  corrections <- 0L
  for (state in state_names) {
    if (state %in% structure$absorbing) {
      survival[[state]] <- rep(1, n_grid)
      next
    }
    edge_idx <- which(trans_list$from == state)
    increments <- matrix(0, nrow = length(edge_idx), ncol = n_intervals)
    for (r in seq_along(edge_idx)) {
      edge <- edge_names[edge_idx[r]]
      aa <- .step_cumulative_hazard(curves[[edge]], grid)
      inc <- diff(aa)
      tol <- pmax(
        hazard_tol_abs,
        hazard_tol_rel * pmax(1, head(aa, -1L), tail(aa, -1L))
      )
      if (any(inc < -tol)) {
        k <- which(inc < -tol)[1L]
        stop("Cumulative hazard decreases for edge '", edge,
             "' on grid interval (", grid[k], ", ", grid[k + 1L], "].")
      }
      tiny <- inc < 0
      corrections <- corrections + sum(tiny)
      inc[tiny] <- 0
      increments[r, ] <- inc
    }
    total <- colSums(increments)
    surv <- c(1, cumprod(exp(-total)))
    survival[[state]] <- surv
    exit_prob <- 1 - exp(-total)
    for (r in seq_along(edge_idx)) {
      vals <- numeric(n_grid)
      positive <- total > 0
      vals[-1L][positive] <- surv[-n_grid][positive] * exit_prob[positive] *
        increments[r, positive] / total[positive]
      dF[[edge_names[edge_idx[r]]]] <- vals
    }
  }

  prob <- array(
    0,
    dim = c(ns, ns, n_grid),
    dimnames = list(
      starting_state = state_names,
      occupied_state = state_names,
      elapsed_time = format(grid, scientific = FALSE, trim = TRUE)
    )
  )
  for (start in state_names) {
    entry <- matrix(0, nrow = ns, ncol = n_grid,
                    dimnames = list(state_names, NULL))
    entry[start, 1L] <- 1
    for (state in structure$topological_order) {
      h <- match(state, state_names)
      prob[start, state, ] <- conv(entry[h, ], survival[[state]], n_grid)
      edge_idx <- which(trans_list$from == state)
      if (length(edge_idx)) {
        for (idx in edge_idx) {
          child <- trans_list$to[idx]
          edge <- edge_names[idx]
          entry[child, ] <- entry[child, ] + conv(entry[h, ], dF[[edge]], n_grid)
        }
      }
    }
  }
  .check_probability_array(prob, prob_tol, "internal grid")
  list(prob = prob, survival = survival, exit_mass = dF,
       corrections = corrections)
}

#' @noRd
.interpolate_probability_array <- function(prob, grid, times,
                                           prob_tol = 1e-8) {
  dims <- dim(prob)
  out <- array(
    0,
    dim = c(dims[1L], dims[2L], length(times)),
    dimnames = list(
      starting_state = dimnames(prob)[[1L]],
      occupied_state = dimnames(prob)[[2L]],
      elapsed_time = format(times, scientific = FALSE, trim = TRUE)
    )
  )
  for (a in seq_len(dims[1L])) {
    for (h in seq_len(dims[2L])) {
      out[a, h, ] <- stats::approx(
        grid, prob[a, h, ], xout = times, method = "linear", rule = 2
      )$y
    }
  }
  .check_probability_array(out, prob_tol, "requested times")
  out
}

#' @noRd
.check_probability_array <- function(prob, prob_tol, location) {
  if (any(!is.finite(prob))) {
    stop("Nonfinite semi-Markov probability on ", location, ".")
  }
  if (min(prob) < -prob_tol || max(prob) > 1 + prob_tol) {
    stop("Semi-Markov probability outside [0,1] tolerance on ", location,
         ": range [", format(min(prob)), ", ", format(max(prob)), "].")
  }
  mass <- apply(prob, c(1L, 3L), sum)
  max_error <- max(abs(mass - 1))
  if (max_error > prob_tol) {
    stop("Semi-Markov probability mass error on ", location, " is ",
         format(max_error), ", exceeding prob_tol = ", prob_tol, ".")
  }
  invisible(TRUE)
}

#' @noRd
.solve_semi_markov <- function(curves, structure, requested_times,
                               grid_step, target_grid_points,
                               max_grid_points, grid_tol,
                               max_grid_refinements, check_grid, prob_tol) {
  tau <- max(requested_times)
  if (tau == 0) {
    one <- .semi_markov_once(curves, structure, 0, prob_tol)
    return(list(prob = one$prob, grid = 0, grid_step = NA_real_,
                converged = TRUE, error = 0, refinements = 0L,
                hazard_corrections = one$corrections))
  }
  integer_control <- function(x, name, minimum) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum ||
        x != as.integer(x)) stop("'", name, "' must be an integer >= ", minimum, ".")
    as.integer(x)
  }
  target_grid_points <- integer_control(target_grid_points,
                                        "target_grid_points", 16L)
  max_grid_points <- integer_control(max_grid_points, "max_grid_points", 33L)
  max_grid_refinements <- integer_control(max_grid_refinements,
                                          "max_grid_refinements", 1L)
  if (!is.numeric(grid_tol) || length(grid_tol) != 1L ||
      !is.finite(grid_tol) || grid_tol <= 0) {
    stop("'grid_tol' must be one positive finite number.")
  }
  if (!is.numeric(prob_tol) || length(prob_tol) != 1L ||
      !is.finite(prob_tol) || prob_tol <= 0) {
    stop("'prob_tol' must be one positive finite number.")
  }
  if (!is.null(grid_step)) {
    if (!is.numeric(grid_step) || length(grid_step) != 1L ||
        !is.finite(grid_step) || grid_step <= 0) {
      stop("'grid_step' must be NULL or one positive finite number.")
    }
    intervals <- max(1L, ceiling(tau / grid_step))
  } else {
    intervals <- target_grid_points
  }
  if (intervals + 1L > max_grid_points) {
    stop("Initial internal grid exceeds 'max_grid_points'.")
  }

  grid <- seq(0, tau, length.out = intervals + 1L)
  coarse <- .semi_markov_once(curves, structure, grid, prob_tol)
  coarse_out <- .interpolate_probability_array(
    coarse$prob, grid, requested_times, prob_tol
  )
  if (!isTRUE(check_grid)) {
    return(list(prob = coarse_out, grid = grid, grid_step = grid[2L],
                converged = NA, error = NA_real_, refinements = 0L,
                hazard_corrections = coarse$corrections))
  }

  last_error <- Inf
  total_corrections <- coarse$corrections
  for (refinement in seq_len(max_grid_refinements)) {
    intervals <- intervals * 2L
    if (intervals + 1L > max_grid_points) break
    fine_grid <- seq(0, tau, length.out = intervals + 1L)
    fine <- .semi_markov_once(curves, structure, fine_grid, prob_tol)
    fine_out <- .interpolate_probability_array(
      fine$prob, fine_grid, requested_times, prob_tol
    )
    last_error <- max(abs(fine_out - coarse_out))
    total_corrections <- total_corrections + fine$corrections
    if (last_error <= grid_tol) {
      return(list(prob = fine_out, grid = fine_grid,
                  grid_step = fine_grid[2L], converged = TRUE,
                  error = last_error, refinements = refinement,
                  hazard_corrections = total_corrections))
    }
    grid <- fine_grid
    coarse_out <- fine_out
  }
  stop("Internal semi-Markov grid did not converge: maximum change ",
       format(last_error), " exceeds grid_tol = ", grid_tol,
       " before max_grid_points/max_grid_refinements.")
}

#' @export
print.trans_prob <- function(x, ...) {
  cat("Entry-Conditioned Semi-Markov State Probabilities\n")
  cat("  Starting state:", x$start_state, "\n")
  cat("  Time scale: clock-reset (fresh entry at duration zero)\n")
  cat("  Elapsed-time range: [", min(x$time), ", ", max(x$time), "]\n", sep = "")
  cat("  Internal grid points:", x$grid_points,
      "(converged:", x$grid_converged, ")\n")
  cat("\nFinal occupied-state probabilities:\n")
  print(round(x$state_occ[nrow(x$state_occ), ], 4))
  invisible(x)
}
