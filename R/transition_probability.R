#' Compute Transition Probability Matrix via Product-Integral
#'
#' Given cause-specific cumulative hazard functions for all transitions,
#' computes the full transition probability matrix P(s,t) using the
#' product-integral formula.
#'
#' @param cum_hazards A list of cumulative hazard data frames, one per
#'   transition. Each should have columns \code{time} and \code{hazard}.
#' @param structure An \code{mstate_structure} object.
#' @param s Numeric, starting time (default 0).
#' @param times Numeric vector of times at which to evaluate P(s,t). If
#'   \code{NULL}, uses all unique event times from the cumulative hazards.
#'
#' @return An object of class \code{"trans_prob"} containing:
#'   \describe{
#'     \item{time}{Evaluation times.}
#'     \item{P}{List of transition probability matrices at each time.}
#'     \item{state_occ}{Matrix of state occupation probabilities.}
#'     \item{structure}{The multistate structure.}
#'     \item{s}{Starting time.}
#'   }
#'
#' @export
compute_trans_prob <- function(cum_hazards, structure, s = 0, times = NULL) {
  state_names <- structure$state_names
  ns <- structure$n_states
  trans_list <- structure$trans_list

  # Collect all unique time points
  all_times <- sort(unique(unlist(
    lapply(cum_hazards, function(ch) ch$time[ch$time > s])
  )))

  if (!is.null(times)) {
    all_times <- sort(unique(c(all_times, times)))
    all_times <- all_times[all_times > s]
  }

  if (length(all_times) == 0) {
    stop("No time points found after s = ", s)
  }

  # Interpolate cumulative hazards to all time points and compute increments
  # For each transition, get cumulative hazard at each time point
  ch_interp <- lapply(cum_hazards, function(ch) {
    stats::approx(ch$time, ch$hazard, xout = all_times,
                  method = "constant", rule = 2, f = 0)$y
  })

  # Compute hazard increments (differences)
  ch_inc <- lapply(ch_interp, function(h) {
    c(h[1], diff(h))
  })

  # Product-integral
  P_list <- vector("list", length(all_times))
  P_current <- diag(ns)
  dimnames(P_current) <- list(state_names, state_names)

  state_occ <- matrix(0, nrow = length(all_times), ncol = ns,
                       dimnames = list(NULL, state_names))

  p0 <- rep(0, ns)
  p0[1] <- 1

  for (k in seq_along(all_times)) {
    # Build increment matrix dA at this time
    dA <- matrix(0, nrow = ns, ncol = ns,
                 dimnames = list(state_names, state_names))

    for (tr_idx in seq_len(nrow(trans_list))) {
      from <- trans_list$from[tr_idx]
      to <- trans_list$to[tr_idx]
      fi <- match(from, state_names)
      ti <- match(to, state_names)
      if (!is.null(ch_inc[[tr_idx]])) {
        dA[fi, ti] <- ch_inc[[tr_idx]][k]
      }
    }

    # Set diagonal
    for (h in seq_len(ns)) {
      dA[h, h] <- -sum(dA[h, -h])
    }

    # Product-integral step
    increment <- diag(ns) + dA
    P_current <- P_current %*% increment
    P_current[P_current < 0] <- 0

    # Normalize rows to sum to 1
    row_sums <- rowSums(P_current)
    for (h in seq_len(ns)) {
      if (row_sums[h] > 0) {
        P_current[h, ] <- P_current[h, ] / row_sums[h]
      }
    }

    P_list[[k]] <- P_current
    state_occ[k, ] <- as.vector(p0 %*% P_current)
  }

  # If specific times requested, subset
  if (!is.null(times)) {
    keep <- all_times %in% times
    if (sum(keep) > 0) {
      # Use nearest available times
      eval_idx <- vapply(times, function(t) {
        which.min(abs(all_times - t))
      }, integer(1))
      eval_idx <- unique(eval_idx)
    } else {
      eval_idx <- seq_along(all_times)
    }
  } else {
    eval_idx <- seq_along(all_times)
  }

  structure(
    list(
      time = all_times[eval_idx],
      P = P_list[eval_idx],
      state_occ = state_occ[eval_idx, , drop = FALSE],
      structure = structure,
      s = s
    ),
    class = "trans_prob"
  )
}

#' @export
print.trans_prob <- function(x, ...) {
  cat("Transition Probabilities\n")
  cat("  Time range: [", min(x$time), ", ", max(x$time), "]\n", sep = "")
  cat("  Time points:", length(x$time), "\n")
  cat("  States:", paste(x$structure$state_names, collapse = ", "), "\n")
  cat("\nFinal P(s, t_max):\n")
  print(round(x$P[[length(x$P)]], 4))
  invisible(x)
}
