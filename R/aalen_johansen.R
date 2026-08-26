#' Aalen-Johansen Point Estimator
#'
#' Computes nonparametric estimates of transition probabilities using the
#' Aalen-Johansen estimator via Nelson-Aalen cumulative hazard increments and
#' product-integral construction.
#'
#' The Aalen-Johansen estimator generalizes the Kaplan-Meier estimator to
#' multistate models. RFmstate exposes it as a descriptive, covariate-free
#' calendar-time population benchmark from the recorded common baseline.
#'
#' @param msdata An \code{msdata} object from \code{\link{prepare_data}}.
#' @param s Numeric, the starting time for transition probabilities
#'   (default 0).
#'
#' @return An object of class \code{"aj_estimate"} containing:
#'   \describe{
#'     \item{time}{Numeric vector of unique event times.}
#'     \item{trans_prob}{List of calendar-time product-integral point-estimate
#'       matrices from the recorded common baseline at each event time.}
#'     \item{state_occ}{Matrix of state occupation probabilities over time.
#'       Rows are time points, columns are states.}
#'     \item{cum_hazard}{List of Nelson-Aalen cumulative hazard matrices.}
#'     \item{hazard_inc}{List of hazard increment matrices at each event time.}
#'     \item{n_risk}{Matrix of at-risk counts over time.}
#'     \item{n_events}{Data frame of event counts per transition.}
#'     \item{structure}{The multistate structure used.}
#'     \item{s}{The starting time.}
#'   }
#'
#' @examples
#' ms <- clinical_states()
#' set.seed(42)
#' dat <- sim_clinical_data(n = 200, structure = ms)
#' msdata <- prepare_data(
#'   data = dat, id = "ID", structure = ms,
#'   time_map = list(
#'     Responded = "time_Responded",
#'     Unresponded = "time_Unresponded",
#'     Stabilized = "time_Stabilized",
#'     Progressed = "time_Progressed",
#'     Death = "time_Death"
#'   ),
#'   censor_col = "time_censored",
#'   covariates = c("age", "sex", "BMI", "treatment")
#' )
#' aj <- aalen_johansen(msdata)
#' print(aj)
#'
#' @export
aalen_johansen <- function(msdata, s = 0) {
  if (!inherits(msdata, "msdata")) {
    stop("'msdata' must be an 'msdata' object from prepare_data().")
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s != 0) {
    stop("Only s = 0 is supported; left truncation and general conditional ",
         "P(s,t) are unavailable.")
  }

  structure <- attr(msdata, "structure")
  state_names <- structure$state_names
  ns <- structure$n_states
  trans_list <- structure$trans_list
  tmat <- structure$trans_matrix

  # Get all unique event times (where transitions occur)
  event_data <- msdata[msdata$status == 1 & msdata$Tstop > s, ]
  all_times <- sort(unique(event_data$Tstop))

  if (length(all_times) == 0) {
    stop("No events found after time s = ", s)
  }

  # Compute at-risk counts and event counts at each time
  # At-risk in state h at time t: patients in state h just before time t
  n_times <- length(all_times)

  # Initialize
  trans_prob_list <- vector("list", n_times)
  cum_haz_matrices <- vector("list", n_times)
  haz_inc_matrices <- vector("list", n_times)
  state_occ <- matrix(0, nrow = n_times, ncol = ns,
                       dimnames = list(NULL, state_names))
  n_risk_mat <- matrix(0, nrow = n_times, ncol = ns,
                       dimnames = list(NULL, state_names))

  # Product-integral: P(s, t_k) = prod_{j: t_j <= t_k} (I + dA(t_j))
  P_current <- diag(ns)
  dimnames(P_current) <- list(state_names, state_names)

  cum_haz <- matrix(0, nrow = ns, ncol = ns,
                    dimnames = list(state_names, state_names))

  for (k in seq_len(n_times)) {
    t_k <- all_times[k]

    # Count at-risk in each state just before t_k
    n_risk <- .count_at_risk(msdata, t_k, state_names)
    n_risk_mat[k, ] <- n_risk

    # Count transitions at t_k
    dN <- .count_transitions(msdata, t_k, state_names, tmat)

    # Compute hazard increments: dA_hj = dN_hj / n_risk_h
    dA <- matrix(0, nrow = ns, ncol = ns,
                 dimnames = list(state_names, state_names))

    for (h in seq_len(ns)) {
      if (n_risk[h] > 0) {
        for (j in seq_len(ns)) {
          if (h != j && !is.na(tmat[h, j])) {
            dA[h, j] <- dN[h, j] / n_risk[h]
          }
        }
        # Diagonal: dA_hh = -sum_j dA_hj
        dA[h, h] <- -sum(dA[h, -h], na.rm = TRUE)
      }
    }

    haz_inc_matrices[[k]] <- dA
    cum_haz <- cum_haz + dA
    cum_haz_matrices[[k]] <- cum_haz

    # Product-integral step: P = P * (I + dA)
    increment <- diag(ns) + dA
    dimnames(increment) <- list(state_names, state_names)
    P_current <- P_current %*% increment

    if (any(!is.finite(P_current)) || min(P_current) < -1e-12 ||
        max(P_current) > 1 + 1e-12 ||
        max(abs(rowSums(P_current) - 1)) > 1e-10) {
      stop("Aalen-Johansen probability invariant failed at time ", t_k, ".")
    }

    trans_prob_list[[k]] <- P_current

    # State occupation probabilities from the recorded common initial state.
    p0 <- rep(0, ns)
    initial_state <- attr(msdata, "initial_state")
    if (is.null(initial_state)) initial_state <- structure$initial_state
    p0[match(initial_state, state_names)] <- 1
    occ <- as.vector(p0 %*% P_current)
    state_occ[k, ] <- occ
  }

  # Event count summary
  event_summary <- .summarize_events(msdata, structure)

  result <- structure(
    list(
      time = all_times,
      trans_prob = trans_prob_list,
      state_occ = state_occ,
      cum_hazard = cum_haz_matrices,
      hazard_inc = haz_inc_matrices,
      n_risk = n_risk_mat,
      n_events = event_summary,
      structure = structure,
      initial_state = initial_state,
      s = 0,
      time_scale = "calendar time from common study origin",
      uncertainty = "No validated variance or confidence interval is returned."
    ),
    class = "aj_estimate"
  )
  result
}

#' Count patients at risk in each state just before time t
#' @noRd
.count_at_risk <- function(msdata, t, state_names) {
  # At risk in state h at time t: have Tstart < t and Tstop >= t,
  # and their 'from' state is h
  at_risk <- msdata[msdata$Tstart < t & msdata$Tstop >= t, ]
  n_risk <- vapply(state_names, function(s) {
    sum(at_risk$from == s)
  }, integer(1))
  n_risk
}

#' Count transitions at time t
#' @noRd
.count_transitions <- function(msdata, t, state_names, tmat) {
  ns <- length(state_names)
  dN <- matrix(0L, nrow = ns, ncol = ns,
               dimnames = list(state_names, state_names))

  events_at_t <- msdata[msdata$Tstop == t & msdata$status == 1, ]
  if (nrow(events_at_t) > 0) {
    for (r in seq_len(nrow(events_at_t))) {
      from <- events_at_t$from[r]
      to <- events_at_t$to[r]
      fi <- match(from, state_names)
      ti <- match(to, state_names)
      if (!is.na(fi) && !is.na(ti)) {
        dN[fi, ti] <- dN[fi, ti] + 1L
      }
    }
  }
  dN
}

#' Summarize event counts per transition
#' @noRd
.summarize_events <- function(msdata, structure) {
  events <- msdata[msdata$status == 1, ]
  trans_list <- structure$trans_list

  counts <- data.frame(
    trans_id = trans_list$trans_id,
    from = trans_list$from,
    to = trans_list$to,
    n_events = 0L,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(counts))) {
    counts$n_events[i] <- sum(
      events$from == counts$from[i] & events$to == counts$to[i],
      na.rm = TRUE
    )
  }

  # Add censoring counts per state
  censored <- msdata[msdata$status == 0, ]
  cens_tab <- table(censored$from)
  counts$n_censored_from <- 0L
  for (i in seq_len(nrow(counts))) {
    if (counts$from[i] %in% names(cens_tab)) {
      counts$n_censored_from[i] <- as.integer(cens_tab[counts$from[i]])
    }
  }

  counts
}

#' @export
print.aj_estimate <- function(x, ...) {
  cat("Aalen-Johansen Estimate\n")
  cat("  Time range: [", min(x$time), ", ", max(x$time), "]\n", sep = "")
  cat("  Event times:", length(x$time), "\n")
  cat("  States:", paste(x$structure$state_names, collapse = ", "), "\n")
  cat("  Common initial state:", x$initial_state, "\n")
  cat("  Uncertainty: point estimates only\n")
  cat("\nEvent counts per transition:\n")
  print(x$n_events[, c("from", "to", "n_events")], row.names = FALSE)
  cat("\nFinal state occupation probabilities:\n")
  final_occ <- x$state_occ[nrow(x$state_occ), ]
  for (i in seq_along(final_occ)) {
    cat("  ", names(final_occ)[i], ": ", round(final_occ[i], 4), "\n", sep = "")
  }
  invisible(x)
}
