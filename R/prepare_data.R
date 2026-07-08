#' Prepare Data for Multistate Analysis
#'
#' Converts wide-format clinical data into long counting-process format
#' suitable for multistate survival analysis.
#'
#' Each patient's trajectory is reconstructed from event times, validated
#' against the allowed transitions, and expanded into start-stop intervals with
#' covariates.
#'
#' @param data A data frame in wide format with one row per patient.
#' @param id Character string, name of the patient ID column.
#' @param structure An \code{mstate_structure} object from
#'   \code{\link{define_multistate}}.
#' @param time_map A named list mapping state names to column names in
#'   \code{data} containing the time of entry into that state (measured from
#'   baseline). The initial state should not be included. Use \code{NA} in the
#'   data for states not visited by a patient.
#' @param censor_col Character string, name of the column containing the right
#'   censoring time (last follow-up time).
#' @param covariates Character vector of covariate column names to carry into
#'   the long-format data.
#' @param initial_state Character string, the starting state for all patients
#'   (default: first state in the structure).
#'
#' @return An object of class \code{"msdata"} (a data frame) with columns:
#'   \describe{
#'     \item{id}{Patient identifier.}
#'     \item{from}{Origin state for this interval.}
#'     \item{to}{Destination state (or \code{NA} if censored).}
#'     \item{Tstart}{Start time of the interval.}
#'     \item{Tstop}{End time of the interval.}
#'     \item{status}{1 if a transition occurred, 0 if censored.}
#'     \item{trans_id}{Integer transition ID (from structure) or \code{NA}.}
#'     \item{duration}{Duration of the interval.}
#'     \item{...}{Covariate columns.}
#'   }
#'   The object also carries an attribute \code{"structure"} (the
#'   \code{mstate_structure}).
#'
#' @examples
#' ms <- clinical_states()
#' set.seed(42)
#' dat <- sim_clinical_data(n = 50, structure = ms)
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
#' head(msdata)
#'
#' @export
prepare_data <- function(data, id, structure, time_map, censor_col,
                         covariates, initial_state = NULL) {
  if (!inherits(structure, "mstate_structure")) {
    stop("'structure' must be an 'mstate_structure' object.")
  }
  if (!(id %in% names(data))) {
    stop("Column '", id, "' not found in data.")
  }
  if (!(censor_col %in% names(data))) {
    stop("Column '", censor_col, "' not found in data.")
  }
  for (col in covariates) {
    if (!(col %in% names(data))) {
      stop("Covariate column '", col, "' not found in data.")
    }
  }
  for (nm in names(time_map)) {
    tcol <- time_map[[nm]]
    if (!(tcol %in% names(data))) {
      stop("Time column '", tcol, "' for state '", nm, "' not found in data.")
    }
  }

  if (is.null(initial_state)) {
    initial_state <- structure$state_names[1]
  }

  # Build per-patient trajectories
  all_rows <- vector("list", nrow(data))

  for (i in seq_len(nrow(data))) {
    patient <- data[i, ]
    pid <- patient[[id]]
    cens_time <- patient[[censor_col]]

    # Collect event times for this patient
    events <- data.frame(
      state = character(0), time = numeric(0),
      stringsAsFactors = FALSE
    )
    for (nm in names(time_map)) {
      t_val <- patient[[time_map[[nm]]]]
      if (!is.na(t_val) && is.finite(t_val) && t_val >= 0) {
        events <- rbind(events, data.frame(
          state = nm, time = t_val, stringsAsFactors = FALSE
        ))
      }
    }

    # Sort events by time
    if (nrow(events) > 0) {
      events <- events[order(events$time), , drop = FALSE]
    }

    # Build trajectory following allowed transitions
    trajectory <- .build_trajectory(
      events, initial_state, cens_time, structure
    )

    # Convert trajectory to counting-process rows
    covs <- patient[covariates]
    rows <- .trajectory_to_rows(pid, trajectory, covs, structure)
    all_rows[[i]] <- rows
  }

  result <- do.call(rbind, all_rows)
  rownames(result) <- NULL

  # Ensure correct types
  result$id <- as.character(result$id)
  result$from <- as.character(result$from)
  result$to <- as.character(result$to)
  result$Tstart <- as.numeric(result$Tstart)
  result$Tstop <- as.numeric(result$Tstop)
  result$status <- as.integer(result$status)
  result$trans_id <- as.integer(result$trans_id)
  result$duration <- result$Tstop - result$Tstart

  attr(result, "structure") <- structure
  class(result) <- c("msdata", "data.frame")
  result
}

#' Build patient trajectory from sorted events
#' @noRd
.build_trajectory <- function(events, initial_state, cens_time, structure) {
  trajectory <- list()
  current_state <- initial_state
  current_time <- 0

  if (nrow(events) == 0) {
    # Patient was censored from initial state
    trajectory[[1]] <- list(
      from = current_state, to = NA_character_,
      Tstart = current_time, Tstop = cens_time, status = 0L
    )
    return(trajectory)
  }

  for (j in seq_len(nrow(events))) {
    next_state <- events$state[j]
    event_time <- events$time[j]

    # Check if this transition is allowed
    allowed <- structure$transitions[[current_state]]
    if (is.null(allowed) || !(next_state %in% allowed)) {
      next
    }

    # Check if censoring happens before this event
    if (!is.na(cens_time) && cens_time < event_time) {
      trajectory[[length(trajectory) + 1]] <- list(
        from = current_state, to = NA_character_,
        Tstart = current_time, Tstop = cens_time, status = 0L
      )
      current_state <- NA_character_
      break
    }

    # Valid transition
    trajectory[[length(trajectory) + 1]] <- list(
      from = current_state, to = next_state,
      Tstart = current_time, Tstop = event_time, status = 1L
    )
    current_time <- event_time
    current_state <- next_state

    # If we reached an absorbing state, stop
    if (current_state %in% structure$absorbing) break
  }

  # If patient is still in a transient state and not censored yet
  if (!is.na(current_state) && !(current_state %in% structure$absorbing)) {
    if (!is.na(cens_time) && cens_time > current_time) {
      trajectory[[length(trajectory) + 1]] <- list(
        from = current_state, to = NA_character_,
        Tstart = current_time, Tstop = cens_time, status = 0L
      )
    }
  }

  trajectory
}

#' Convert trajectory to data frame rows
#' @noRd
.trajectory_to_rows <- function(pid, trajectory, covs, structure) {
  if (length(trajectory) == 0) return(NULL)

  rows <- vector("list", length(trajectory))
  for (k in seq_along(trajectory)) {
    tr <- trajectory[[k]]
    # Get transition ID
    if (tr$status == 1L) {
      from_idx <- match(tr$from, structure$state_names)
      to_idx <- match(tr$to, structure$state_names)
      trans_id <- structure$trans_matrix[from_idx, to_idx]
    } else {
      trans_id <- NA_integer_
    }

    row <- data.frame(
      id = pid,
      from = tr$from,
      to = ifelse(tr$status == 1L, tr$to, NA_character_),
      Tstart = tr$Tstart,
      Tstop = tr$Tstop,
      status = tr$status,
      trans_id = trans_id,
      stringsAsFactors = FALSE
    )
    row <- cbind(row, covs)
    rows[[k]] <- row
  }
  do.call(rbind, rows)
}

#' @export
print.msdata <- function(x, ...) {
  structure <- attr(x, "structure")
  cat("Multistate Data (msdata)\n")
  cat("  Patients:", length(unique(x$id)), "\n")
  cat("  Intervals:", nrow(x), "\n")
  cat("  Transitions observed:", sum(x$status == 1), "\n")
  cat("  Censored intervals:", sum(x$status == 0), "\n")
  if (!is.null(structure)) {
    cat("  States:", paste(structure$state_names, collapse = ", "), "\n")
  }
  cat("\nTransition counts:\n")
  trans_events <- x[x$status == 1, ]
  if (nrow(trans_events) > 0) {
    tab <- table(
      from = trans_events$from,
      to = trans_events$to
    )
    print(tab)
  }
  invisible(x)
}
