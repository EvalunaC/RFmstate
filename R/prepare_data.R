#' Prepare Data for Multistate Analysis
#'
#' Converts one-row-per-subject wide data into validated, long-form sojourn
#' records for an acyclic, non-recurrent multistate process. Entry times are
#' calendar times from a common origin; forest response times are the resulting
#' clock-reset durations within states.
#'
#' @param data A data frame with one row per subject.
#' @param id Name of the unique subject-ID column.
#' @param structure An \code{mstate_structure} from
#'   \code{\link{define_multistate}}.
#' @param time_map Named list mapping every noninitial state exactly once to its
#'   first and only calendar-time entry column. The initial state is excluded.
#' @param censor_col Name of the external right-censoring column.
#' @param covariates Character vector of baseline, time-fixed covariates.
#' @param initial_state Common initial state. It must equal the unique graph
#'   root; the default is stored in \code{structure}.
#'
#' @return An \code{msdata} data frame containing \code{id}, \code{from},
#'   \code{to}, \code{Tstart}, \code{Tstop}, \code{status},
#'   \code{trans_id}, \code{duration}, and retained covariates. Metadata store
#'   the graph, input mapping, initial state, predictor schema, counts, and
#'   observed duration support.
#'
#' @details Delayed entry, recurrent visits, tied entry times, time-dependent
#'   covariates, and missing covariates are unsupported. Invalid trajectories
#'   produce an error; no event or interval is silently discarded.
#'
#' @examples
#' ms <- clinical_states()
#' dat <- sim_clinical_data(n = 50, structure = ms, seed = 42)
#' msdata <- prepare_data(
#'   dat, id = "ID", structure = ms,
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
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!inherits(structure, "mstate_structure")) {
    stop("'structure' must be an 'mstate_structure' object.")
  }
  .scalar_column_name(id, "id")
  .scalar_column_name(censor_col, "censor_col")
  if (!(id %in% names(data))) stop("Column '", id, "' not found in data.")
  if (!(censor_col %in% names(data))) {
    stop("Column '", censor_col, "' not found in data.")
  }

  if (is.null(initial_state)) initial_state <- structure$initial_state
  if (!is.character(initial_state) || length(initial_state) != 1L ||
      is.na(initial_state) || !(initial_state %in% structure$state_names)) {
    stop("'initial_state' must be one valid state name.")
  }
  if (!identical(initial_state, structure$initial_state)) {
    stop("'initial_state' must equal the graph's unique common initial state '",
         structure$initial_state, "'.")
  }

  if (!is.list(time_map) || is.null(names(time_map)) ||
      anyNA(names(time_map)) || any(!nzchar(names(time_map))) ||
      anyDuplicated(names(time_map))) {
    stop("'time_map' must be a named list with unique state names.")
  }
  expected_states <- setdiff(structure$state_names, initial_state)
  missing_states <- setdiff(expected_states, names(time_map))
  extra_states <- setdiff(names(time_map), expected_states)
  if (length(missing_states) || length(extra_states)) {
    stop("'time_map' must map every noninitial state exactly once. Missing: ",
         if (length(missing_states)) paste(missing_states, collapse = ", ") else "none",
         "; unexpected: ",
         if (length(extra_states)) paste(extra_states, collapse = ", ") else "none",
         ".")
  }
  map_cols <- unlist(time_map, use.names = FALSE)
  if (!is.character(map_cols) || anyNA(map_cols) || any(!nzchar(map_cols)) ||
      anyDuplicated(map_cols)) {
    stop("Every 'time_map' value must be a unique, nonmissing column name.")
  }
  missing_time_cols <- setdiff(map_cols, names(data))
  if (length(missing_time_cols)) {
    stop("Time column(s) not found in data: ",
         paste(missing_time_cols, collapse = ", "), ".")
  }

  if (missing(covariates) || is.null(covariates)) covariates <- character(0)
  if (!is.character(covariates) || anyNA(covariates) ||
      anyDuplicated(covariates)) {
    stop("'covariates' must be a unique character vector.")
  }
  missing_covs <- setdiff(covariates, names(data))
  if (length(missing_covs)) {
    stop("Covariate column(s) not found: ", paste(missing_covs, collapse = ", "), ".")
  }

  ids <- data[[id]]
  if (anyNA(ids) || (is.character(ids) && any(!nzchar(ids)))) {
    stop("Subject IDs must be nonmissing and nonempty.")
  }
  if (anyDuplicated(ids)) stop("Subject IDs must be unique.")
  coerced_ids <- as.character(ids)
  if (anyDuplicated(coerced_ids)) {
    stop("Subject IDs collide after internal character representation.")
  }

  work_data <- data
  schema <- vector("list", length(covariates))
  names(schema) <- covariates
  for (nm in covariates) {
    x <- work_data[[nm]]
    if (anyNA(x)) stop("Covariate '", nm, "' contains missing values.")
    if (is.character(x)) {
      x <- factor(x, levels = unique(x))
      work_data[[nm]] <- x
    }
    if (is.numeric(x) && any(!is.finite(x))) {
      stop("Numeric covariate '", nm, "' contains nonfinite values.")
    }
    if (!(is.numeric(x) || is.integer(x) || is.logical(x) || is.factor(x))) {
      stop("Unsupported covariate class for '", nm, "': ",
           paste(class(x), collapse = "/"), ".")
    }
    schema[[nm]] <- list(
      class = class(x),
      levels = if (is.factor(x)) levels(x) else NULL,
      ordered = is.ordered(x),
      range = if (is.numeric(x)) range(x) else NULL
    )
  }

  all_rows <- vector("list", nrow(work_data))
  for (i in seq_len(nrow(work_data))) {
    patient <- work_data[i, , drop = FALSE]
    pid <- patient[[id]][[1L]]
    label <- as.character(pid)
    censor_time <- patient[[censor_col]][[1L]]
    if (!is.na(censor_time) &&
        (!is.numeric(censor_time) || length(censor_time) != 1L ||
         !is.finite(censor_time) || censor_time < 0)) {
      stop("Subject '", label, "' has an invalid censoring time.")
    }

    values <- vapply(expected_states, function(state) {
      value <- patient[[time_map[[state]]]][[1L]]
      if (is.na(value)) return(NA_real_)
      if (!is.numeric(value) || !is.finite(value) || value < 0) {
        stop("Subject '", label, "' has an invalid entry time for state '",
             state, "'.")
      }
      as.numeric(value)
    }, numeric(1))
    visited <- !is.na(values)
    events <- data.frame(
      state = expected_states[visited],
      time = values[visited],
      stringsAsFactors = FALSE
    )
    if (nrow(events)) {
      events <- events[order(events$time), , drop = FALSE]
      if (anyDuplicated(events$time)) {
        tied <- paste(events$state[duplicated(events$time) |
                                     duplicated(events$time, fromLast = TRUE)],
                      collapse = ", ")
        stop("Subject '", label, "' has simultaneous state-entry times: ",
             tied, ".")
      }
      if (!is.na(censor_time) && any(events$time > censor_time)) {
        bad <- events$state[which(events$time > censor_time)[1L]]
        stop("Subject '", label, "' has event entry into '", bad,
             "' after censoring.")
      }
    }

    trajectory <- .validated_trajectory(
      events = events,
      initial_state = initial_state,
      censor_time = censor_time,
      structure = structure,
      subject_label = label
    )
    covs <- patient[, covariates, drop = FALSE]
    all_rows[[i]] <- .trajectory_to_rows(pid, trajectory, covs, structure)
  }

  result <- do.call(rbind, all_rows)
  rownames(result) <- NULL
  result$from <- factor(result$from, levels = structure$state_names)
  result$to <- factor(result$to, levels = structure$state_names)
  result$Tstart <- as.numeric(result$Tstart)
  result$Tstop <- as.numeric(result$Tstop)
  result$status <- as.integer(result$status)
  result$trans_id <- as.integer(result$trans_id)
  result$duration <- result$Tstop - result$Tstart
  if (any(!is.finite(result$duration)) || any(result$duration <= 0)) {
    bad <- which(!is.finite(result$duration) | result$duration <= 0)[1L]
    stop("Internal validation failed: nonpositive duration for subject '",
         as.character(result$id[bad]), "' in state '",
         as.character(result$from[bad]), "'.")
  }

  edge_counts <- structure$trans_list
  edge_counts$n_events <- vapply(seq_len(nrow(edge_counts)), function(k) {
    sum(result$status == 1L &
          as.character(result$from) == edge_counts$from[k] &
          as.character(result$to) == edge_counts$to[k], na.rm = TRUE)
  }, integer(1))
  origin_support <- stats::setNames(vapply(structure$transient, function(state) {
    vals <- result$duration[as.character(result$from) == state]
    if (length(vals)) max(vals) else NA_real_
  }, numeric(1)), structure$transient)
  origin_sojourns <- stats::setNames(vapply(structure$transient, function(state) {
    sum(as.character(result$from) == state)
  }, integer(1)), structure$transient)
  origin_person_time <- stats::setNames(vapply(structure$transient, function(state) {
    sum(result$duration[as.character(result$from) == state])
  }, numeric(1)), structure$transient)

  metadata <- list(
    initial_state = initial_state,
    original_id_column = id,
    original_id_class = class(data[[id]]),
    time_map = time_map,
    censor_col = censor_col,
    covariates = covariates,
    covariate_schema = schema,
    n_subjects = nrow(data),
    n_intervals = nrow(result),
    n_events = sum(result$status == 1L),
    n_censored = sum(result$status == 0L),
    edge_counts = edge_counts,
    origin_sojourns = origin_sojourns,
    origin_person_time = origin_person_time,
    max_duration_by_origin = origin_support,
    graph_validation = structure$graph_validation,
    display_order = structure$state_names,
    topological_order = structure$topological_order
  )
  attr(result, "structure") <- structure
  attr(result, "initial_state") <- initial_state
  attr(result, "covariates") <- covariates
  attr(result, "covariate_schema") <- schema
  attr(result, "metadata") <- metadata
  class(result) <- c("msdata", "data.frame")
  result
}

#' @noRd
.scalar_column_name <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("'", label, "' must name exactly one column.")
  }
}

#' Construct one validated subject trajectory
#' @noRd
.validated_trajectory <- function(events, initial_state, censor_time,
                                  structure, subject_label) {
  trajectory <- list()
  current_state <- initial_state
  current_time <- 0

  if (nrow(events)) {
    for (j in seq_len(nrow(events))) {
      next_state <- events$state[j]
      event_time <- events$time[j]
      if (current_state %in% structure$absorbing) {
        stop("Subject '", subject_label, "' has event entry into '", next_state,
             "' after absorption in '", current_state, "'.")
      }
      if (event_time <= current_time) {
        stop("Subject '", subject_label,
             "' has nonincreasing entry times at transition '", current_state,
             " -> ", next_state, "'.")
      }
      allowed <- structure$transitions[[current_state]]
      if (!(next_state %in% allowed)) {
        stop("Subject '", subject_label, "' has forbidden observed transition '",
             current_state, " -> ", next_state, "'.")
      }
      trajectory[[length(trajectory) + 1L]] <- list(
        from = current_state, to = next_state,
        Tstart = current_time, Tstop = event_time, status = 1L
      )
      current_state <- next_state
      current_time <- event_time
    }
  }

  if (current_state %in% structure$absorbing) {
    if (!is.na(censor_time) && censor_time < current_time) {
      stop("Subject '", subject_label,
           "' has censoring before the observed absorbing-state entry.")
    }
    return(trajectory)
  }

  if (is.na(censor_time)) {
    stop("Subject '", subject_label,
         "' does not reach absorption and requires a finite censoring time.")
  }
  if (censor_time <= current_time) {
    stop("Subject '", subject_label,
         "' has censoring time that is not strictly after the last state entry.")
  }
  trajectory[[length(trajectory) + 1L]] <- list(
    from = current_state, to = NA_character_,
    Tstart = current_time, Tstop = censor_time, status = 0L
  )
  trajectory
}

#' Convert a trajectory to long-form rows
#' @noRd
.trajectory_to_rows <- function(pid, trajectory, covs, structure) {
  rows <- lapply(trajectory, function(tr) {
    trans_id <- if (tr$status == 1L) {
      structure$trans_matrix[tr$from, tr$to]
    } else {
      NA_integer_
    }
    row <- data.frame(
      id = pid,
      from = tr$from,
      to = if (tr$status == 1L) tr$to else NA_character_,
      Tstart = tr$Tstart,
      Tstop = tr$Tstop,
      status = tr$status,
      trans_id = trans_id,
      stringsAsFactors = FALSE
    )
    cbind(row, covs)
  })
  do.call(rbind, rows)
}

#' @export
print.msdata <- function(x, ...) {
  structure <- attr(x, "structure")
  cat("Multistate Data (msdata)\n")
  cat("  Patients:", length(unique(x$id)), "\n")
  cat("  Intervals:", nrow(x), "\n")
  cat("  Transitions observed:", sum(x$status == 1L), "\n")
  cat("  Censored intervals:", sum(x$status == 0L), "\n")
  if (!is.null(structure)) {
    cat("  Initial state:", attr(x, "initial_state"), "\n")
    cat("  States:", paste(structure$state_names, collapse = ", "), "\n")
  }
  cat("\nTransition counts:\n")
  events <- x[x$status == 1L, , drop = FALSE]
  if (nrow(events)) {
    print(table(
      from = factor(events$from, levels = structure$state_names),
      to = factor(events$to, levels = structure$state_names)
    ))
  } else {
    cat("  No observed transitions.\n")
  }
  invisible(x)
}

#' Return ordinary data rows from an msdata object
#'
#' @param x An \code{msdata} object.
#' @param n Number of rows to return.
#' @param ... Ignored.
#' @return An ordinary \code{data.frame} containing the first \code{n} rows.
#' @export
head.msdata <- function(x, n = 6L, ...) {
  n <- max(0L, min(as.integer(n), nrow(x)))
  out <- x[seq_len(n), , drop = FALSE]
  class(out) <- "data.frame"
  out
}
