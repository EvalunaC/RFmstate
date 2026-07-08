#' Define Multistate Structure
#'
#' Defines the state space, absorbing states, and allowed transitions for a
#' multistate model.
#'
#' @param state_names Character vector of state names.
#' @param absorbing Character vector of absorbing state names (must be a subset
#'   of \code{state_names}).
#' @param transitions A named list where each element name is an origin state
#'   and the value is a character vector of destination states reachable from
#'   that origin. Absorbing states should not appear as list names.
#'
#' @return An object of class \code{"mstate_structure"} containing:
#'   \describe{
#'     \item{state_names}{Character vector of all state names.}
#'     \item{n_states}{Integer, number of states.}
#'     \item{absorbing}{Character vector of absorbing states.}
#'     \item{transient}{Character vector of transient (non-absorbing) states.}
#'     \item{transitions}{Named list of allowed transitions.}
#'     \item{trans_matrix}{Integer matrix where entry \code{[i,j]} is the
#'       transition number for allowed transition i->j, or \code{NA} if not
#'       allowed.}
#'     \item{n_transitions}{Total number of allowed transitions.}
#'     \item{trans_list}{Data frame listing all transitions with columns
#'       \code{trans_id}, \code{from}, \code{to}.}
#'   }
#'
#' @examples
#' ms <- define_multistate(
#'   state_names = c("Baseline", "Responded", "Progressed", "Death"),
#'   absorbing = "Death",
#'   transitions = list(
#'     Baseline = c("Responded", "Progressed", "Death"),
#'     Responded = c("Progressed", "Death"),
#'     Progressed = c("Death")
#'   )
#' )
#' print(ms)
#'
#' @export
define_multistate <- function(state_names, absorbing, transitions) {
  # Validate inputs
  if (!is.character(state_names) || length(state_names) < 2) {
    stop("'state_names' must be a character vector with at least 2 states.")
  }
  if (any(duplicated(state_names))) {
    stop("'state_names' must not contain duplicates.")
  }
  if (!is.character(absorbing) || length(absorbing) < 1) {
    stop("'absorbing' must be a character vector with at least 1 absorbing state.")
  }
  if (!all(absorbing %in% state_names)) {
    stop("All absorbing states must be in 'state_names'.")
  }
  if (!is.list(transitions)) {
    stop("'transitions' must be a named list.")
  }

  transient <- setdiff(state_names, absorbing)

  # Validate transitions
  for (nm in names(transitions)) {
    if (!(nm %in% transient)) {
      stop("Transition origin '", nm,
           "' must be a transient (non-absorbing) state.")
    }
    dests <- transitions[[nm]]
    if (!all(dests %in% state_names)) {
      bad <- setdiff(dests, state_names)
      stop("Unknown destination state(s) in transitions from '", nm, "': ",
           paste(bad, collapse = ", "))
    }
    if (nm %in% dests) {
      stop("Self-transitions not allowed: state '", nm, "'.")
    }
  }

  # Check all transient states have outgoing transitions
  missing_origins <- setdiff(transient, names(transitions))
  if (length(missing_origins) > 0) {
    stop("Transient state(s) missing from 'transitions': ",
         paste(missing_origins, collapse = ", "))
  }

  # Build transition matrix
  ns <- length(state_names)
  tmat <- matrix(NA_integer_, nrow = ns, ncol = ns,
                 dimnames = list(state_names, state_names))

  trans_id <- 0L
  trans_list <- data.frame(
    trans_id = integer(0),
    from = character(0),
    to = character(0),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(state_names)) {
    from <- state_names[i]
    if (from %in% absorbing) next
    dests <- transitions[[from]]
    for (to in dests) {
      j <- match(to, state_names)
      trans_id <- trans_id + 1L
      tmat[i, j] <- trans_id
      trans_list <- rbind(trans_list, data.frame(
        trans_id = trans_id, from = from, to = to,
        stringsAsFactors = FALSE
      ))
    }
  }

  structure(
    list(
      state_names = state_names,
      n_states = ns,
      absorbing = absorbing,
      transient = transient,
      transitions = transitions,
      trans_matrix = tmat,
      n_transitions = trans_id,
      trans_list = trans_list
    ),
    class = "mstate_structure"
  )
}

#' @export
print.mstate_structure <- function(x, ...) {
  cat("Multistate Structure\n")
  cat("  States:", paste(x$state_names, collapse = " -> "), "\n")
  cat("  Absorbing:", paste(x$absorbing, collapse = ", "), "\n")
  cat("  Transitions:", x$n_transitions, "\n")
  for (i in seq_len(nrow(x$trans_list))) {
    cat("    ", x$trans_list$trans_id[i], ": ",
        x$trans_list$from[i], " -> ", x$trans_list$to[i], "\n", sep = "")
  }
  invisible(x)
}

#' Create Clinical Trial Multistate Structure
#'
#' A convenience function that creates the standard clinical trial multistate
#' structure with states: Baseline, Responded, Unresponded, Stabilized,
#' Progressed, Death.
#'
#' @return An \code{mstate_structure} object.
#'
#' @examples
#' ms <- clinical_states()
#' print(ms)
#'
#' @export
clinical_states <- function() {
  define_multistate(
    state_names = c("Baseline", "Responded", "Unresponded",
                    "Stabilized", "Progressed", "Death"),
    absorbing = "Death",
    transitions = list(
      Baseline = c("Responded", "Unresponded", "Death"),
      Responded = c("Stabilized", "Progressed", "Death"),
      Unresponded = c("Stabilized", "Progressed", "Death"),
      Stabilized = c("Progressed", "Death"),
      Progressed = c("Death")
    )
  )
}
