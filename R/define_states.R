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
#'     \item{initial_state}{The unique graph root used as the common initial
#'       state.}
#'     \item{topological_order}{State order used by the semi-Markov solver.}
#'   }
#'
#' @section Limitations:
#' Only directed acyclic, non-recurrent structures with exactly one root and at
#' least one absorbing state are supported. Every transient state must be
#' reachable from the root and must have a directed path to absorption. Cycles,
#' self-transitions, recurrent visits, disconnected states, and multiple
#' subject-specific initial states are rejected.
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
  if (!is.character(state_names) || length(state_names) < 2) {
    stop("'state_names' must be a character vector with at least 2 states.")
  }
  if (anyNA(state_names) || any(!nzchar(state_names)) ||
      anyDuplicated(state_names)) {
    stop("'state_names' must contain unique, nonmissing, nonempty names.")
  }
  if (!is.character(absorbing) || length(absorbing) < 1L ||
      anyNA(absorbing) || anyDuplicated(absorbing)) {
    stop("'absorbing' must be a character vector with at least 1 absorbing state.")
  }
  if (!all(absorbing %in% state_names)) {
    stop("All absorbing states must be in 'state_names'.")
  }
  if (!is.list(transitions) || is.null(names(transitions)) ||
      anyNA(names(transitions)) || any(!nzchar(names(transitions))) ||
      anyDuplicated(names(transitions))) {
    stop("'transitions' must be a named list.")
  }

  transient <- setdiff(state_names, absorbing)

  extra_origins <- setdiff(names(transitions), transient)
  if (length(extra_origins)) {
    stop("Transition origin(s) must be transient states: ",
         paste(extra_origins, collapse = ", "), ".")
  }
  missing_origins <- setdiff(transient, names(transitions))
  if (length(missing_origins)) {
    stop("Transient state(s) missing from 'transitions': ",
         paste(missing_origins, collapse = ", "))
  }

  for (nm in names(transitions)) {
    dests <- transitions[[nm]]
    if (!is.character(dests) || !length(dests) || anyNA(dests) ||
        any(!nzchar(dests))) {
      stop("Destinations from '", nm,
           "' must be a nonempty character vector without missing values.")
    }
    if (anyDuplicated(dests)) {
      stop("Duplicated directed edge(s) from '", nm, "' are not allowed.")
    }
    if (!all(dests %in% state_names)) {
      bad <- setdiff(dests, state_names)
      stop("Unknown destination state(s) in transitions from '", nm, "': ",
           paste(bad, collapse = ", "))
    }
    if (nm %in% dests) {
      stop("Self-transitions not allowed: state '", nm, "'.")
    }
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

  graph_info <- .validate_dag(
    state_names = state_names,
    absorbing = absorbing,
    trans_list = trans_list
  )

  structure(
    list(
      state_names = state_names,
      n_states = ns,
      absorbing = absorbing,
      transient = transient,
      transitions = transitions,
      trans_matrix = tmat,
      n_transitions = trans_id,
      trans_list = trans_list,
      initial_state = graph_info$initial_state,
      topological_order = graph_info$topological_order,
      graph_validation = graph_info
    ),
    class = "mstate_structure"
  )
}

#' Validate the supported single-root DAG contract
#' @noRd
.validate_dag <- function(state_names, absorbing, trans_list) {
  indegree <- stats::setNames(integer(length(state_names)), state_names)
  for (to in trans_list$to) indegree[to] <- indegree[to] + 1L

  roots <- state_names[indegree == 0L]
  if (length(roots) != 1L) {
    stop("The supported cohort graph must have exactly one common initial ",
         "state with no incoming edge; found ", length(roots), ": ",
         paste(roots, collapse = ", "), ".")
  }
  initial_state <- roots[[1L]]
  if (initial_state %in% absorbing) {
    stop("The common initial state must be transient.")
  }

  work_indegree <- indegree
  queue <- state_names[work_indegree == 0L]
  topo <- character(0)
  while (length(queue)) {
    node <- queue[[1L]]
    queue <- queue[-1L]
    topo <- c(topo, node)
    children <- trans_list$to[trans_list$from == node]
    for (child in children) {
      work_indegree[child] <- work_indegree[child] - 1L
      if (work_indegree[child] == 0L) {
        candidates <- c(queue, child)
        queue <- state_names[state_names %in% candidates]
      }
    }
  }
  if (length(topo) != length(state_names)) {
    cyclic <- state_names[work_indegree > 0L]
    stop("Directed cycles/recurrent state structures are unsupported; ",
         "cycle involves: ", paste(cyclic, collapse = ", "), ".")
  }

  reachable <- initial_state
  frontier <- initial_state
  while (length(frontier)) {
    children <- unique(trans_list$to[trans_list$from %in% frontier])
    new <- setdiff(children, reachable)
    reachable <- c(reachable, new)
    frontier <- new
  }
  unreachable <- setdiff(state_names, reachable)
  if (length(unreachable)) {
    stop("State(s) unreachable from common initial state '", initial_state,
         "': ", paste(unreachable, collapse = ", "), ".")
  }

  can_absorb <- absorbing
  repeat {
    parents <- unique(trans_list$from[trans_list$to %in% can_absorb])
    enlarged <- union(can_absorb, parents)
    if (length(enlarged) == length(can_absorb)) break
    can_absorb <- enlarged
  }
  no_absorbing_path <- setdiff(setdiff(state_names, absorbing), can_absorb)
  if (length(no_absorbing_path)) {
    stop("Transient state(s) without a directed path to absorption: ",
         paste(no_absorbing_path, collapse = ", "), ".")
  }

  list(
    valid = TRUE,
    initial_state = initial_state,
    topological_order = topo,
    indegree = indegree,
    reachable = reachable
  )
}

#' Print a Multistate Structure
#'
#' Prints the user-defined display order, absorbing states, common initial
#' state, computational topological order, and numbered allowed transitions.
#'
#' @param x An \code{mstate_structure} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @section Limitations:
#' The printout summarizes the validated graph; it does not imply support for
#' cycles, recurrent visits, left truncation, or time-dependent covariates.
#'
#' @examples
#' ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
#' print(ms)
#'
#' @export
print.mstate_structure <- function(x, ...) {
  cat("Multistate Structure\n")
  cat("  States:", paste(x$state_names, collapse = " -> "), "\n")
  cat("  Absorbing:", paste(x$absorbing, collapse = ", "), "\n")
  cat("  Common initial state:", x$initial_state, "\n")
  cat("  Computational order:",
      paste(x$topological_order, collapse = " -> "), "\n")
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
#' @section Limitations:
#' This is a fixed demonstration DAG, not a claim that the package supports
#' cyclic or recurrent clinical histories. Use \code{define_multistate()} for
#' another supported single-root DAG.
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
