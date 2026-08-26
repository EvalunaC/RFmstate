#' Fit Clock-Reset Random-Forest Multistate Models
#'
#' Fits one ranger cause-specific survival forest for every declared edge of a
#' validated acyclic, non-recurrent multistate structure. Forest response time
#' is duration since fresh entry into the origin state.
#'
#' @param msdata An \code{msdata} object from \code{\link{prepare_data}}.
#' @param covariates Predictor names. \code{NULL} uses the stored preparation
#'   schema.
#' @param num.trees Number of trees per edge forest.
#' @param mtry Number of predictors considered at a split.
#' @param min.node.size Ranger minimum node size.
#' @param min_events Explicit technical safeguard for target events per edge.
#'   It is not a universal adequacy threshold.
#' @param sparse_warning Warn for edges below this event count; set
#'   \code{Inf} to disable.
#' @param importance Ranger importance mode.
#' @param seed Reproducibility seed.
#' @param ... Supported ranger controls: \code{replace},
#'   \code{sample.fraction}, \code{splitrule}, \code{num.random.splits},
#'   \code{respect.unordered.factors}, \code{num.threads},
#'   \code{save.memory}, \code{max.depth}, \code{always.split.variables},
#'   \code{alpha}, \code{minprop}, and \code{verbose}. Unnamed, conflicting,
#'   and unknown arguments are rejected.
#'
#' @return An \code{rfmstate} object containing all edge models, schemas,
#'   event/support metadata, OOB information, and the approved clock-reset
#'   time-scale contract.
#'
#' @details Every competing exit from an origin state remains an observed exit
#'   for that sojourn, but is coded as a non-target outcome in the binary
#'   cause-specific forest for a particular edge. An unestimable declared edge
#'   stops the whole fit; it is never omitted or represented by zero hazard.
#'
#' @examples
#' \donttest{
#' ms <- clinical_states()
#' dat <- sim_clinical_data(300, structure = ms, seed = 42)
#' long <- prepare_data(
#'   dat, "ID", ms,
#'   list(Responded="time_Responded", Unresponded="time_Unresponded",
#'        Stabilized="time_Stabilized", Progressed="time_Progressed",
#'        Death="time_Death"),
#'   "time_censored", c("age", "sex", "BMI", "treatment")
#' )
#' fit <- rfmstate(long, num.trees = 100, min_events = 3, seed = 42)
#' }
#'
#' @export
rfmstate <- function(msdata, covariates = NULL, num.trees = 1000L,
                     mtry = NULL, min.node.size = 15L, min_events = 5L,
                     sparse_warning = 50L,
                     importance = "permutation", seed = NULL, ...) {
  if (!inherits(msdata, "msdata")) {
    stop("'msdata' must be an 'msdata' object from prepare_data().")
  }
  structure_obj <- attr(msdata, "structure")
  metadata <- attr(msdata, "metadata")
  if (!inherits(structure_obj, "mstate_structure") || is.null(metadata)) {
    stop("'msdata' lacks validated structure/metadata; rerun prepare_data().")
  }
  cl <- match.call()
  covariates <- covariates %||% attr(msdata, "covariates")
  if (!is.character(covariates) || !length(covariates) || anyNA(covariates) ||
      anyDuplicated(covariates)) {
    stop("At least one unique, nonmissing fitting covariate is required.")
  }
  missing_covs <- setdiff(covariates, names(msdata))
  if (length(missing_covs)) {
    stop("Covariates not found in msdata: ", paste(missing_covs, collapse = ", "))
  }
  for (nm in covariates) {
    if (anyNA(msdata[[nm]])) stop("Covariate '", nm, "' contains missing values.")
    if (is.numeric(msdata[[nm]]) && any(!is.finite(msdata[[nm]]))) {
      stop("Numeric covariate '", nm, "' contains nonfinite values.")
    }
  }

  .positive_integer(num.trees, "num.trees")
  .positive_integer(min.node.size, "min.node.size")
  .positive_integer(min_events, "min_events")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed))) {
    stop("'seed' must be NULL or one finite number.")
  }
  if (!(importance %in% c("permutation", "impurity", "impurity_corrected",
                          "none"))) {
    stop("Unsupported ranger importance mode.")
  }
  if (is.null(mtry)) mtry <- max(1L, floor(sqrt(length(covariates))))
  .positive_integer(mtry, "mtry")
  mtry <- min(as.integer(mtry), length(covariates))

  dots <- list(...)
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || anyNA(dot_names) ||
                       any(!nzchar(dot_names)) || anyDuplicated(dot_names))) {
    stop("All ranger arguments in '...' must be uniquely named.")
  }
  whitelist <- c(
    "replace", "sample.fraction", "splitrule", "num.random.splits",
    "respect.unordered.factors", "num.threads", "save.memory",
    "max.depth", "always.split.variables", "alpha", "minprop", "verbose"
  )
  unknown <- setdiff(dot_names, whitelist)
  if (length(unknown)) {
    stop("Unsupported or RFmstate-controlled ranger argument(s): ",
         paste(unknown, collapse = ", "), ".")
  }
  ranger_formals <- names(formals(ranger::ranger))
  unavailable <- setdiff(dot_names, ranger_formals)
  if (length(unavailable)) {
    stop("Forwarded ranger argument(s) unavailable in ranger ",
         as.character(utils::packageVersion("ranger")), ": ",
         paste(unavailable, collapse = ", "), ".")
  }

  origin_data <- stats::setNames(vector("list", length(structure_obj$transient)),
                                 structure_obj$transient)
  failures <- character(0)
  sparse <- character(0)
  edge_precheck <- vector("list", nrow(structure_obj$trans_list))
  names(edge_precheck) <- paste0(structure_obj$trans_list$from, "->",
                                 structure_obj$trans_list$to)
  for (state in structure_obj$transient) {
    rows <- msdata[as.character(msdata$from) == state, , drop = FALSE]
    cs_data <- .build_cause_specific_data(
      rows, structure_obj$transitions[[state]], covariates
    )
    origin_data[[state]] <- cs_data
    for (dest in structure_obj$transitions[[state]]) {
      edge <- paste0(state, "->", dest)
      dest_idx <- match(dest, structure_obj$transitions[[state]])
      target <- cs_data$event_type == dest_idx
      n_events <- sum(target)
      n_distinct <- length(unique(cs_data$duration[target]))
      details <- list(
        edge = edge,
        n_sojourns = nrow(cs_data),
        n_events = n_events,
        n_distinct_event_times = n_distinct,
        n_competing_exits = sum(cs_data$event_type > 0L & !target),
        n_external_censored = sum(cs_data$event_type == 0L)
      )
      edge_precheck[[edge]] <- details
      criteria <- character(0)
      if (!nrow(cs_data)) criteria <- c(criteria, "no origin-state sojourns")
      if (n_events < min_events) {
        criteria <- c(criteria, paste0("target events ", n_events,
                                      " < min_events ", min_events))
      }
      if (nrow(cs_data) <= min.node.size) {
        criteria <- c(criteria, paste0("origin sojourns ", nrow(cs_data),
                                      " <= min.node.size ", min.node.size))
      }
      if (length(criteria)) {
        failures <- c(
          failures,
          paste0(edge, " [sojourns=", nrow(cs_data), ", target_events=",
                 n_events, ", distinct_event_times=", n_distinct,
                 "]: ", paste(criteria, collapse = "; "))
        )
      } else if (is.finite(sparse_warning) && n_events < sparse_warning) {
        sparse <- c(sparse, paste0(edge, " (", n_events, " events)"))
      }
    }
  }
  if (length(failures)) {
    stop("Every declared edge must be estimable. Failed edge(s):\n- ",
         paste(failures, collapse = "\n- "),
         "\nRevise the declared graph, combine clinically equivalent states, ",
         "collect more data, or change the documented min_events safeguard.")
  }
  if (length(sparse)) {
    warning("Sparse transition(s); interpret edge diagnostics cautiously: ",
            paste(sparse, collapse = ", "), ". The threshold is descriptive, ",
            "not a universal adequacy rule.", call. = FALSE)
  }

  models <- list()
  event_times <- list()
  edge_metadata <- list()
  exact_args <- c(list(
    num.trees = as.integer(num.trees), mtry = mtry,
    min.node.size = as.integer(min.node.size), importance = importance,
    seed = seed, write.forest = TRUE, oob.error = TRUE
  ), dots)
  for (state in structure_obj$transient) {
    cs_data <- origin_data[[state]]
    models[[state]] <- list()
    event_times[[state]] <- numeric(0)
    for (dest in structure_obj$transitions[[state]]) {
      edge <- paste0(state, "->", dest)
      dest_idx <- match(dest, structure_obj$transitions[[state]])
      status <- as.integer(cs_data$event_type == dest_idx)
      fit_data <- data.frame(
        time = cs_data$duration,
        status = status,
        cs_data[, covariates, drop = FALSE],
        check.names = FALSE
      )
      formula <- stats::as.formula(
        paste("survival::Surv(time, status) ~",
              paste(sprintf("`%s`", covariates), collapse = " + "))
      )
      call_args <- c(list(
        formula = formula,
        data = fit_data,
        num.trees = as.integer(num.trees),
        mtry = mtry,
        min.node.size = as.integer(min.node.size),
        importance = importance,
        seed = seed,
        write.forest = TRUE,
        oob.error = TRUE
      ), dots)
      rf_fit <- do.call(ranger::ranger, call_args)
      models[[state]][[dest]] <- rf_fit
      event_times[[state]] <- sort(unique(c(
        event_times[[state]], rf_fit$unique.death.times
      )))
      info <- edge_precheck[[edge]]
      edge_metadata[[edge]] <- c(info, list(
        backend_event_times = rf_fit$unique.death.times,
        prediction_error = rf_fit$prediction.error,
        oob_concordance = 1 - rf_fit$prediction.error,
        ranger_arguments = exact_args
      ))
    }
  }

  expected_edges <- paste0(structure_obj$trans_list$from, "->",
                           structure_obj$trans_list$to)
  fitted_edges <- unlist(lapply(names(models), function(from) {
    paste0(from, "->", names(models[[from]]))
  }), use.names = FALSE)
  missing_models <- setdiff(expected_edges, fitted_edges)
  if (length(missing_models)) {
    stop("Internal fit failure: missing declared edge model(s): ",
         paste(missing_models, collapse = ", "), ".")
  }

  schema <- attr(msdata, "covariate_schema")[covariates]
  support <- metadata$max_duration_by_origin[structure_obj$transient]
  structure(
    list(
      models = models,
      structure = structure_obj,
      covariates = covariates,
      predictor_schema = schema,
      factor_levels = lapply(schema, `[[`, "levels"),
      origin_data = origin_data,
      event_times = event_times,
      edge_metadata = edge_metadata,
      event_counts = metadata$edge_counts,
      max_duration_by_origin = support,
      initial_state = attr(msdata, "initial_state"),
      time_scale = "clock-reset",
      process_model = "semi-Markov",
      msdata = msdata,
      call = cl,
      params = list(
        num.trees = as.integer(num.trees), mtry = mtry,
        min.node.size = as.integer(min.node.size),
        min_events = as.integer(min_events), sparse_warning = sparse_warning,
        importance = importance, seed = seed,
        ranger_args = dots
      ),
      package_versions = list(
        RFmstate = as.character(utils::packageVersion("RFmstate")),
        ranger = as.character(utils::packageVersion("ranger")),
        R = as.character(getRversion())
      )
    ),
    class = "rfmstate"
  )
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' @noRd
.positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
      x != as.integer(x)) {
    stop("'", name, "' must be a positive integer.")
  }
  invisible(as.integer(x))
}

#' @noRd
.build_cause_specific_data <- function(state_data, dests, covariates) {
  duration <- state_data$Tstop - state_data$Tstart
  if (anyNA(duration) || any(!is.finite(duration)) || any(duration <= 0)) {
    bad <- which(is.na(duration) | !is.finite(duration) | duration <= 0)[1L]
    stop("Invalid nonpositive sojourn duration for subject '",
         as.character(state_data$id[bad]), "' in state '",
         as.character(state_data$from[bad]), "'.")
  }
  event_type <- integer(nrow(state_data))
  for (d_idx in seq_along(dests)) {
    event_type[state_data$status == 1L &
                 as.character(state_data$to) == dests[d_idx]] <- d_idx
  }
  observed_unmapped <- state_data$status == 1L & event_type == 0L
  if (any(observed_unmapped)) {
    bad <- which(observed_unmapped)[1L]
    stop("Observed exit '", as.character(state_data$from[bad]), " -> ",
         as.character(state_data$to[bad]), "' is not a declared edge.")
  }
  data.frame(
    id = state_data$id,
    duration = duration,
    event_type = event_type,
    state_data[, covariates, drop = FALSE],
    check.names = FALSE
  )
}

#' @export
print.rfmstate <- function(x, ...) {
  cat("Clock-Reset Semi-Markov Random-Forest Model\n")
  cat("  Common initial state:", x$initial_state, "\n")
  cat("  Time scale: duration since fresh state entry\n")
  cat("  Covariates:", paste(x$covariates, collapse = ", "), "\n")
  cat("  Trees per edge:", x$params$num.trees, "\n")
  cat("  min_events safeguard:", x$params$min_events, "\n")
  cat("\nFitted edge models:\n")
  for (edge in names(x$edge_metadata)) {
    info <- x$edge_metadata[[edge]]
    cat("  ", edge, ": ", info$n_sojourns, " sojourns, ",
        info$n_events, " target events, OOB C = ",
        format(round(info$oob_concordance, 4), nsmall = 4), "\n", sep = "")
  }
  invisible(x)
}
