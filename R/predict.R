#' Predict Entry-Conditioned State Probabilities
#'
#' Predicts clock-reset semi-Markov state-occupation probabilities conditional
#' on fresh entry into a selected state at elapsed duration zero.
#'
#' @param object A fitted \code{rfmstate} model.
#' @param newdata Prediction profiles. \code{NULL} uses training profiles and
#'   labels the result apparent/training.
#' @param times Finite nonnegative elapsed durations. \code{NULL} uses fitted
#'   edge event times that fall within the conservative support horizon.
#' @param s Legacy start-time argument; only zero is supported.
#' @param start_state Fresh-entry starting state. Defaults to the common initial
#'   state.
#' @param grid_step Optional positive initial internal grid step.
#' @param target_grid_points Initial regular-grid interval count when
#'   \code{grid_step} is \code{NULL}.
#' @param max_grid_points Maximum allowed internal grid-point count.
#' @param grid_tol Maximum change permitted between successive refinements.
#' @param max_grid_refinements Maximum number of grid-step halvings.
#' @param check_grid Whether convergence is required. Keep \code{TRUE} for
#'   inferential or reported work.
#' @param extrapolate \code{"error"} (default) or explicit \code{"flat"}
#'   cumulative-hazard sensitivity extension.
#' @param ... Ignored; no ranger arguments are accepted during prediction.
#'
#' @return An \code{rfmstate_pred} object. \code{entry_prob} and the temporary
#'   alias \code{P} have dimensions profile by selected starting state by
#'   occupied state by elapsed time; the starting-state dimension has length
#'   one. \code{state_occ} is that selected starting-state slice.
#'   The object also contains requested \code{time}, per-profile edge
#'   \code{cum_hazard}, validated \code{newdata}, used columns, support and
#'   extrapolation metadata, grid-convergence records, conditioning statement,
#'   predictor-support extrapolation records, fit-specific schema identifier,
#'   state structure, and package versions.
#'
#' @details These are not general Markov \eqn{P(s,t)} matrices. Predictions do
#'   not condition on an already elapsed sojourn and carry no confidence bands.
#'   By default the horizon cannot exceed the minimum observed sojourn support
#'   among reachable transient states.
#'
#' @section Limitations:
#' Prediction is conditional on fresh entry into \code{start_state} at elapsed
#' duration zero. General \eqn{P(s,t)} prediction for \code{s > 0}, left
#' truncation, ongoing-sojourn/landmark prediction, recurrent histories, and
#' confidence intervals are unavailable. Missing covariates, nonfinite numeric
#' values, incompatible classes, and unseen factor levels are rejected.
#' Numeric values outside the fit-specific range warn and are recorded as
#' predictor-support extrapolation; they are never truncated.
#' Times beyond support fail unless \code{extrapolate = "flat"}; that option
#' assumes zero additional hazard and is a sensitivity analysis only.
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
#' horizon <- min(fit$max_duration_by_origin)
#' predict(fit, data.frame(age=60, sex=1, BMI=26, treatment=1),
#'         times = c(0, horizon / 2))
#' }
#'
#' @export
predict.rfmstate <- function(object, newdata = NULL, times = NULL, s = 0,
                             start_state = NULL, grid_step = NULL,
                             target_grid_points = 1024L,
                             max_grid_points = 131073L,
                             grid_tol = 5e-4,
                             max_grid_refinements = 7L,
                             check_grid = TRUE,
                             extrapolate = c("error", "flat"), ...) {
  if (!inherits(object, "rfmstate")) stop("'object' must be an rfmstate fit.")
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s != 0) {
    stop("Only s = 0 fresh-entry prediction is supported.")
  }
  structure_obj <- object$structure
  if (is.null(start_state)) start_state <- object$initial_state
  if (!is.character(start_state) || length(start_state) != 1L ||
      is.na(start_state) || !(start_state %in% structure_obj$state_names)) {
    stop("'start_state' must be one state in the fitted structure.")
  }
  extrapolate <- match.arg(extrapolate)
  prediction_type <- if (is.null(newdata)) "apparent/training" else "new profiles"
  if (is.null(newdata)) {
    first <- !duplicated(object$msdata$id)
    newdata <- object$msdata[first, object$covariates, drop = FALSE]
    class(newdata) <- "data.frame"
  }
  newdata <- .validate_prediction_data(newdata, object)
  predictor_extrapolation <- attr(newdata, "predictor_extrapolation")
  n_subjects <- nrow(newdata)

  reachable <- intersect(.reachable_states(structure_obj, start_state),
                         structure_obj$transient)
  support <- if (start_state %in% structure_obj$absorbing) {
    Inf
  } else {
    min(object$max_duration_by_origin[reachable])
  }
  if (is.null(times)) {
    times <- sort(unique(c(0, unlist(object$event_times))))
    if (is.finite(support)) times <- times[times <= support]
  }
  if (!is.numeric(times) || !length(times) || anyNA(times) ||
      any(!is.finite(times)) || any(times < 0)) {
    stop("'times' must contain finite, nonnegative elapsed durations.")
  }
  times <- sort(unique(as.numeric(times)))
  if (is.finite(support) && max(times) > support && extrapolate == "error") {
    stop("Requested horizon exceeds conservative support ", format(support),
         " from starting state '", start_state, "'.")
  }

  ns <- structure_obj$n_states
  nt <- length(times)
  state_names <- structure_obj$state_names
  entry_array <- array(
    0,
    dim = c(n_subjects, 1L, ns, nt),
    dimnames = list(
      profile = rownames(newdata),
      starting_state = start_state,
      occupied_state = state_names,
      elapsed_time = format(times, scientific = FALSE, trim = TRUE)
    )
  )
  occ_array <- array(
    0,
    dim = c(n_subjects, ns, nt),
    dimnames = list(
      profile = rownames(newdata), occupied_state = state_names,
      elapsed_time = format(times, scientific = FALSE, trim = TRUE)
    )
  )
  hazard_profiles <- vector("list", n_subjects)
  grid_metadata <- vector("list", n_subjects)
  for (i in seq_len(n_subjects)) {
    profile <- newdata[i, object$covariates, drop = FALSE]
    curves <- .predict_subject_hazards(object, profile)
    hazard_profiles[[i]] <- curves
    tp <- compute_trans_prob(
      curves, structure_obj, s = 0, times = times,
      start_state = start_state, grid_step = grid_step,
      target_grid_points = target_grid_points,
      max_grid_points = max_grid_points, grid_tol = grid_tol,
      max_grid_refinements = max_grid_refinements,
      check_grid = check_grid, extrapolate = extrapolate
    )
    entry_array[i, 1L, , ] <- tp$entry_prob[1L, , ]
    occ_array[i, , ] <- t(tp$state_occ)
    grid_metadata[[i]] <- list(
      grid = tp$grid, grid_step = tp$grid_step,
      grid_points = tp$grid_points, grid_converged = tp$grid_converged,
      grid_checked = tp$grid_checked,
      grid_error = tp$grid_error, refinements = tp$grid_refinements,
      hazard_roundoff_corrections = tp$hazard_roundoff_corrections
    )
  }

  structure(
    list(
      time = times,
      entry_prob = entry_array,
      P = entry_array,
      state_occ = occ_array,
      cum_hazard = hazard_profiles,
      structure = structure_obj,
      newdata = newdata,
      used_columns = object$covariates,
      predictor_extrapolation = predictor_extrapolation,
      n_subjects = n_subjects,
      start_state = start_state,
      initial_state = object$initial_state,
      time_scale = "clock-reset",
      process_assumption = "semi-Markov",
      process_model = "semi-Markov",
      history_summary = paste(
        "current state, duration since entry, and recorded baseline covariates"
      ),
      prediction_condition = paste(
        "fresh entry into the selected starting state at duration zero"
      ),
      conditioning = "fresh entry at elapsed duration zero",
      prediction_type = prediction_type,
      support_horizon = support,
      support_by_origin = object$max_duration_by_origin[reachable],
      predictor_schema_id = object$predictor_schema_id,
      extrapolation = if (max(times) > support) {
        "flat cumulative hazard"
      } else {
        "none"
      },
      grid_metadata = grid_metadata,
      package_versions = object$package_versions
    ),
    class = "rfmstate_pred"
  )
}

#' @noRd
.validate_prediction_data <- function(newdata, object) {
  if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.")
  missing_covs <- setdiff(object$covariates, names(newdata))
  if (length(missing_covs)) {
    stop("Missing covariates in newdata: ", paste(missing_covs, collapse = ", "))
  }
  out <- newdata[, object$covariates, drop = FALSE]
  extrapolation <- data.frame(
    predictor = character(0), profile = character(0), value = numeric(0),
    training_min = numeric(0), training_max = numeric(0),
    stringsAsFactors = FALSE
  )
  for (nm in object$covariates) {
    x <- out[[nm]]
    schema <- object$predictor_schema[[nm]]
    if (anyNA(x)) stop("Prediction covariate '", nm, "' contains missing values.")
    trained_factor <- !is.null(schema$levels)
    if (trained_factor) {
      if (!is.factor(x)) {
        stop("Prediction covariate '", nm,
             "' must be a factor compatible with training data.")
      }
      unseen <- setdiff(unique(as.character(x)), schema$levels)
      if (length(unseen)) {
        stop("Unseen factor level(s) for '", nm, "': ",
             paste(unseen, collapse = ", "), ".")
      }
      out[[nm]] <- factor(as.character(x), levels = schema$levels,
                          ordered = isTRUE(schema$ordered))
    } else if ("character" %in% schema$class) {
      stop("Stored character predictor schema is unsupported; refit after prepare_data().")
    } else if (any(schema$class %in% c("numeric", "integer"))) {
      if (!is.numeric(x) || any(!is.finite(x))) {
        stop("Prediction covariate '", nm, "' must be finite numeric data.")
      }
      outside <- x < schema$range[1L] | x > schema$range[2L]
      if (any(outside)) {
        extrapolation <- rbind(extrapolation, data.frame(
          predictor = nm,
          profile = rownames(out)[outside] %||% as.character(which(outside)),
          value = as.numeric(x[outside]),
          training_min = schema$range[1L],
          training_max = schema$range[2L],
          stringsAsFactors = FALSE
        ))
      }
    } else if ("logical" %in% schema$class && !is.logical(x)) {
      stop("Prediction covariate '", nm, "' must be logical.")
    }
  }
  if (nrow(extrapolation)) {
    warning(
      "Prediction contains numeric value(s) outside fit-specific range for: ",
      paste(unique(extrapolation$predictor), collapse = ", "),
      ". Values were not truncated.",
      call. = FALSE
    )
  }
  attr(out, "predictor_extrapolation") <- extrapolation
  out
}

#' Predict all edge cumulative hazards for one profile
#' @noRd
.predict_subject_hazards <- function(object, profile) {
  trans_list <- object$structure$trans_list
  expected <- paste0(trans_list$from, "->", trans_list$to)
  curves <- stats::setNames(vector("list", length(expected)), expected)
  for (k in seq_len(nrow(trans_list))) {
    from <- trans_list$from[k]
    to <- trans_list$to[k]
    edge <- expected[k]
    model <- object$models[[from]][[to]]
    if (is.null(model)) {
      stop("Missing fitted model for declared edge '", edge, "'.")
    }
    pred <- stats::predict(model, data = profile)
    times <- pred$unique.death.times
    chf <- pred$chf
    if (is.null(chf)) {
      stop("Ranger did not return predicted cumulative hazard for edge '",
           edge, "'.")
    }
    values <- if (is.matrix(chf)) chf[1L, ] else as.numeric(chf)
    if (length(values) != length(times) || any(!is.finite(values)) ||
        any(values < 0)) {
      stop("Invalid ranger cumulative-hazard prediction for edge '", edge, "'.")
    }
    support <- object$max_duration_by_origin[[from]]
    curve_times <- c(0, times)
    curve_values <- c(0, values)
    keep <- !duplicated(curve_times, fromLast = TRUE)
    curve_times <- curve_times[keep]
    curve_values <- curve_values[keep]
    if (support > max(curve_times)) {
      curve_times <- c(curve_times, support)
      curve_values <- c(curve_values, tail(curve_values, 1L))
    }
    curves[[edge]] <- data.frame(time = curve_times, hazard = curve_values)
  }
  curves
}

#' @rdname print_rfmstate_objects
#' @export
print.rfmstate_pred <- function(x, ...) {
  cat("Entry-Conditioned RFmstate Predictions\n")
  cat("  Profiles:", x$n_subjects, "\n")
  cat("  Starting state:", x$start_state, "\n")
  cat("  Conditioning: fresh entry at duration zero\n")
  cat("  Time scale: clock-reset\n")
  cat("  Prediction type:", x$prediction_type, "\n")
  cat("  Elapsed-time range: [", min(x$time), ", ", max(x$time), "]\n", sep = "")
  cat("  Extrapolation:", x$extrapolation, "\n")
  invisible(x)
}
