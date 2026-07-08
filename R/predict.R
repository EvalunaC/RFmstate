#' Predict Transition Probabilities for New Data
#'
#' Predicts patient-specific transition probability matrices and state
#' occupation probabilities using fitted random forest multistate models.
#'
#' @param object A fitted \code{rfmstate} model.
#' @param newdata A data frame with the same covariates used in fitting.
#'   If \code{NULL}, predictions are made for the training data.
#' @param times Numeric vector of times at which to compute transition
#'   probabilities. If \code{NULL}, uses all unique event times.
#' @param s Numeric, starting time (default 0).
#' @param ... Ignored.
#'
#' @return An object of class \code{"rfmstate_pred"} containing:
#'   \describe{
#'     \item{time}{Evaluation times.}
#'     \item{P}{Array of transition probability matrices (n_subjects x
#'       n_states x n_states x n_times).}
#'     \item{state_occ}{Array of state occupation probabilities (n_subjects x
#'       n_states x n_times).}
#'     \item{cum_hazard}{List of per-subject cumulative hazard matrices.}
#'     \item{structure}{The multistate structure.}
#'     \item{newdata}{The prediction data.}
#'   }
#'
#' @examples
#' \donttest{
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
#' fit <- rfmstate(msdata, covariates = c("age", "sex", "BMI", "treatment"),
#'                 num.trees = 100)
#' newpat <- data.frame(age = c(50, 70), sex = c(0, 1),
#'                      BMI = c(25, 30), treatment = c(1, 0))
#' pred <- predict(fit, newdata = newpat, times = c(30, 90, 180, 365))
#' }
#'
#' @export
predict.rfmstate <- function(object, newdata = NULL, times = NULL,
                             s = 0, ...) {
  structure <- object$structure
  state_names <- structure$state_names
  ns <- structure$n_states
  covariates <- object$covariates

  if (is.null(newdata)) {
    # Use unique covariate profiles from training data
    first_rows <- object$msdata[!duplicated(object$msdata$id), ]
    newdata <- first_rows[, covariates, drop = FALSE]
  }

  # Validate newdata
  missing_covs <- setdiff(covariates, names(newdata))
  if (length(missing_covs) > 0) {
    stop("Missing covariates in newdata: ",
         paste(missing_covs, collapse = ", "))
  }

  n_subj <- nrow(newdata)

  # Get all unique event times across all origin states
  all_event_times <- sort(unique(unlist(object$event_times)))
  if (!is.null(times)) {
    eval_times <- sort(unique(times[times > s]))
  } else {
    eval_times <- all_event_times[all_event_times > s]
  }

  n_times <- length(eval_times)

  # For each subject, compute predicted cause-specific cumulative hazards
  # from each origin-state model, then compute transition probabilities
  P_array <- array(0, dim = c(n_subj, ns, ns, n_times),
                   dimnames = list(NULL, state_names, state_names, NULL))
  occ_array <- array(0, dim = c(n_subj, ns, n_times),
                     dimnames = list(NULL, state_names, NULL))

  for (subj in seq_len(n_subj)) {
    subj_data <- newdata[subj, covariates, drop = FALSE]

    # Get cumulative hazards for each transition from RF predictions
    cum_hazards <- .predict_subject_hazards(
      object, subj_data, eval_times, structure
    )

    # Compute transition probabilities via product-integral
    tp <- .product_integral(cum_hazards, eval_times, structure, s)

    P_array[subj, , , ] <- tp$P
    occ_array[subj, , ] <- tp$occ
  }

  structure(
    list(
      time = eval_times,
      P = P_array,
      state_occ = occ_array,
      structure = structure,
      newdata = newdata,
      n_subjects = n_subj
    ),
    class = "rfmstate_pred"
  )
}

#' Predict cumulative hazards for one subject across all transitions
#' @noRd
.predict_subject_hazards <- function(object, subj_data, times, structure) {
  trans_list <- structure$trans_list
  cum_hazards <- vector("list", nrow(trans_list))

  for (tr_idx in seq_len(nrow(trans_list))) {
    from <- trans_list$from[tr_idx]
    to <- trans_list$to[tr_idx]

    # Check if we have a model for this transition
    if (is.null(object$models[[from]]) ||
        is.null(object$models[[from]][[to]])) {
      # No model: assume zero hazard
      cum_hazards[[tr_idx]] <- data.frame(
        time = times, hazard = rep(0, length(times))
      )
      next
    }

    rf_model <- object$models[[from]][[to]]

    # Predict survival function for this subject
    pred <- stats::predict(rf_model, data = subj_data)

    # ranger returns survival function at unique death times
    surv_times <- pred$unique.death.times
    # Handle both single subject (vector) and multiple subjects (matrix)
    surv_mat <- pred$survival
    if (is.matrix(surv_mat)) {
      surv_fn <- surv_mat[1, ]
    } else {
      surv_fn <- surv_mat
    }

    # Convert survival to cumulative hazard: H(t) = -log(S(t))
    surv_fn[surv_fn <= 0] <- .Machine$double.eps
    chf <- -log(surv_fn)

    # Interpolate to evaluation times
    chf_interp <- stats::approx(
      surv_times, chf, xout = times,
      method = "constant", rule = 2, f = 0
    )$y

    cum_hazards[[tr_idx]] <- data.frame(
      time = times, hazard = chf_interp
    )
  }

  cum_hazards
}

#' Product-integral for one subject
#' @noRd
.product_integral <- function(cum_hazards, times, structure, s) {
  state_names <- structure$state_names
  ns <- structure$n_states
  trans_list <- structure$trans_list
  n_times <- length(times)

  P_array <- array(0, dim = c(ns, ns, n_times))
  occ_mat <- matrix(0, nrow = ns, ncol = n_times)

  P_current <- diag(ns)
  p0 <- rep(0, ns)
  p0[1] <- 1

  # Get hazard increments
  haz_inc <- lapply(cum_hazards, function(ch) {
    c(ch$hazard[1], diff(ch$hazard))
  })

  for (k in seq_len(n_times)) {
    # Build increment matrix
    dA <- matrix(0, nrow = ns, ncol = ns)

    for (tr_idx in seq_len(nrow(trans_list))) {
      fi <- match(trans_list$from[tr_idx], state_names)
      ti <- match(trans_list$to[tr_idx], state_names)
      dA[fi, ti] <- haz_inc[[tr_idx]][k]
    }

    # Diagonal
    for (h in seq_len(ns)) {
      dA[h, h] <- -sum(dA[h, -h])
    }

    # Product-integral step
    P_current <- P_current %*% (diag(ns) + dA)
    P_current[P_current < 0] <- 0

    # Normalize rows
    rs <- rowSums(P_current)
    for (h in seq_len(ns)) {
      if (rs[h] > 0) P_current[h, ] <- P_current[h, ] / rs[h]
    }

    P_array[, , k] <- P_current
    occ_mat[, k] <- as.vector(p0 %*% P_current)
  }

  list(P = P_array, occ = occ_mat)
}

#' @export
print.rfmstate_pred <- function(x, ...) {
  cat("RF Multistate Predictions\n")
  cat("  Subjects:", x$n_subjects, "\n")
  cat("  Time points:", length(x$time), "\n")
  cat("  Time range: [", min(x$time), ", ", max(x$time), "]\n", sep = "")
  cat("  States:", paste(x$structure$state_names, collapse = ", "), "\n")
  invisible(x)
}
