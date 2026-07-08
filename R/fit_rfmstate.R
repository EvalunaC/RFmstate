#' Fit Random Forest Multistate Model
#'
#' Fits transition-specific cause-specific random survival forests for
#' multistate survival analysis. For each transient origin state, a competing
#' risks model is fit using random forests, where the competing events are
#' the possible transitions to destination states.
#'
#' @param msdata An \code{msdata} object from \code{\link{prepare_data}}.
#' @param covariates Character vector of covariate column names to use as
#'   predictors. If \code{NULL}, all non-structural columns are used.
#' @param num.trees Integer, number of trees per forest (default 1000).
#' @param mtry Integer, number of variables to try at each split. Default
#'   \code{NULL} uses \code{floor(sqrt(p))} where p is number of covariates.
#' @param min.node.size Integer, minimum node size (default 15).
#' @param importance Character, variable importance mode. One of
#'   \code{"permutation"} (default), \code{"impurity"}, or \code{"none"}.
#' @param seed Integer, random seed for reproducibility (default \code{NULL}).
#' @param ... Additional arguments passed to \code{\link[ranger]{ranger}}.
#'
#' @return An object of class \code{"rfmstate"} containing:
#'   \describe{
#'     \item{models}{Named list of fitted \code{ranger} objects, one per
#'       origin state.}
#'     \item{structure}{The multistate structure.}
#'     \item{covariates}{Character vector of covariate names used.}
#'     \item{origin_data}{Named list of per-origin-state data subsets.}
#'     \item{event_times}{Named list of unique event times per origin state.}
#'     \item{call}{The matched call.}
#'     \item{params}{List of tuning parameters used.}
#'   }
#'
#' @details
#' For each transient state \eqn{h}, the method:
#' \enumerate{
#'   \item Subsets all intervals where the patient is in state \eqn{h}.
#'   \item Defines time as the duration in state \eqn{h} (Tstop - Tstart).
#'   \item Codes competing events: 0 = censored, 1, 2, ... for each possible
#'     destination state.
#'   \item Fits a cause-specific random survival forest using
#'     \code{\link[ranger]{ranger}} with \code{survival} tree type.
#' }
#'
#' Transition probabilities are then computed by combining per-origin-state
#' predicted cumulative hazards via the product-integral formula.
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
#' fit <- rfmstate(msdata, covariates = c("age", "sex", "BMI", "treatment"))
#' print(fit)
#' }
#'
#' @export
rfmstate <- function(msdata, covariates = NULL, num.trees = 1000L,
                     mtry = NULL, min.node.size = 15L,
                     importance = "permutation", seed = NULL, ...) {
  if (!inherits(msdata, "msdata")) {
    stop("'msdata' must be an 'msdata' object from prepare_data().")
  }

  cl <- match.call()
  structure <- attr(msdata, "structure")

  # Determine covariates
  if (is.null(covariates)) {
    structural_cols <- c("id", "from", "to", "Tstart", "Tstop",
                         "status", "trans_id", "duration")
    covariates <- setdiff(names(msdata), structural_cols)
  }

  # Validate covariates exist
  missing_covs <- setdiff(covariates, names(msdata))
  if (length(missing_covs) > 0) {
    stop("Covariates not found in msdata: ",
         paste(missing_covs, collapse = ", "))
  }

  if (is.null(mtry)) {
    mtry <- max(1L, floor(sqrt(length(covariates))))
  }

  transient_states <- structure$transient
  models <- list()
  origin_data <- list()
  event_times <- list()

  for (state_h in transient_states) {
    # Get destinations from this state
    dests <- structure$transitions[[state_h]]
    if (is.null(dests) || length(dests) == 0) next

    # Subset data: intervals where patient is in state_h
    state_data <- msdata[msdata$from == state_h, ]
    if (nrow(state_data) == 0) {
      message("No data for origin state '", state_h,
              "'. Skipping.")
      next
    }

    # Build cause-specific dataset
    cs_data <- .build_cause_specific_data(state_data, dests, covariates)

    if (is.null(cs_data) || nrow(cs_data) < 2 * min.node.size) {
      message("Insufficient data for origin state '", state_h,
              "' (n=", nrow(cs_data), "). Skipping.")
      next
    }

    origin_data[[state_h]] <- cs_data

    # Fit cause-specific random forests
    # For each cause (destination), fit a separate survival forest
    # treating other events as censored
    state_models <- list()

    for (d_idx in seq_along(dests)) {
      dest <- dests[d_idx]

      # Create cause-specific dataset: event = 1 if transition to dest,
      # 0 otherwise (censored or other cause)
      cs_status <- as.integer(cs_data$event_type == d_idx)

      # Need at least some events
      if (sum(cs_status) < 3) {
        message("  Too few events for ", state_h, " -> ", dest,
                " (n_events=", sum(cs_status), "). Skipping.")
        next
      }

      # Fit ranger survival forest
      fit_data <- data.frame(
        time = cs_data$duration,
        status = cs_status,
        cs_data[, covariates, drop = FALSE],
        check.names = FALSE
      )

      # Remove rows with zero or negative duration
      fit_data <- fit_data[fit_data$time > 0, ]

      if (nrow(fit_data) < 2 * min.node.size) next

      # Build formula explicitly to avoid using time/status as predictors
      surv_formula <- stats::as.formula(
        paste("survival::Surv(time, status) ~",
              paste(covariates, collapse = " + "))
      )

      rf_fit <- ranger::ranger(
        formula = surv_formula,
        data = fit_data,
        num.trees = num.trees,
        mtry = min(mtry, length(covariates)),
        min.node.size = min.node.size,
        importance = importance,
        seed = seed,
        ...
      )

      state_models[[dest]] <- rf_fit
    }

    if (length(state_models) > 0) {
      models[[state_h]] <- state_models
      event_times[[state_h]] <- sort(unique(
        cs_data$duration[cs_data$event_type > 0]
      ))
    }
  }

  if (length(models) == 0) {
    stop("No models could be fit. Check data quality and sample sizes.")
  }

  result <- structure(
    list(
      models = models,
      structure = structure,
      covariates = covariates,
      origin_data = origin_data,
      event_times = event_times,
      msdata = msdata,
      call = cl,
      params = list(
        num.trees = num.trees,
        mtry = mtry,
        min.node.size = min.node.size,
        importance = importance
      )
    ),
    class = "rfmstate"
  )
  result
}

#' Build cause-specific data for one origin state
#' @noRd
.build_cause_specific_data <- function(state_data, dests, covariates) {
  # Duration in this state
  duration <- state_data$Tstop - state_data$Tstart

  # Event type: 0 = censored, 1..K = transition to dest 1..K
  event_type <- integer(nrow(state_data))
  for (d_idx in seq_along(dests)) {
    mask <- state_data$status == 1 & state_data$to == dests[d_idx]
    mask[is.na(mask)] <- FALSE
    event_type[mask] <- d_idx
  }

  result <- data.frame(
    id = state_data$id,
    duration = duration,
    event_type = event_type,
    state_data[, covariates, drop = FALSE],
    stringsAsFactors = FALSE
  )

  # Remove invalid rows
 result <- result[result$duration > 0 & !is.na(result$duration), ]
  result
}

#' @export
print.rfmstate <- function(x, ...) {
  cat("Random Forest Multistate Model\n")
  cat("Call: ")
  print(x$call)
  cat("\nCovariates:", paste(x$covariates, collapse = ", "), "\n")
  cat("Parameters:\n")
  cat("  num.trees:", x$params$num.trees, "\n")
  cat("  mtry:", x$params$mtry, "\n")
  cat("  min.node.size:", x$params$min.node.size, "\n")

  cat("\nModels fitted per origin state:\n")
  for (state_h in names(x$models)) {
    dests <- names(x$models[[state_h]])
    n_data <- nrow(x$origin_data[[state_h]])
    cat("  ", state_h, " (n=", n_data, "): -> ",
        paste(dests, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
