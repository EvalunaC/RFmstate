#' Summary of Random Forest Multistate Model
#'
#' Provides a comprehensive summary of the fitted model including per-origin
#' state information, OOB prediction error, and transition event counts.
#'
#' @param object A fitted \code{rfmstate} model.
#' @param ... Ignored.
#'
#' @return An object of class \code{"summary.rfmstate"}, printed invisibly.
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
#' summary(fit)
#' }
#'
#' @export
summary.rfmstate <- function(object, ...) {
  structure_obj <- object$structure
  trans_list <- structure_obj$trans_list

  # Per-transition summary
  trans_summary <- data.frame(
    transition = character(0),
    from = character(0),
    to = character(0),
    n_total = integer(0),
    n_events = integer(0),
    n_censored = integer(0),
    oob_error = numeric(0),
    oob_concordance = numeric(0),
    stringsAsFactors = FALSE
  )

  for (state_h in names(object$models)) {
    for (dest in names(object$models[[state_h]])) {
      rf_model <- object$models[[state_h]][[dest]]
      cs_data <- object$origin_data[[state_h]]

      n_total <- nrow(cs_data)
      dest_idx <- match(dest, structure_obj$transitions[[state_h]])
      n_events <- sum(cs_data$event_type == dest_idx, na.rm = TRUE)
      n_censored <- n_total - n_events

      # OOB prediction error
      oob_err <- rf_model$prediction.error

      trans_summary <- rbind(trans_summary, data.frame(
        transition = paste(state_h, "->", dest),
        from = state_h,
        to = dest,
        n_total = n_total,
        n_events = n_events,
        n_censored = n_censored,
        oob_error = round(oob_err, 4),
        oob_concordance = round(1 - oob_err, 4),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Overall summary
  total_patients <- length(unique(object$msdata$id))
  total_events <- sum(object$msdata$status == 1)
  total_intervals <- nrow(object$msdata)

  result <- structure(
    list(
      call = object$call,
      covariates = object$covariates,
      params = object$params,
      n_patients = total_patients,
      n_events = total_events,
      n_intervals = total_intervals,
      trans_summary = trans_summary,
      structure = structure_obj
    ),
    class = "summary.rfmstate"
  )

  result
}

#' @export
print.summary.rfmstate <- function(x, ...) {
  cat("Random Forest Multistate Model Summary\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  cat("\nCall: ")
  print(x$call)

  cat("\nData:\n")
  cat("  Patients:", x$n_patients, "\n")
  cat("  Total transitions:", x$n_events, "\n")
  cat("  Total intervals:", x$n_intervals, "\n")

  cat("\nCovariates:", paste(x$covariates, collapse = ", "), "\n")
  cat("  Time scale: clock-reset duration; semi-Markov assembly\n")

  cat("\nForest parameters:\n")
  cat("  Trees:", x$params$num.trees, "\n")
  cat("  mtry:", x$params$mtry, "\n")
  cat("  Min node size:", x$params$min.node.size, "\n")
  cat("  min_events safeguard:", x$params$min_events, "\n")
  cat("  Forwarded ranger arguments:",
      if (length(x$params$ranger_args)) {
        paste(names(x$params$ranger_args), collapse = ", ")
      } else {
        "none"
      }, "\n")

  cat("\nTransition-specific models:\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-25s %6s %6s %6s %10s %10s\n",
              "Transition", "Total", "Events", "Cens", "OOB Error", "OOB C"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in seq_len(nrow(x$trans_summary))) {
    row <- x$trans_summary[i, ]
    cat(sprintf("%-25s %6d %6d %6d %10.4f %10.4f\n",
                row$transition, row$n_total, row$n_events,
                row$n_censored, row$oob_error, row$oob_concordance))
  }
  cat(paste(rep("-", 70), collapse = ""), "\n")

  invisible(x)
}
