#' Summary of Random Forest Multistate Model
#'
#' Provides a comprehensive summary of the fitted model including per-origin
#' state information, OOB prediction error, and transition event counts.
#'
#' @param object A fitted \code{rfmstate} model.
#' @param ... Ignored.
#'
#' @return An object of class \code{"summary.rfmstate"} containing the stored
#'   call, covariates, forest parameters, patient/event/interval totals,
#'   structure, and a per-edge table of risk-set size, target events,
#'   competing/non-target exits, external censoring, ranger OOB error,
#'   ranger OOB concordance, and verified OOB coverage.
#'
#' @details The per-edge OOB concordance is \code{1 - prediction.error} from
#' ranger for that edge's binary cause-specific endpoint. The printed summary
#' also lists every forwarded ranger argument.
#'
#' @section Limitations:
#' Edge OOB values do not validate assembled full-state probabilities and may
#' be unstable for sparse transitions. Use \code{diagnose(..., method = "cv")}
#' for patient-level held-out full-state Brier scoring.
#'
#' @examples
#' \donttest{
#' ms <- clinical_states()
#' dat <- sim_clinical_data(n = 200, structure = ms, seed = 42)
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
#' fit <- rfmstate(msdata, num.trees = 100, seed = 42)
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
    n_target_events = integer(0),
    n_competing_exits = integer(0),
    n_external_censored = integer(0),
    oob_error = numeric(0),
    oob_concordance = numeric(0),
    oob_fraction = numeric(0),
    stringsAsFactors = FALSE
  )

  for (state_h in names(object$models)) {
    for (dest in names(object$models[[state_h]])) {
      rf_model <- object$models[[state_h]][[dest]]
      info <- object$edge_metadata[[paste0(state_h, "->", dest)]]
      n_total <- info$n_sojourns

      # OOB prediction error
      oob_err <- rf_model$prediction.error

      trans_summary <- rbind(trans_summary, data.frame(
        transition = paste(state_h, "->", dest),
        from = state_h,
        to = dest,
        n_total = n_total,
        n_target_events = info$n_events,
        n_competing_exits = info$n_competing_exits,
        n_external_censored = info$n_external_censored,
        oob_error = round(oob_err, 4),
        oob_concordance = round(1 - oob_err, 4),
        oob_fraction = info$oob_coverage$fraction_with_oob,
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

#' @rdname print_rfmstate_objects
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
  cat("  Effective OOB sampling: replace =",
      x$params$effective_ranger_args$replace,
      "; sample.fraction =",
      x$params$effective_ranger_args$sample.fraction,
      "; oob.error = TRUE; keep.inbag = TRUE\n")

  cat("\nTransition-specific models:\n")
  cat(paste(rep("-", 100), collapse = ""), "\n")
  cat(sprintf("%-20s %6s %7s %7s %7s %10s %10s %8s\n",
              "Transition", "Total", "Target", "Compete", "Extern",
              "OOB Error", "OOB C", "OOB Frac"))
  cat(paste(rep("-", 100), collapse = ""), "\n")

  for (i in seq_len(nrow(x$trans_summary))) {
    row <- x$trans_summary[i, ]
    cat(sprintf("%-20s %6d %7d %7d %7d %10.4f %10.4f %8.3f\n",
                row$transition, row$n_total, row$n_target_events,
                row$n_competing_exits, row$n_external_censored,
                row$oob_error, row$oob_concordance, row$oob_fraction))
  }
  cat(paste(rep("-", 100), collapse = ""), "\n")

  invisible(x)
}
