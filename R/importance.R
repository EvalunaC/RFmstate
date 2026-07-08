#' Feature Importance per Transition
#'
#' Extracts and organizes variable importance scores from the fitted random
#' forest models for each transition.
#'
#' @param object A fitted \code{rfmstate} model (must have been fit with
#'   \code{importance != "none"}).
#' @param ... Ignored.
#'
#' @return An object of class \code{"rfmstate_importance"} containing:
#'   \describe{
#'     \item{importance}{Data frame with columns \code{variable},
#'       \code{from}, \code{to}, \code{importance}.}
#'     \item{importance_matrix}{Matrix with variables as rows and transitions
#'       as columns.}
#'     \item{covariates}{Covariate names.}
#'     \item{transitions}{Character vector of transition labels.}
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
#' imp <- importance(fit)
#' print(imp)
#' }
#'
#' @export
importance <- function(object, ...) {
  UseMethod("importance")
}

#' @export
importance.rfmstate <- function(object, ...) {
  covariates <- object$covariates

  # Collect importance from each model
  imp_list <- list()
  transitions <- character(0)

  for (state_h in names(object$models)) {
    for (dest in names(object$models[[state_h]])) {
      rf_model <- object$models[[state_h]][[dest]]
      trans_label <- paste(state_h, "->", dest)
      transitions <- c(transitions, trans_label)

      vi <- rf_model$variable.importance
      if (is.null(vi)) {
        vi <- stats::setNames(rep(NA_real_, length(covariates)), covariates)
      }

      for (var in covariates) {
        val <- if (var %in% names(vi)) vi[var] else NA_real_
        imp_list[[length(imp_list) + 1]] <- data.frame(
          variable = var,
          from = state_h,
          to = dest,
          transition = trans_label,
          importance = as.numeric(val),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  imp_df <- do.call(rbind, imp_list)
  rownames(imp_df) <- NULL

  # Create matrix form
  imp_mat <- matrix(
    NA_real_,
    nrow = length(covariates),
    ncol = length(transitions),
    dimnames = list(covariates, transitions)
  )
  for (i in seq_len(nrow(imp_df))) {
    var <- imp_df$variable[i]
    trans <- imp_df$transition[i]
    imp_mat[var, trans] <- imp_df$importance[i]
  }

  structure(
    list(
      importance = imp_df,
      importance_matrix = imp_mat,
      covariates = covariates,
      transitions = transitions
    ),
    class = "rfmstate_importance"
  )
}

#' @export
print.rfmstate_importance <- function(x, ...) {
  cat("Feature Importance per Transition\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  # Print matrix form
  imp_mat <- x$importance_matrix
  # Round for display
  imp_display <- round(imp_mat, 4)
  print(imp_display)

  cat("\nTop variables per transition:\n")
  for (trans in x$transitions) {
    vals <- imp_mat[, trans]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      top <- names(sort(vals, decreasing = TRUE))[1]
      cat("  ", trans, ": ", top,
          " (", round(vals[top], 4), ")\n", sep = "")
    }
  }

  invisible(x)
}
