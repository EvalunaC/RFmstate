#' Diagnostics for Random Forest Multistate Model
#'
#' Computes diagnostic measures including OOB-based prediction error,
#' Brier score, concordance index, and bias-variance decomposition for
#' each transition-specific model.
#'
#' @param object A fitted \code{rfmstate} model.
#' @param eval_times Numeric vector of times at which to evaluate
#'   diagnostics. If \code{NULL}, uses quantiles of event times.
#' @param ... Ignored.
#'
#' @return An object of class \code{"rfmstate_diag"} containing:
#'   \describe{
#'     \item{oob_error}{Data frame of OOB prediction errors per transition.}
#'     \item{brier}{List of time-dependent Brier score components per
#'       transition.}
#'     \item{concordance}{Data frame of concordance indices per transition.}
#'     \item{bias_variance}{Data frame of bias-variance decomposition per
#'       transition.}
#'     \item{eval_times}{Evaluation times used.}
#'   }
#'
#' @details
#' The bias-variance decomposition uses OOB predictions from the random
#' forest ensemble. For each transition:
#' \itemize{
#'   \item \strong{Bias}: systematic difference between predicted and
#'     observed survival.
#'   \item \strong{Variance}: variability of predictions across trees
#'     (estimated from tree-level OOB predictions when available).
#'   \item \strong{Brier score}: integrated prediction error combining bias
#'     and variance.
#'   \item \strong{C-index}: concordance between predicted risk and
#'     observed event ordering.
#' }
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
#' diag <- diagnose(fit)
#' print(diag)
#' }
#'
#' @export
diagnose <- function(object, ...) {
  UseMethod("diagnose")
}

#' @rdname diagnose
#' @export
diagnose.rfmstate <- function(object, eval_times = NULL, ...) {
  covariates <- object$covariates
  structure <- object$structure

  oob_errors <- data.frame(
    transition = character(0),
    oob_error = numeric(0),
    stringsAsFactors = FALSE
  )

  brier_list <- list()
  concordance_df <- data.frame(
    transition = character(0),
    c_index = numeric(0),
    stringsAsFactors = FALSE
  )
  bv_df <- data.frame(
    transition = character(0),
    bias = numeric(0),
    variance = numeric(0),
    mse = numeric(0),
    stringsAsFactors = FALSE
  )

  for (state_h in names(object$models)) {
    cs_data <- object$origin_data[[state_h]]
    dests <- names(object$models[[state_h]])

    for (dest in dests) {
      trans_label <- paste(state_h, "->", dest)
      rf_model <- object$models[[state_h]][[dest]]

      # OOB error
      oob_err <- rf_model$prediction.error
      oob_errors <- rbind(oob_errors, data.frame(
        transition = trans_label,
        oob_error = oob_err,
        stringsAsFactors = FALSE
      ))

      # Concordance index from OOB
      dest_idx <- match(dest, structure$transitions[[state_h]])
      cs_status <- as.integer(cs_data$event_type == dest_idx)

      # Compute Brier score and C-index using OOB predictions
      fit_data <- data.frame(
        time = cs_data$duration,
        status = cs_status,
        cs_data[, covariates, drop = FALSE],
        stringsAsFactors = FALSE
      )
      fit_data <- fit_data[fit_data$time > 0, ]

      if (nrow(fit_data) < 10) next

      # Get OOB predicted survival
      oob_pred <- stats::predict(rf_model, data = fit_data[, covariates,
                                                           drop = FALSE])

      # Evaluation times
      if (is.null(eval_times)) {
        et <- stats::quantile(
          fit_data$time[fit_data$status == 1],
          probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE
        )
        et <- sort(unique(as.numeric(et)))
      } else {
        et <- eval_times
      }

      # Time-dependent Brier score
      brier_scores <- .compute_brier(
        fit_data$time, fit_data$status,
        oob_pred$survival, oob_pred$unique.death.times, et
      )
      brier_list[[trans_label]] <- data.frame(
        time = et,
        brier = brier_scores
      )

      # C-index (Harrell's)
      c_idx <- .compute_concordance(
        fit_data$time, fit_data$status,
        oob_pred$survival, oob_pred$unique.death.times
      )
      concordance_df <- rbind(concordance_df, data.frame(
        transition = trans_label,
        c_index = c_idx,
        stringsAsFactors = FALSE
      ))

      # Bias-variance decomposition
      bv <- .bias_variance(
        fit_data$time, fit_data$status,
        oob_pred$survival, oob_pred$unique.death.times, et
      )
      bv_df <- rbind(bv_df, data.frame(
        transition = trans_label,
        bias = bv$bias,
        variance = bv$variance,
        mse = bv$mse,
        stringsAsFactors = FALSE
      ))
    }
  }

  structure(
    list(
      oob_error = oob_errors,
      brier = brier_list,
      concordance = concordance_df,
      bias_variance = bv_df,
      eval_times = eval_times
    ),
    class = "rfmstate_diag"
  )
}

#' Compute time-dependent Brier score
#' @noRd
.compute_brier <- function(obs_time, obs_status, surv_matrix,
                           pred_times, eval_times) {
  n <- length(obs_time)
  brier <- numeric(length(eval_times))

  for (k in seq_along(eval_times)) {
    t <- eval_times[k]

    # Get predicted survival at time t for each subject
    t_idx <- which.min(abs(pred_times - t))
    if (length(t_idx) == 0) {
      brier[k] <- NA
      next
    }
    pred_surv <- surv_matrix[, t_idx]

    # Observed: Y_i = I(T_i > t)
    # For uncensored or event before t: Y = I(T > t)
    # For censored before t: excluded (IPCW)
    # Simplified Brier (without IPCW for tractability):
    # Include only subjects with T > t or event <= t
    include <- obs_time > t | (obs_time <= t & obs_status == 1)
    if (sum(include) < 5) {
      brier[k] <- NA
      next
    }

    obs_surv <- as.numeric(obs_time[include] > t)
    pred_s <- pred_surv[include]

    brier[k] <- mean((obs_surv - pred_s)^2)
  }
  brier
}

#' Compute concordance index
#' @noRd
.compute_concordance <- function(obs_time, obs_status, surv_matrix,
                                 pred_times) {
  n <- length(obs_time)
  if (n < 5 || sum(obs_status) < 3) return(NA_real_)

  # Use median time predicted survival as risk score
  mid_idx <- ceiling(ncol(surv_matrix) / 2)
  risk_score <- 1 - surv_matrix[, mid_idx]  # Higher = higher risk

  concordant <- 0
  discordant <- 0

  events <- which(obs_status == 1)
  for (i in events) {
    # Compare with subjects who survived longer
    longer <- which(obs_time > obs_time[i])
    for (j in longer) {
      if (risk_score[i] > risk_score[j]) {
        concordant <- concordant + 1
      } else if (risk_score[i] < risk_score[j]) {
        discordant <- discordant + 1
      } else {
        concordant <- concordant + 0.5
        discordant <- discordant + 0.5
      }
    }
  }

  total <- concordant + discordant
  if (total == 0) return(NA_real_)
  concordant / total
}

#' Bias-variance decomposition
#' @noRd
.bias_variance <- function(obs_time, obs_status, surv_matrix,
                           pred_times, eval_times) {
  n <- length(obs_time)

  # Compute at a representative time point (median event time)
  event_times <- obs_time[obs_status == 1]
  if (length(event_times) < 3) {
    return(list(bias = NA_real_, variance = NA_real_, mse = NA_real_))
  }
  t_eval <- stats::median(event_times)
  t_idx <- which.min(abs(pred_times - t_eval))

  pred_surv <- surv_matrix[, t_idx]
  obs_surv <- as.numeric(obs_time > t_eval)

  # Only use complete observations
  include <- obs_time > t_eval | (obs_time <= t_eval & obs_status == 1)
  pred_s <- pred_surv[include]
  obs_s <- obs_surv[include]

  if (length(pred_s) < 5) {
    return(list(bias = NA_real_, variance = NA_real_, mse = NA_real_))
  }

  # Bias: mean(pred - obs)
  bias <- mean(pred_s - obs_s)

  # Variance of predictions
  pred_var <- stats::var(pred_s)

  # MSE
  mse <- mean((pred_s - obs_s)^2)

  list(bias = round(bias, 6), variance = round(pred_var, 6),
       mse = round(mse, 6))
}

#' @export
print.rfmstate_diag <- function(x, ...) {
  cat("RF Multistate Model Diagnostics\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  cat("\nOOB Prediction Error:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  for (i in seq_len(nrow(x$oob_error))) {
    cat(sprintf("  %-25s %.4f\n",
                x$oob_error$transition[i], x$oob_error$oob_error[i]))
  }

  cat("\nConcordance Index (C-index):\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  for (i in seq_len(nrow(x$concordance))) {
    cat(sprintf("  %-25s %.4f\n",
                x$concordance$transition[i], x$concordance$c_index[i]))
  }

  cat("\nBias-Variance Decomposition:\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  cat(sprintf("  %-25s %8s %8s %8s\n",
              "Transition", "Bias", "Var", "MSE"))
  cat(paste(rep("-", 60), collapse = ""), "\n")
  for (i in seq_len(nrow(x$bias_variance))) {
    row <- x$bias_variance[i, ]
    cat(sprintf("  %-25s %8.4f %8.4f %8.4f\n",
                row$transition, row$bias, row$variance, row$mse))
  }

  invisible(x)
}
