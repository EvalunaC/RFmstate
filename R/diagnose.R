#' Diagnostics and Patient-Level Validation for RFmstate
#'
#' Reports genuine ranger edge-specific OOB concordance. Optionally, complete
#' patient-level cross-validation refits every edge forest and returns IPCW
#' full-state Brier scores and an integrated Brier score.
#'
#' @param object A fitted \code{rfmstate} model.
#' @param eval_times Prespecified calendar-time validation horizons. For
#'   \code{method="cv"}, defaults to nine points inside 80 percent of the
#'   fitted conservative support.
#' @param method \code{"edge_oob"} (default) or \code{"cv"}.
#' @param folds Number of patient-level folds for cross-validation.
#' @param repeats Number of repeated fold assignments.
#' @param seed Fold-assignment and refit seed.
#' @param g_min Minimum allowed training-fold censoring survival.
#' @param ... Ignored.
#'
#' @return An \code{rfmstate_diag} object. Edge OOB tables are always present.
#'   With cross-validation, \code{brier} is a full-state IPCW score table and
#'   \code{ibs} is its trapezoidal integral over the reported interval.
#'
#' @details Ranger OOB concordance is \code{1 - prediction.error} for each
#'   binary cause-specific edge endpoint. It does not validate the assembled
#'   multistate probability vector. Full-state scores use subject-level held-out
#'   predictions and training-fold Kaplan--Meier censoring estimates. The
#'   stronger marginal independent-censoring assumption applies to these KM
#'   weights. No bias--variance decomposition is provided.
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
#' diagnose(fit)
#' }
#'
#' @export
diagnose <- function(object, ...) UseMethod("diagnose")

#' @rdname diagnose
#' @export
diagnose.rfmstate <- function(object, eval_times = NULL,
                              method = c("edge_oob", "cv"), folds = 5L,
                              repeats = 1L, seed = 2026L, g_min = 0.05, ...) {
  method <- match.arg(method)
  if (!inherits(object, "rfmstate")) stop("'object' must be an rfmstate fit.")
  edge_rows <- lapply(names(object$edge_metadata), function(edge) {
    info <- object$edge_metadata[[edge]]
    data.frame(
      transition = edge,
      n_sojourns = info$n_sojourns,
      n_events = info$n_events,
      prediction_error = info$prediction_error,
      oob_concordance = 1 - info$prediction_error,
      endpoint = "binary cause-specific next-exit endpoint",
      stringsAsFactors = FALSE
    )
  })
  edge_oob <- do.call(rbind, edge_rows)
  rownames(edge_oob) <- NULL

  brier <- NULL
  ibs <- NA_real_
  fold_info <- NULL
  integration_interval <- NULL
  if (method == "cv") {
    .positive_integer(folds, "folds")
    .positive_integer(repeats, "repeats")
    if (!is.numeric(g_min) || length(g_min) != 1L || !is.finite(g_min) ||
        g_min <= 0 || g_min >= 1) {
      stop("'g_min' must lie strictly between zero and one.")
    }
    ids <- unique(object$msdata$id)
    if (folds > length(ids)) stop("'folds' cannot exceed the number of subjects.")
    if (is.null(eval_times)) {
      upper <- 0.8 * min(object$max_duration_by_origin)
      eval_times <- seq(0, upper, length.out = 9L)
    }
    if (!is.numeric(eval_times) || length(eval_times) < 2L ||
        anyNA(eval_times) || any(!is.finite(eval_times)) ||
        any(eval_times < 0)) {
      stop("'eval_times' must contain at least two finite nonnegative horizons.")
    }
    eval_times <- sort(unique(as.numeric(eval_times)))
    validation <- .cross_validated_state_brier(
      object, eval_times = eval_times, folds = as.integer(folds),
      repeats = as.integer(repeats), seed = seed, g_min = g_min
    )
    brier <- validation$brier
    ibs <- .trapezoid_integral(brier$time, brier$brier) /
      diff(range(brier$time))
    integration_interval <- range(brier$time)
    fold_info <- validation$fold_info
  }

  structure(
    list(
      method = method,
      edge_oob = edge_oob,
      oob_error = edge_oob[, c("transition", "prediction_error")],
      concordance = data.frame(
        transition = edge_oob$transition,
        c_index = edge_oob$oob_concordance,
        label = "genuine ranger edge OOB concordance",
        stringsAsFactors = FALSE
      ),
      brier = brier,
      ibs = ibs,
      eval_times = if (is.null(brier)) eval_times else brier$time,
      integration_interval = integration_interval,
      folds = fold_info,
      censoring_model = if (method == "cv") "training-fold marginal KM" else NULL,
      g_min = if (method == "cv") g_min else NULL,
      validation_label = if (method == "cv") {
        "patient-level cross-validated full-state probabilities"
      } else {
        "edge-level ranger OOB only"
      }
    ),
    class = "rfmstate_diag"
  )
}

#' @noRd
.subset_msdata <- function(x, keep) {
  out <- x[keep, , drop = FALSE]
  for (nm in c("structure", "initial_state", "covariates",
               "covariate_schema", "metadata")) {
    attr(out, nm) <- attr(x, nm)
  }
  meta <- attr(out, "metadata")
  meta$n_subjects <- length(unique(out$id))
  meta$n_intervals <- nrow(out)
  meta$n_events <- sum(out$status == 1L)
  meta$n_censored <- sum(out$status == 0L)
  states <- attr(out, "structure")$transient
  meta$max_duration_by_origin <- stats::setNames(vapply(states, function(state) {
    value <- out$duration[as.character(out$from) == state]
    if (length(value)) max(value) else NA_real_
  }, numeric(1)), states)
  edges <- attr(out, "structure")$trans_list
  edges$n_events <- vapply(seq_len(nrow(edges)), function(k) {
    sum(out$status == 1L & as.character(out$from) == edges$from[k] &
          as.character(out$to) == edges$to[k], na.rm = TRUE)
  }, integer(1))
  meta$edge_counts <- edges
  attr(out, "metadata") <- meta
  class(out) <- c("msdata", "data.frame")
  out
}

#' @noRd
.cross_validated_state_brier <- function(object, eval_times, folds, repeats,
                                         seed, g_min) {
  ids <- unique(object$msdata$id)
  total_n <- length(ids) * repeats
  accum <- numeric(length(eval_times))
  evaluable <- integer(length(eval_times))
  g_min_seen <- rep(Inf, length(eval_times))
  g_max_seen <- rep(-Inf, length(eval_times))
  fold_records <- list()
  record_index <- 0L

  for (repeat_id in seq_len(repeats)) {
    set.seed(seed + repeat_id - 1L)
    shuffled <- sample(ids, length(ids), replace = FALSE)
    fold_id <- rep(seq_len(folds), length.out = length(ids))
    assignment <- stats::setNames(fold_id, as.character(shuffled))
    for (fold in seq_len(folds)) {
      test_ids <- shuffled[fold_id == fold]
      train_rows <- !(object$msdata$id %in% test_ids)
      test_rows <- object$msdata$id %in% test_ids
      train <- .subset_msdata(object$msdata, train_rows)
      test <- .subset_msdata(object$msdata, test_rows)
      fit_args <- c(list(
        msdata = train,
        covariates = object$covariates,
        num.trees = object$params$num.trees,
        mtry = object$params$mtry,
        min.node.size = object$params$min.node.size,
        min_events = object$params$min_events,
        sparse_warning = Inf,
        importance = object$params$importance,
        seed = seed + repeat_id * 1000L + fold
      ), object$params$ranger_args)
      fold_fit <- do.call(rfmstate, fit_args)
      support <- min(fold_fit$max_duration_by_origin)
      if (max(eval_times) > support) {
        stop("Cross-validation fold ", fold, " repeat ", repeat_id,
             " has support ", format(support),
             " below the prespecified evaluation horizon ", max(eval_times), ".")
      }
      first <- !duplicated(test$id)
      profiles <- test[first, object$covariates, drop = FALSE]
      class(profiles) <- "data.frame"
      pred <- predict(
        fold_fit, newdata = profiles, times = eval_times,
        target_grid_points = 512L, max_grid_points = 8193L,
        grid_tol = 5e-4, check_grid = TRUE
      )
      km <- .fit_censoring_km(train)
      scored <- .score_full_state_fold(test, pred, km, eval_times, g_min)
      accum <- accum + scored$contribution_sum
      evaluable <- evaluable + scored$n_evaluable
      g_min_seen <- pmin(g_min_seen, scored$g)
      g_max_seen <- pmax(g_max_seen, scored$g)
      record_index <- record_index + 1L
      fold_records[[record_index]] <- data.frame(
        repeat_id = repeat_id, fold = fold,
        n_train = length(unique(train$id)),
        n_validation = length(test_ids), support = support,
        seed = seed + repeat_id * 1000L + fold,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    brier = data.frame(
      time = eval_times,
      brier = accum / total_n,
      n_evaluable = evaluable,
      n_validation = total_n,
      G_min = g_min_seen,
      G_max = g_max_seen,
      stringsAsFactors = FALSE
    ),
    fold_info = do.call(rbind, fold_records)
  )
}

#' @noRd
.fit_censoring_km <- function(msdata) {
  terminal <- do.call(rbind, lapply(split(msdata, msdata$id), function(rows) {
    rows <- rows[order(rows$Tstop), , drop = FALSE]
    last <- rows[nrow(rows), , drop = FALSE]
    data.frame(time = last$Tstop, censored = as.integer(last$status == 0L))
  }))
  survival::survfit(survival::Surv(time, censored) ~ 1, data = terminal)
}

#' @noRd
.km_value <- function(km, time, left_limit = FALSE) {
  index <- if (left_limit) sum(km$time < time) else sum(km$time <= time)
  if (index == 0L) 1 else km$surv[index]
}

#' @noRd
.score_full_state_fold <- function(msdata, pred, km, eval_times, g_min) {
  ids <- unique(msdata$id)
  state_names <- pred$structure$state_names
  contribution <- numeric(length(eval_times))
  evaluable <- integer(length(eval_times))
  g_values <- numeric(length(eval_times))
  rows_by_id <- split(msdata, msdata$id)

  for (k in seq_along(eval_times)) {
    time <- eval_times[k]
    g_time <- .km_value(km, time)
    if (!is.finite(g_time) || g_time < g_min) {
      stop("Training-fold censoring survival G(", time, ") = ",
           format(g_time), " is below g_min = ", g_min, ".")
    }
    g_values[k] <- g_time
    for (i in seq_along(ids)) {
      rows <- rows_by_id[[as.character(ids[i])]]
      rows <- rows[order(rows$Tstart), , drop = FALSE]
      absorb <- rows[rows$status == 1L &
                       as.character(rows$to) %in% pred$structure$absorbing,
                     , drop = FALSE]
      absorb_time <- if (nrow(absorb)) min(absorb$Tstop) else Inf
      terminal <- rows[nrow(rows), , drop = FALSE]
      censor_time <- if (terminal$status == 0L) terminal$Tstop else Inf
      weight <- 0
      observed_state <- NA_character_
      if (absorb_time <= time && absorb_time <= censor_time) {
        g_absorb <- .km_value(km, absorb_time, left_limit = TRUE)
        if (!is.finite(g_absorb) || g_absorb < g_min) {
          stop("Training-fold censoring survival before absorption is below g_min.")
        }
        weight <- 1 / g_absorb
        observed_state <- as.character(absorb$to[which.min(absorb$Tstop)])
      } else if (min(absorb_time, censor_time) > time) {
        weight <- 1 / g_time
        prior <- rows[rows$status == 1L & rows$Tstop <= time, , drop = FALSE]
        observed_state <- if (nrow(prior)) {
          as.character(prior$to[nrow(prior)])
        } else {
          pred$initial_state
        }
      }
      if (weight > 0) {
        truth <- as.numeric(state_names == observed_state)
        probability <- pred$state_occ[i, , k]
        contribution[k] <- contribution[k] +
          weight * sum((truth - probability)^2)
        evaluable[k] <- evaluable[k] + 1L
      }
    }
  }
  list(contribution_sum = contribution, n_evaluable = evaluable, g = g_values)
}

#' IPCW survival Brier reference used in focused tests
#' @noRd
.ipcw_survival_brier <- function(obs_time, event, pred_survival, eval_times,
                                 g_min = 0.05) {
  if (!is.matrix(pred_survival) ||
      !identical(dim(pred_survival), c(length(obs_time), length(eval_times)))) {
    stop("'pred_survival' must be subject by evaluation-time matrix.")
  }
  km <- survival::survfit(survival::Surv(obs_time, 1L - event) ~ 1)
  score <- numeric(length(eval_times))
  for (k in seq_along(eval_times)) {
    time <- eval_times[k]
    g_time <- .km_value(km, time)
    if (g_time < g_min) stop("Censoring survival is below g_min.")
    terms <- numeric(length(obs_time))
    failed <- obs_time <= time & event == 1L
    at_risk <- obs_time > time
    if (any(failed)) {
      g_event <- vapply(obs_time[failed], function(value) {
        .km_value(km, value, left_limit = TRUE)
      }, numeric(1))
      terms[failed] <- pred_survival[failed, k]^2 / g_event
    }
    terms[at_risk] <- (1 - pred_survival[at_risk, k])^2 / g_time
    score[k] <- mean(terms)
  }
  score
}

#' @noRd
.trapezoid_integral <- function(x, y) {
  if (length(x) < 2L || length(y) != length(x) || anyNA(y)) return(NA_real_)
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

#' @export
print.rfmstate_diag <- function(x, ...) {
  cat("RFmstate Diagnostics\n")
  cat("  Validation label:", x$validation_label, "\n\n")
  cat("Genuine ranger edge OOB concordance:\n")
  print(x$edge_oob[, c("transition", "n_events", "prediction_error",
                       "oob_concordance")], row.names = FALSE)
  if (!is.null(x$brier)) {
    cat("\nPatient-level cross-validated full-state IPCW Brier score:\n")
    print(x$brier, row.names = FALSE)
    cat("  IBS [", x$integration_interval[1L], ", ",
        x$integration_interval[2L], "] = ", format(x$ibs), "\n", sep = "")
  }
  invisible(x)
}
