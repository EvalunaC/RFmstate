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
#' @param folds Positive integer number of patient-level folds for
#'   cross-validation.
#' @param repeats Positive integer number of repeated fold assignments.
#' @param seed Nonnegative integer fold-assignment and refit seed.
#' @param g_min Minimum allowed training-fold censoring survival, strictly
#'   between zero and one.
#' @param ... Ignored; validation settings must use the documented arguments.
#'
#' @return An \code{rfmstate_diag} object. Edge OOB tables
#'   (\code{edge_oob}, \code{oob_error}, and \code{concordance}) are always
#'   present. With cross-validation, \code{brier} is a full-state IPCW score
#'   table, \code{ibs} is its trapezoidal integral over
#'   \code{integration_interval}, and fold assignments/support/seeds are stored
#'   in \code{assignments} and \code{fold_summary}; \code{folds} is a
#'   backward-compatible alias for \code{fold_summary}. The object also records
#'   evaluation times, whether their grid was prespecified or exploratory,
#'   censoring model, \code{g_min}, method, and an explicit validation label.
#'
#' @details Ranger OOB concordance is \code{1 - prediction.error} for each
#'   binary cause-specific edge endpoint. It does not validate the assembled
#'   multistate probability vector. Full-state scores use subject-level held-out
#'   predictions and training-fold Kaplan--Meier censoring estimates. The
#'   stronger marginal independent-censoring assumption applies to these KM
#'   weights. Fold assignment occurs before schema construction; unordered
#'   factor levels and numeric ranges are reconstructed from training subjects
#'   only. A held-out-only factor level is an explicit fold failure. No
#'   bias--variance decomposition is provided.
#'
#' @section Limitations:
#' Edge OOB concordance applies only to separate binary cause-specific
#' endpoints and is not a full-pipeline validation score. Cross-validation
#' currently uses a training-fold marginal Kaplan--Meier censoring model and
#' therefore requires marginal independent censoring for the score. Every fold
#' must fit every declared edge and cover all evaluation horizons; failures,
#' sparse edges, or censoring survival below \code{g_min} stop validation. No
#' calibration model, prediction interval, or bias--variance decomposition is
#' returned. Automatically generated evaluation times are exploratory and are
#' labeled as such; confirmatory work should supply prespecified times.
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
      n_target_events = info$n_events,
      n_competing_exits = info$n_competing_exits,
      n_external_censored = info$n_external_censored,
      prediction_error = info$prediction_error,
      oob_concordance = 1 - info$prediction_error,
      oob_fraction = info$oob_coverage$fraction_with_oob,
      min_oob_trees = info$oob_coverage$min_oob_trees,
      median_oob_trees = info$oob_coverage$median_oob_trees,
      max_oob_trees = info$oob_coverage$max_oob_trees,
      replace = info$ranger_arguments$replace,
      sample_fraction = info$ranger_arguments$sample.fraction,
      endpoint = "binary cause-specific next-exit endpoint",
      stringsAsFactors = FALSE
    )
  })
  edge_oob <- do.call(rbind, edge_rows)
  rownames(edge_oob) <- NULL

  brier <- NULL
  ibs <- NA_real_
  fold_info <- NULL
  assignments <- NULL
  integration_interval <- NULL
  evaluation_grid <- NULL
  if (method == "cv") {
    .positive_integer(folds, "folds")
    .positive_integer(repeats, "repeats")
    if (!is.numeric(g_min) || length(g_min) != 1L || !is.finite(g_min) ||
        g_min <= 0 || g_min >= 1) {
      stop("'g_min' must lie strictly between zero and one.")
    }
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed != as.integer(seed)) {
      stop("'seed' must be one nonnegative integer.")
    }
    seed <- as.integer(seed)
    ids <- unique(object$msdata$id)
    if (folds > length(ids)) stop("'folds' cannot exceed the number of subjects.")
    if (is.null(eval_times)) {
      upper <- 0.8 * min(object$max_duration_by_origin)
      eval_times <- seq(0, upper, length.out = 9L)
      evaluation_grid <- "automatic exploratory"
    } else {
      evaluation_grid <- "prespecified"
    }
    if (!is.numeric(eval_times) || length(eval_times) < 2L ||
        anyNA(eval_times) || any(!is.finite(eval_times)) ||
        any(eval_times < 0)) {
      stop("'eval_times' must contain at least two finite nonnegative horizons.")
    }
    eval_times <- sort(unique(as.numeric(eval_times)))
    if (length(eval_times) < 2L) {
      stop("'eval_times' must contain at least two distinct horizons.")
    }
    validation <- .cross_validated_state_brier(
      object, eval_times = eval_times, folds = as.integer(folds),
      repeats = as.integer(repeats), seed = seed, g_min = g_min
    )
    brier <- validation$brier
    ibs <- .trapezoid_integral(brier$time, brier$brier) /
      diff(range(brier$time))
    integration_interval <- range(brier$time)
    fold_info <- validation$fold_info
    assignments <- validation$assignments
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
      fold_summary = fold_info,
      assignments = assignments,
      evaluation_grid = evaluation_grid,
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
               "predictor_contract", "source_role_table", "metadata")) {
    attr(out, nm) <- attr(x, nm)
  }
  attr(out, "covariate_schema") <- NULL
  meta <- attr(out, "metadata")
  meta$n_subjects <- length(unique(out$id))
  meta$original_id_values <- unique(out$id)
  meta$n_intervals <- nrow(out)
  meta$n_events <- sum(out$status == 1L)
  meta$n_censored <- sum(out$status == 0L)
  meta$covariate_schema <- NULL
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
  edges$n_competing_exits <- vapply(seq_len(nrow(edges)), function(k) {
    sum(out$status == 1L & as.character(out$from) == edges$from[k] &
          as.character(out$to) != edges$to[k], na.rm = TRUE)
  }, integer(1))
  edges$n_external_censored <- vapply(seq_len(nrow(edges)), function(k) {
    sum(out$status == 0L & as.character(out$from) == edges$from[k],
        na.rm = TRUE)
  }, integer(1))
  meta$edge_counts <- edges
  meta$origin_sojourns <- stats::setNames(vapply(states, function(state) {
    sum(as.character(out$from) == state)
  }, integer(1)), states)
  meta$origin_person_time <- stats::setNames(vapply(states, function(state) {
    sum(out$duration[as.character(out$from) == state])
  }, numeric(1)), states)
  meta$origin_external_censoring <- stats::setNames(vapply(states, function(state) {
    sum(out$status == 0L & as.character(out$from) == state)
  }, integer(1)), states)
  attr(out, "metadata") <- meta
  class(out) <- c("msdata", "data.frame")
  out
}

#' @noRd
.cross_validated_state_brier <- function(object, eval_times, folds, repeats,
                                         seed, g_min) {
  .with_local_seed(seed, {
    ids <- unique(object$msdata$id)
    total_n <- length(ids) * repeats
    accum <- numeric(length(eval_times))
    evaluable <- integer(length(eval_times))
    g_min_seen <- rep(Inf, length(eval_times))
    g_max_seen <- rep(-Inf, length(eval_times))
    fold_records <- list()
    assignment_records <- list()
    record_index <- 0L

    for (repeat_id in seq_len(repeats)) {
      assignment_seed <- seed + repeat_id - 1L
      set.seed(assignment_seed)
      shuffled <- sample(ids, length(ids), replace = FALSE)
      fold_id <- rep(seq_len(folds), length.out = length(ids))
      assignment_records[[repeat_id]] <- data.frame(
        id = shuffled,
        repeat_id = repeat_id,
        fold = fold_id,
        assignment_seed = assignment_seed,
        refit_seed = seed + repeat_id * 1000L + fold_id,
        stringsAsFactors = FALSE
      )
      for (fold in seq_len(folds)) {
        test_ids <- shuffled[fold_id == fold]
        train_rows <- !(object$msdata$id %in% test_ids)
        test_rows <- object$msdata$id %in% test_ids
        train <- .subset_msdata(object$msdata, train_rows)
        test <- .subset_msdata(object$msdata, test_rows)
        refit_seed <- seed + repeat_id * 1000L + fold
        fit_args <- c(list(
          msdata = train,
          covariates = object$covariates,
          num.trees = object$params$num.trees,
          mtry = object$params$mtry,
          min.node.size = object$params$min.node.size,
          min_events = object$params$min_events,
          sparse_warning = Inf,
          importance = object$params$importance,
          seed = refit_seed
        ), object$params$ranger_args)
        fold_fit <- do.call(rfmstate, fit_args)
        support <- min(fold_fit$max_duration_by_origin)
        if (max(eval_times) > support) {
          stop("Cross-validation repeat ", repeat_id, " fold ", fold,
               " has support ", format(support),
               " below the prespecified evaluation horizon ", max(eval_times), ".")
        }
        first <- !duplicated(test$id)
        profile_ids <- test$id[first]
        profiles <- test[first, object$covariates, drop = FALSE]
        class(profiles) <- "data.frame"
        .validate_cv_profiles(
          profiles, profile_ids, fold_fit, repeat_id = repeat_id, fold = fold
        )
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
        record <- data.frame(
          repeat_id = repeat_id, fold = fold,
          n_train = length(unique(train$id)),
          n_validation = length(test_ids), support = support,
          assignment_seed = assignment_seed,
          refit_seed = refit_seed,
          min_censoring_survival = min(scored$g),
          fit_status = "success",
          failure_reason = NA_character_,
          stringsAsFactors = FALSE
        )
        for (edge in names(fold_fit$edge_metadata)) {
          column <- paste0("events_", make.names(edge))
          record[[column]] <- fold_fit$edge_metadata[[edge]]$n_events
        }
        for (state in names(fold_fit$max_duration_by_origin)) {
          column <- paste0("support_", make.names(state))
          record[[column]] <- fold_fit$max_duration_by_origin[[state]]
        }
        fold_records[[record_index]] <- record
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
      fold_info = do.call(rbind, fold_records),
      assignments = do.call(rbind, assignment_records)
    )
  })
}

#' Reject held-out levels before invoking a fold prediction
#' @noRd
.validate_cv_profiles <- function(profiles, ids, object, repeat_id, fold) {
  for (nm in object$covariates) {
    schema <- object$predictor_schema[[nm]]
    if (!is.null(schema$levels)) {
      values <- as.character(profiles[[nm]])
      unseen <- setdiff(unique(values), schema$levels)
      if (length(unseen)) {
        affected <- unique(as.character(ids[values %in% unseen]))
        stop(
          "Cross-validation repeat ", repeat_id, " fold ", fold,
          " has subject(s) ", paste(affected, collapse = ", "),
          " with predictor '", nm, "' held-out-only level(s): ",
          paste(unseen, collapse = ", "),
          ". Use fewer folds or a scientifically prespecified category pooling rule."
        )
      }
    }
  }
  invisible(TRUE)
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

#' @rdname print_rfmstate_objects
#' @export
print.rfmstate_diag <- function(x, ...) {
  cat("RFmstate Diagnostics\n")
  cat("  Validation label:", x$validation_label, "\n\n")
  cat("Genuine ranger edge OOB concordance:\n")
  print(x$edge_oob[, c(
    "transition", "n_target_events", "n_competing_exits",
    "n_external_censored", "prediction_error", "oob_concordance",
    "oob_fraction", "replace", "sample_fraction"
  )], row.names = FALSE)
  if (!is.null(x$brier)) {
    cat("\nPatient-level cross-validated full-state IPCW Brier score:\n")
    print(x$brier, row.names = FALSE)
    cat("  IBS [", x$integration_interval[1L], ", ",
        x$integration_interval[2L], "] = ", format(x$ibs), "\n", sep = "")
  }
  invisible(x)
}
