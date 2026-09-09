#' Fit Clock-Reset Random-Forest Multistate Models
#'
#' Fits one ranger cause-specific survival forest for every declared edge of a
#' validated acyclic, non-recurrent multistate structure. Forest response time
#' is duration since fresh entry into the origin state.
#'
#' @param msdata An \code{msdata} object from \code{\link{prepare_data}}.
#' @param covariates Predictor names. \code{NULL} uses the complete stored
#'   predictor contract; an explicit value must be a nonempty unique subset of
#'   that contract.
#' @param num.trees Positive integer number of trees fitted for every edge.
#' @param mtry Positive integer number of predictors considered at each split.
#'   \code{NULL} uses \code{floor(sqrt(p))}, bounded below by one.
#' @param min.node.size Positive integer ranger minimum terminal-node size.
#' @param min_events Explicit technical safeguard for target events per edge.
#'   It is not a universal adequacy threshold.
#' @param sparse_warning Positive integer descriptive event-count threshold for
#'   an edge warning; set \code{Inf} to disable. It is not a fit threshold.
#' @param importance One of \code{"permutation"}, \code{"impurity"},
#'   \code{"impurity_corrected"}, or \code{"none"}. Permutation importance
#'   is recommended for interpretation under the limitations below.
#' @param seed \code{NULL} or one nonnegative integer seed passed to every edge
#'   forest without changing the caller's RNG state.
#' @param ... Supported ranger controls: \code{replace},
#'   \code{sample.fraction}, \code{splitrule}, \code{num.random.splits},
#'   \code{respect.unordered.factors}, \code{num.threads},
#'   \code{save.memory}, \code{max.depth}, \code{always.split.variables},
#'   \code{alpha}, \code{minprop}, and \code{verbose}. Unnamed, conflicting,
#'   and unknown arguments are rejected. Configurations without genuine OOB
#'   observations are rejected.
#'
#' @return An \code{rfmstate} object containing \code{models} for every
#'   declared edge; \code{structure}; selected \code{covariates} and
#'   \code{predictor_schema}; origin-state fitting data; backend event-time
#'   grids; per-edge event, support, ranger-argument, and genuine OOB metadata;
#'   the validated \code{msdata}; fit \code{params}; package versions; verified
#'   OOB coverage; and the clock-reset semi-Markov time-scale contract.
#'
#' @details Every competing exit from an origin state remains an observed exit
#'   for that sojourn, but is coded as a non-target outcome in the binary
#'   cause-specific forest for a particular edge. An unestimable declared edge
#'   stops the whole fit; it is never omitted or represented by zero hazard.
#'   The ranger model frame is constructed anew as \code{.rfm_time},
#'   \code{.rfm_event}, and the contract-approved predictors. Predictor factor
#'   levels and numeric ranges are learned from the actual fitting rows, not
#'   copied from a preparation-time or full-data schema.
#'
#' @section Ranger argument contract:
#' RFmstate controls and rejects duplicate specification of \code{formula},
#' \code{data}, \code{num.trees}, \code{mtry}, \code{min.node.size},
#' \code{importance}, \code{seed}, \code{write.forest}, \code{oob.error}, and
#' \code{keep.inbag}.
#' Unnamed and unknown arguments also fail. The only names accepted through
#' \code{...} are:
#' \describe{
#'   \item{replace, sample.fraction}{Bootstrap/subsampling controls.}
#'   \item{splitrule, num.random.splits, alpha, minprop}{Survival split
#'     controls supported by the installed ranger version.}
#'   \item{respect.unordered.factors}{Unordered-factor handling.}
#'   \item{num.threads, save.memory, verbose}{Computation controls.}
#'   \item{max.depth}{Maximum tree depth.}
#'   \item{always.split.variables}{Predictors always considered for splitting.}
#' }
#' These arguments retain ranger's definitions and are checked against the
#' installed ranger formal arguments before fitting. \code{case.weights},
#' \code{class.weights}, \code{split.select.weights}, response controls, and
#' every other unlisted ranger argument are rejected in this release.
#' Effective sampling defaults are resolved and stored. In particular,
#' \code{replace = FALSE, sample.fraction = 1} is rejected because it leaves no
#' OOB observations. RFmstate forces \code{oob.error = TRUE} and
#' \code{keep.inbag = TRUE}, verifies finite ranger OOB error after every edge
#' fit, and stores per-sojourn OOB-tree coverage.
#'
#' @section Limitations:
#' The fit supports baseline, time-fixed, complete covariates in single-root
#' acyclic non-recurrent data. It does not support left truncation,
#' time-dependent covariates, recurrent visits/cycles, missing fitting
#' covariates, subject-specific frailty, clock-forward hazards, or confidence
#' intervals. Structural, outcome-time, censoring, ID, and arbitrary
#' long-format columns cannot be added as predictors. Every declared edge must
#' meet \code{min_events}; this technical
#' safeguard is not a universal adequacy threshold.
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
  predictor_contract <- attr(msdata, "predictor_contract")
  if (is.null(predictor_contract) ||
      is.null(predictor_contract$allowed_names)) {
    stop("'msdata' lacks a predictor contract; rerun prepare_data().")
  }
  allowed_predictors <- predictor_contract$allowed_names
  covariates <- covariates %||% allowed_predictors
  if (!is.character(covariates) || !length(covariates) || anyNA(covariates) ||
      anyDuplicated(covariates)) {
    stop("At least one unique, nonmissing fitting covariate is required.")
  }
  outside_contract <- setdiff(covariates, allowed_predictors)
  if (length(outside_contract)) {
    stop("Covariate(s) are outside the approved baseline predictor contract: ",
         paste(outside_contract, collapse = ", "), ".")
  }
  missing_covs <- setdiff(covariates, names(msdata))
  if (length(missing_covs)) {
    stop("Predictor-contract columns missing from msdata: ",
         paste(missing_covs, collapse = ", "), ".")
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
  if (!is.numeric(sparse_warning) || length(sparse_warning) != 1L ||
      is.na(sparse_warning) || sparse_warning < 1 ||
      (is.finite(sparse_warning) && sparse_warning != as.integer(sparse_warning))) {
    stop("'sparse_warning' must be a positive integer or Inf.")
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed) || seed < 0 ||
                         seed != as.integer(seed))) {
    stop("'seed' must be NULL or one nonnegative integer.")
  }
  if (!is.null(seed)) seed <- as.integer(seed)
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
  sampling <- .resolve_ranger_sampling_args(dots)
  if ("always.split.variables" %in% dot_names) {
    forced <- dots$always.split.variables
    if (!is.character(forced) || anyNA(forced) ||
        any(!forced %in% covariates)) {
      stop("'always.split.variables' may contain only approved selected predictors.")
    }
  }

  fit_predictors <- .build_fit_schema(msdata, covariates, predictor_contract)
  fit_msdata <- msdata
  for (nm in covariates) fit_msdata[[nm]] <- fit_predictors$data[[nm]]
  schema <- fit_predictors$schema

  origin_data <- stats::setNames(vector("list", length(structure_obj$transient)),
                                 structure_obj$transient)
  failures <- character(0)
  sparse <- character(0)
  edge_precheck <- vector("list", nrow(structure_obj$trans_list))
  names(edge_precheck) <- paste0(structure_obj$trans_list$from, "->",
                                 structure_obj$trans_list$to)
  for (state in structure_obj$transient) {
    rows <- fit_msdata[as.character(fit_msdata$from) == state, , drop = FALSE]
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
    seed = seed, write.forest = TRUE, oob.error = TRUE, keep.inbag = TRUE,
    replace = sampling$replace, sample.fraction = sampling$sample.fraction
  ), dots[setdiff(names(dots), c("replace", "sample.fraction"))])
  for (state in structure_obj$transient) {
    cs_data <- origin_data[[state]]
    models[[state]] <- list()
    event_times[[state]] <- numeric(0)
    for (dest in structure_obj$transitions[[state]]) {
      edge <- paste0(state, "->", dest)
      dest_idx <- match(dest, structure_obj$transitions[[state]])
      status <- as.integer(cs_data$event_type == dest_idx)
      fit_data <- data.frame(
        .rfm_time = cs_data$duration,
        .rfm_event = status,
        cs_data[, covariates, drop = FALSE],
        check.names = FALSE
      )
      formula <- stats::as.formula(
        paste("survival::Surv(.rfm_time, .rfm_event) ~",
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
        oob.error = TRUE,
        keep.inbag = TRUE,
        replace = sampling$replace,
        sample.fraction = sampling$sample.fraction
      ), dots[setdiff(names(dots), c("replace", "sample.fraction"))])
      rf_fit <- .with_local_seed(seed, do.call(ranger::ranger, call_args))
      oob <- .validate_oob_result(rf_fit, nrow(fit_data), edge)
      models[[state]][[dest]] <- rf_fit
      event_times[[state]] <- sort(unique(c(
        event_times[[state]], rf_fit$unique.death.times
      )))
      info <- edge_precheck[[edge]]
      edge_metadata[[edge]] <- c(info, list(
        backend_event_times = rf_fit$unique.death.times,
        prediction_error = rf_fit$prediction.error,
        oob_concordance = 1 - rf_fit$prediction.error,
        oob_coverage = oob,
        model_frame_names = names(fit_data),
        predictor_names = covariates,
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

  support <- metadata$max_duration_by_origin[structure_obj$transient]
  schema_id <- .schema_identifier(schema)
  structure(
    list(
      models = models,
      structure = structure_obj,
      covariates = covariates,
      predictor_contract = predictor_contract,
      predictor_schema = schema,
      predictor_schema_id = schema_id,
      factor_levels = lapply(schema, `[[`, "levels"),
      origin_data = origin_data,
      event_times = event_times,
      edge_metadata = edge_metadata,
      event_counts = metadata$edge_counts,
      max_duration_by_origin = support,
      initial_state = attr(msdata, "initial_state"),
      time_scale = "clock-reset",
      process_assumption = "semi-Markov",
      process_model = "semi-Markov",
      history_summary = paste(
        "current state, duration since entry, and recorded baseline covariates"
      ),
      prediction_condition = paste(
        "fresh entry into the selected starting state at duration zero"
      ),
      msdata = msdata,
      call = cl,
      params = list(
        num.trees = as.integer(num.trees), mtry = mtry,
        min.node.size = as.integer(min.node.size),
        min_events = as.integer(min_events), sparse_warning = sparse_warning,
        importance = importance, seed = seed,
        ranger_args = dots,
        effective_ranger_args = exact_args
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

#' Resolve and validate ranger sampling settings required for genuine OOB use
#' @noRd
.resolve_ranger_sampling_args <- function(dots) {
  replace <- if ("replace" %in% names(dots)) dots$replace else TRUE
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop("'replace' must be one nonmissing logical value.")
  }
  sample_fraction <- if ("sample.fraction" %in% names(dots)) {
    dots$sample.fraction
  } else if (replace) {
    1
  } else {
    0.632
  }
  if (!is.numeric(sample_fraction) || length(sample_fraction) != 1L ||
      is.na(sample_fraction) || !is.finite(sample_fraction) ||
      sample_fraction <= 0 || sample_fraction > 1) {
    stop("'sample.fraction' must be one finite number in (0, 1].")
  }
  if (!replace && identical(as.numeric(sample_fraction), 1)) {
    stop("The ranger configuration replace = FALSE and sample.fraction = 1 ",
         "leaves no OOB observations and conflicts with RFmstate's verified ",
         "OOB diagnostics.")
  }
  list(replace = replace, sample.fraction = as.numeric(sample_fraction))
}

#' Build a schema from the rows used by one actual fit
#' @noRd
.build_fit_schema <- function(data, covariates, predictor_contract) {
  out <- data[, covariates, drop = FALSE]
  schema <- stats::setNames(vector("list", length(covariates)), covariates)
  for (nm in covariates) {
    x <- out[[nm]]
    if (anyNA(x)) stop("Covariate '", nm, "' contains missing values.")
    rows_by_subject <- split(seq_along(x), data$id)
    varies_within_subject <- vapply(rows_by_subject, function(index) {
      length(unique(x[index])) > 1L
    }, logical(1))
    if (any(varies_within_subject)) {
      offenders <- names(rows_by_subject)[varies_within_subject]
      stop(
        "Baseline predictor '", nm, "' is not constant within subject(s): ",
        paste(offenders, collapse = ", "), "."
      )
    }
    original_class <- predictor_contract$original_classes[[nm]] %||% class(x)
    converted_from_character <- "character" %in% original_class
    if (is.character(x)) {
      x <- factor(x, levels = unique(x))
    } else if (is.factor(x)) {
      values <- as.character(x)
      observed <- unique(values)
      fit_levels <- if (is.ordered(x)) {
        levels(x)[levels(x) %in% observed]
      } else {
        observed
      }
      x <- factor(values, levels = fit_levels, ordered = is.ordered(x))
    }
    if (is.numeric(x) && any(!is.finite(x))) {
      stop("Numeric covariate '", nm, "' contains nonfinite values.")
    }
    if (!(is.numeric(x) || is.integer(x) || is.logical(x) || is.factor(x))) {
      stop("Unsupported covariate class for '", nm, "': ",
           paste(class(x), collapse = "/"), ".")
    }
    out[[nm]] <- x
    schema[[nm]] <- list(
      class = class(x),
      original_class = original_class,
      levels = if (is.factor(x)) levels(x) else NULL,
      ordered = is.ordered(x),
      range = if (is.numeric(x)) range(x) else NULL,
      converted_from_character = converted_from_character,
      n_distinct = as.integer(length(unique(x))),
      scope = "fit-specific schema learned from fitting rows only"
    )
  }
  list(data = out, schema = schema)
}

#' Deterministic human-readable schema identifier
#' @noRd
.schema_identifier <- function(schema) {
  fields <- vapply(names(schema), function(nm) {
    item <- schema[[nm]]
    paste(
      nm,
      paste(item$class, collapse = "/"),
      paste(item$levels %||% character(0), collapse = "/"),
      paste(item$range %||% numeric(0), collapse = "/"),
      item$n_distinct,
      sep = ":"
    )
  }, character(1))
  paste0("fit-schema|", paste(fields, collapse = "|"))
}

#' Evaluate seeded code without changing the caller's RNG state
#' @noRd
.with_local_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

#' Verify ranger OOB statistics and per-sojourn tree coverage
#' @noRd
.validate_oob_result <- function(fit, n_observations, edge) {
  error <- fit$prediction.error
  if (!is.numeric(error) || length(error) != 1L || !is.finite(error) ||
      error < 0 || error > 1) {
    stop("Ranger did not return a finite valid OOB prediction error for edge '",
         edge, "'.")
  }
  concordance <- 1 - error
  if (!is.finite(concordance)) {
    stop("Ranger OOB concordance is unavailable for edge '", edge, "'.")
  }
  inbag <- fit$inbag.counts
  if (is.null(inbag)) {
    stop("Ranger did not retain in-bag counts required to verify OOB coverage ",
         "for edge '", edge, "'.")
  }
  matrix_inbag <- if (is.list(inbag)) {
    do.call(cbind, inbag)
  } else {
    as.matrix(inbag)
  }
  if (nrow(matrix_inbag) != n_observations &&
      ncol(matrix_inbag) == n_observations) {
    matrix_inbag <- t(matrix_inbag)
  }
  if (nrow(matrix_inbag) != n_observations || anyNA(matrix_inbag)) {
    stop("Ranger returned invalid in-bag counts for OOB verification on edge '",
         edge, "'.")
  }
  oob_trees <- rowSums(matrix_inbag == 0)
  fraction <- mean(oob_trees > 0)
  if (!is.finite(fraction) || fraction <= 0) {
    stop("Ranger produced no verified OOB observations for edge '", edge, "'.")
  }
  list(
    min_oob_trees = min(oob_trees),
    median_oob_trees = stats::median(oob_trees),
    max_oob_trees = max(oob_trees),
    fraction_with_oob = fraction,
    n_observations = n_observations,
    num_trees = ncol(matrix_inbag)
  )
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

#' @rdname print_rfmstate_objects
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
        info$n_events, " target, ", info$n_competing_exits, " competing, ",
        info$n_external_censored, " externally censored; OOB C = ",
        format(round(info$oob_concordance, 4), nsmall = 4),
        ", OOB coverage = ",
        format(round(info$oob_coverage$fraction_with_oob, 3), nsmall = 3),
        "\n", sep = "")
  }
  invisible(x)
}
