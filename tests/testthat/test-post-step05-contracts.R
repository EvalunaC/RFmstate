post_step05_fixture <- function(n = 60L, factor_predictor = FALSE) {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  dat <- data.frame(
    subject = seq_len(n),
    x = seq(-1, 1, length.out = n),
    group = c(rep("common", n - 1L), "heldout_only"),
    time_B = seq_len(n),
    censor = NA_real_,
    duration = seq_len(n) + 100,
    from = "not-a-predictor",
    to = "not-a-predictor",
    Tstart = 0,
    Tstop = seq_len(n),
    status = 1L,
    trans_id = 1L,
    .rfm_hidden = seq_len(n),
    stringsAsFactors = FALSE
  )
  covariates <- if (factor_predictor) c("x", "group") else "x"
  long <- prepare_data(
    dat, "subject", ms, list(B = "time_B"), "censor", covariates
  )
  list(ms = ms, dat = dat, long = long)
}

test_that("prepare_data enforces source roles and stores a predictor contract", {
  fixture <- post_step05_fixture()
  prohibited <- c(
    "subject", "time_B", "censor", "duration", "from", "to",
    "Tstart", "Tstop", "status", "trans_id", ".rfm_hidden"
  )
  for (nm in prohibited) {
    expect_error(
      prepare_data(
        fixture$dat, "subject", fixture$ms, list(B = "time_B"),
        "censor", nm
      ),
      "not baseline covariates|conflicting role",
      info = nm
    )
  }
  expect_error(
    prepare_data(
      fixture$dat, "subject", fixture$ms,
      list(B = c("time_B", "Tstop")), "censor", "x"
    ),
    "exactly one"
  )
  expect_error(
    prepare_data(
      fixture$dat, "subject", fixture$ms, list(B = "time_B"),
      "subject", "x"
    ),
    "conflicting role"
  )
  expect_error(
    prepare_data(
      fixture$dat, "subject", fixture$ms, list(B = "time_B"),
      "censor", NULL
    ),
    "explicit.*character vector"
  )

  contract <- attr(fixture$long, "predictor_contract")
  expect_equal(contract$allowed_names, "x")
  expect_equal(contract$id_source, "subject")
  expect_equal(contract$censor_source, "censor")
  expect_equal(unname(contract$time_sources), "time_B")
  expect_true(is.data.frame(attr(fixture$long, "source_role_table")))
})

test_that("rfmstate accepts only contract predictors and uses a minimal frame", {
  fixture <- post_step05_fixture()
  for (nm in c("id", "from", "to", "Tstart", "Tstop", "duration",
               "status", "trans_id")) {
    expect_error(
      rfmstate(
        fixture$long, covariates = nm, num.trees = 10,
        min.node.size = 3, min_events = 3, sparse_warning = Inf
      ),
      "predictor contract|approved baseline",
      info = nm
    )
  }
  fit <- rfmstate(
    fixture$long, covariates = "x", num.trees = 20,
    min.node.size = 3, min_events = 3, sparse_warning = Inf, seed = 31
  )
  expect_equal(
    fit$edge_metadata[["A->B"]]$model_frame_names,
    c(".rfm_time", ".rfm_event", "x")
  )
  expect_equal(fit$predictor_schema$x$n_distinct, 60L)
  expect_equal(fit$process_assumption, "semi-Markov")
  expect_match(fit$prediction_condition, "fresh entry")
})

test_that("fit-specific schemas recheck baseline constancy within subject", {
  ms <- define_multistate(c("A", "B", "C"), "C", list(A = "B", B = "C"))
  dat <- data.frame(
    subject = seq_len(60), x = seq_len(60),
    time_B = seq_len(60), time_C = seq_len(60) + 1,
    censor = NA_real_
  )
  long <- prepare_data(
    dat, "subject", ms, list(B = "time_B", C = "time_C"), "censor", "x"
  )
  long$x[long$id == 1 & as.character(long$from) == "B"] <- -999
  expect_error(
    rfmstate(
      long, num.trees = 10, min.node.size = 3,
      min_events = 3, sparse_warning = Inf
    ),
    "Baseline predictor 'x'.*not constant within subject.*1"
  )
})

test_that("fit schemas are training-local and CV rejects validation-only levels", {
  fixture <- post_step05_fixture(factor_predictor = TRUE)
  train <- RFmstate:::.subset_msdata(
    fixture$long, fixture$long$id != max(fixture$long$id)
  )
  expect_null(attr(train, "covariate_schema"))
  fold_fit <- rfmstate(
    train, num.trees = 20, min.node.size = 3,
    min_events = 3, sparse_warning = Inf, seed = 32
  )
  expect_equal(fold_fit$predictor_schema$group$levels, "common")
  heldout <- data.frame(x = 0, group = factor("heldout_only"))
  expect_error(predict(fold_fit, heldout, times = c(0, 1)), "Unseen factor")

  full_fit <- rfmstate(
    fixture$long, num.trees = 20, min.node.size = 3,
    min_events = 3, sparse_warning = Inf, seed = 33
  )
  expect_error(
    diagnose(
      full_fit, method = "cv", folds = 2, repeats = 1,
      eval_times = c(0, 1), seed = 34
    ),
    "repeat.*fold.*subject.*group.*heldout_only"
  )
})

test_that("CV returns exact assignments and preserves caller RNG state", {
  fixture <- post_step05_fixture()
  fit <- rfmstate(
    fixture$long, num.trees = 10, min.node.size = 3,
    min_events = 3, sparse_warning = Inf, seed = 35, num.threads = 1
  )
  set.seed(1001)
  before <- .Random.seed
  diag <- suppressWarnings(diagnose(
    fit, method = "cv", folds = 2, repeats = 1,
    eval_times = c(0, 1), seed = 36
  ))
  expect_identical(.Random.seed, before)
  expect_true(all(c(
    "id", "repeat_id", "fold", "assignment_seed", "refit_seed"
  ) %in% names(diag$assignments)))
  expect_equal(nrow(diag$assignments), 60L)
  expect_true(all(c(
    "repeat_id", "fold", "n_train", "n_validation", "support",
    "refit_seed", "min_censoring_survival", "fit_status"
  ) %in% names(diag$fold_summary)))
})

test_that("ranger sampling contract guarantees verified OOB output", {
  fixture <- post_step05_fixture()
  expect_error(
    rfmstate(
      fixture$long, num.trees = 10, min.node.size = 3,
      min_events = 3, sparse_warning = Inf,
      replace = FALSE, sample.fraction = 1
    ),
    "OOB"
  )
  for (value in list(0, -0.1, 1.1, NA_real_, c(0.5, 0.6), "0.5")) {
    expect_error(
      rfmstate(
        fixture$long, num.trees = 10, min.node.size = 3,
        min_events = 3, sparse_warning = Inf, sample.fraction = value
      ),
      "sample.fraction"
    )
  }
  expect_error(
    rfmstate(
      fixture$long, num.trees = 10, min.node.size = 3,
      min_events = 3, sparse_warning = Inf, keep.inbag = FALSE
    ),
    "controlled ranger argument"
  )
  fit <- rfmstate(
    fixture$long, num.trees = 20, min.node.size = 3,
    min_events = 3, sparse_warning = Inf, seed = 37
  )
  info <- fit$edge_metadata[["A->B"]]
  expect_true(is.finite(info$prediction_error))
  expect_equal(info$oob_concordance, 1 - info$prediction_error)
  expect_gt(info$oob_coverage$fraction_with_oob, 0)
  expect_true(isTRUE(info$ranger_arguments$keep.inbag))
})

test_that("solver validates controls and returns only the selected start", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  valid <- list(
    "A->B" = data.frame(time = c(0, 1, 2), hazard = c(0, 0.2, 0.4))
  )
  for (value in list(NA, 0, 1, "yes", c(TRUE, FALSE), list(TRUE))) {
    expect_error(
      compute_trans_prob(valid, ms, times = c(0, 1), check_grid = value),
      "check_grid"
    )
  }
  bad_origin <- valid
  bad_origin[[1]]$hazard[1] <- 0.01
  expect_error(
    compute_trans_prob(bad_origin, ms, times = c(0, 1)),
    "duration zero|H\\(0\\)|time zero"
  )
  expect_warning(
    unchecked <- compute_trans_prob(
      valid, ms, times = c(0, 1), check_grid = FALSE
    ),
    "not passed.*convergence|grid.*not"
  )
  expect_false(unchecked$grid_checked)
  expect_true(is.na(unchecked$grid_converged))

  absorbing <- suppressWarnings(compute_trans_prob(
    valid, ms, times = c(0, 1), start_state = "B", check_grid = FALSE
  ))
  expect_equal(dim(absorbing$entry_prob), c(1L, 2L, 2L))
  expect_equal(dimnames(absorbing$entry_prob)[[1]], "B")
  expect_equal(unname(absorbing$entry_prob["B", "B", ]), c(1, 1))

  tiny <- list(
    "A->B" = data.frame(
      time = c(0, 1, 2), hazard = c(0, 0.2, 0.2 - 1e-12)
    )
  )
  corrected <- suppressWarnings(compute_trans_prob(
    tiny, ms, times = c(0, 1), check_grid = FALSE
  ))
  expect_gt(corrected$hazard_roundoff_corrections, 0)
})

test_that("prediction exposes one start row and normalized method metadata", {
  fixture <- post_step05_fixture()
  fit <- rfmstate(
    fixture$long, num.trees = 20, min.node.size = 3,
    min_events = 3, sparse_warning = Inf, seed = 38
  )
  pred <- predict(
    fit, data.frame(x = 0), times = c(0, 1),
    target_grid_points = 64, max_grid_points = 1025
  )
  expect_equal(dim(pred$entry_prob), c(1L, 1L, 2L, 2L))
  expect_equal(dimnames(pred$entry_prob)[[2]], "A")
  expect_equal(pred$process_assumption, "semi-Markov")
  expect_match(pred$history_summary, "duration since entry")
  expect_match(pred$prediction_condition, "fresh entry")
})

test_that("summaries separate target, competing, and external censoring", {
  ms <- define_multistate(
    c("A", "B", "C"), c("B", "C"), list(A = c("B", "C"))
  )
  dat <- data.frame(
    id = 1:60, x = seq_len(60),
    time_B = c(seq_len(30), rep(NA_real_, 30)),
    time_C = c(rep(NA_real_, 30), seq_len(30)),
    censor = NA_real_
  )
  long <- prepare_data(
    dat, "id", ms, list(B = "time_B", C = "time_C"), "censor", "x"
  )
  fit <- rfmstate(
    long, num.trees = 20, min.node.size = 3,
    min_events = 3, sparse_warning = Inf, seed = 39
  )
  out <- summary(fit)$trans_summary
  expect_true(all(c(
    "n_target_events", "n_competing_exits", "n_external_censored"
  ) %in% names(out)))
  expect_equal(out$n_target_events, c(30L, 30L))
  expect_equal(out$n_competing_exits, c(30L, 30L))
  expect_equal(out$n_external_censored, c(0L, 0L))
})

test_that("seeded public helpers restore caller RNG state", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  set.seed(1002)
  before <- .Random.seed
  invisible(sim_clinical_data(10, structure = ms, seed = 40))
  expect_identical(.Random.seed, before)

  fixture <- post_step05_fixture()
  set.seed(1003)
  before_fit <- .Random.seed
  invisible(rfmstate(
    fixture$long, num.trees = 10, min.node.size = 3,
    min_events = 3, sparse_warning = Inf, seed = 41
  ))
  expect_identical(.Random.seed, before_fit)
})
