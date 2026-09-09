two_state_fixture <- function(n = 80L, factor_predictor = FALSE) {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  dat <- data.frame(
    id = seq_len(n),
    x = seq(-1, 1, length.out = n),
    group = rep(c("control", "treated"), length.out = n),
    time_B = seq(1, n),
    censor = NA_real_,
    stringsAsFactors = FALSE
  )
  covariates <- if (factor_predictor) c("x", "group") else "x"
  long <- prepare_data(dat, "id", ms, list(B = "time_B"), "censor",
                       covariates)
  list(ms = ms, data = dat, long = long, covariates = covariates)
}

test_that("fit creates every edge and rejects ranger contract violations", {
  fixture <- two_state_fixture()
  fit <- rfmstate(fixture$long, num.trees = 20, min.node.size = 5,
                  min_events = 5, sparse_warning = Inf, seed = 7,
                  num.threads = 1)
  expect_named(fit$models$A, "B")
  expect_equal(names(fit$edge_metadata), "A->B")
  expect_equal(fit$time_scale, "clock-reset")
  expect_equal(fit$initial_state, "A")
  expect_error(
    rfmstate(fixture$long, num.trees = 10, min.node.size = 5,
             min_events = 5, sparse_warning = Inf,
             case.weights = rep(1, nrow(fixture$long))),
    "Unsupported or RFmstate-controlled"
  )
  expect_error(
    rfmstate(fixture$long, num.trees = 10, min.node.size = 5,
             min_events = 5, sparse_warning = 1.5),
    "sparse_warning"
  )
})

test_that("sparse declared edges stop the complete fit", {
  ms <- define_multistate(
    c("A", "B", "C"), c("B", "C"), list(A = c("B", "C"))
  )
  dat <- data.frame(
    id = 1:20, x = rnorm(20),
    time_B = seq_len(20), time_C = NA_real_, censor = NA_real_
  )
  long <- prepare_data(dat, "id", ms,
                       list(B = "time_B", C = "time_C"), "censor", "x")
  expect_error(
    rfmstate(long, num.trees = 10, min.node.size = 3,
             min_events = 2, sparse_warning = Inf),
    "A->C.*target_events=0"
  )
})

test_that("prediction enforces schema and support", {
  fixture <- two_state_fixture(factor_predictor = TRUE)
  fit <- rfmstate(fixture$long, num.trees = 20, min.node.size = 5,
                  min_events = 5, sparse_warning = Inf, seed = 8)
  valid <- data.frame(
    x = 0,
    group = factor("control", levels = c("control", "treated"))
  )
  pred <- predict(fit, valid, times = c(0, 10),
                  target_grid_points = 128, max_grid_points = 2049)
  expect_s3_class(pred, "rfmstate_pred")
  expect_equal(
    unname(apply(pred$entry_prob[1, , , , drop = FALSE], c(2, 4), sum)),
    matrix(1, nrow = 1, ncol = 2), tolerance = 1e-8
  )
  expect_error(predict(fit, valid["group"], times = 1), "Missing covariates")
  bad <- valid
  bad$group <- factor("new", levels = "new")
  expect_error(predict(fit, bad, times = 1), "Unseen factor")
  wrong_class <- transform(valid, group = as.character(group))
  expect_error(predict(fit, wrong_class, times = 1), "must be a factor")
  expect_error(predict(fit, valid, times = 100), "support")
})

test_that("edge concordance is genuine ranger OOB and no bias-variance remains", {
  fixture <- two_state_fixture()
  fit <- rfmstate(fixture$long, num.trees = 30, min.node.size = 5,
                  min_events = 5, sparse_warning = Inf, seed = 9)
  diagnostics <- diagnose(fit)
  model <- fit$models$A$B
  expect_equal(diagnostics$edge_oob$oob_concordance,
               1 - model$prediction.error)
  expect_false("bias_variance" %in% names(diagnostics))
  expect_error(plot(diagnostics, type = "bias_variance"), "arg")
})

test_that("IPCW Brier equals ordinary Brier without censoring", {
  obs_time <- c(1, 2, 3, 4)
  event <- rep(1L, 4)
  eval_times <- c(1.5, 3.5)
  pred <- matrix(c(0.8, 0.7, 0.6, 0.5,
                   0.4, 0.3, 0.2, 0.1), nrow = 4)
  score <- RFmstate:::.ipcw_survival_brier(
    obs_time, event, pred, eval_times, g_min = 0.05
  )
  ordinary <- vapply(seq_along(eval_times), function(k) {
    mean((as.numeric(obs_time > eval_times[k]) - pred[, k])^2)
  }, numeric(1))
  expect_equal(score, ordinary, tolerance = 1e-12)
})

test_that("patient-level CV returns full-state IPCW Brier and IBS", {
  fixture <- two_state_fixture(60)
  fit <- rfmstate(fixture$long, num.trees = 10, min.node.size = 3,
                  min_events = 3, sparse_warning = Inf, seed = 10,
                  num.threads = 1)
  diagnostics <- suppressWarnings(diagnose(
    fit, method = "cv", folds = 2, repeats = 1,
    eval_times = c(0, 5, 10), seed = 11
  ))
  expect_equal(nrow(diagnostics$brier), 3)
  expect_true(all(is.finite(diagnostics$brier$brier)))
  expect_true(is.finite(diagnostics$ibs))
  expect_equal(diagnostics$censoring_model, "training-fold marginal KM")
  expect_error(
    diagnose(fit, method = "cv", folds = 2,
             eval_times = c(1, 1), seed = 11),
    "distinct horizons"
  )
  expect_error(
    diagnose(fit, method = "cv", folds = 2,
             eval_times = c(0, 1), seed = -1),
    "nonnegative integer"
  )
})
