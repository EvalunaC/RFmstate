test_that("rfmstate fits models", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 200, structure = ms)
  msdata <- prepare_data(
    data = dat, id = "ID", structure = ms,
    time_map = list(
      Responded = "time_Responded",
      Unresponded = "time_Unresponded",
      Stabilized = "time_Stabilized",
      Progressed = "time_Progressed",
      Death = "time_Death"
    ),
    censor_col = "time_censored",
    covariates = c("age", "sex", "BMI", "treatment")
  )

  fit <- suppressWarnings(rfmstate(
    msdata, covariates = c("age", "sex", "BMI", "treatment"),
    num.trees = 50, min_events = 3, seed = 42
  ))

  expect_s3_class(fit, "rfmstate")
  expect_true(length(fit$models) > 0)
  expect_equal(fit$covariates, c("age", "sex", "BMI", "treatment"))
})

test_that("predict.rfmstate returns valid predictions", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 200, structure = ms)
  msdata <- prepare_data(
    data = dat, id = "ID", structure = ms,
    time_map = list(
      Responded = "time_Responded",
      Unresponded = "time_Unresponded",
      Stabilized = "time_Stabilized",
      Progressed = "time_Progressed",
      Death = "time_Death"
    ),
    censor_col = "time_censored",
    covariates = c("age", "sex", "BMI", "treatment")
  )

  fit <- suppressWarnings(rfmstate(
    msdata, covariates = c("age", "sex", "BMI", "treatment"),
    num.trees = 50, min_events = 3, seed = 42
  ))

  newdata <- data.frame(age = 60, sex = 1, BMI = 25, treatment = 1)
  horizon <- min(fit$max_duration_by_origin)
  prediction_times <- c(0, horizon / 3, 2 * horizon / 3)
  pred <- predict(fit, newdata = newdata, times = prediction_times,
                  target_grid_points = 1024)

  expect_s3_class(pred, "rfmstate_pred")
  expect_equal(pred$n_subjects, 1)
  expect_equal(length(pred$time), 3)

  # State occupation probabilities should be valid
  for (k in seq_along(pred$time)) {
    occ <- pred$state_occ[1, , k]
    expect_true(all(occ >= -1e-8))
    expect_true(abs(sum(occ) - 1) < 1e-8)
  }
})

test_that("summary.rfmstate works", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 200, structure = ms)
  msdata <- prepare_data(
    data = dat, id = "ID", structure = ms,
    time_map = list(
      Responded = "time_Responded",
      Unresponded = "time_Unresponded",
      Stabilized = "time_Stabilized",
      Progressed = "time_Progressed",
      Death = "time_Death"
    ),
    censor_col = "time_censored",
    covariates = c("age", "sex", "BMI", "treatment")
  )

  fit <- suppressWarnings(rfmstate(
    msdata, covariates = c("age", "sex", "BMI", "treatment"),
    num.trees = 50, min_events = 3, seed = 42
  ))

  s <- summary(fit)
  expect_s3_class(s, "summary.rfmstate")
  expect_true(nrow(s$trans_summary) > 0)
  expect_output(print(s), "Random Forest Multistate")
})

test_that("importance.rfmstate works", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 200, structure = ms)
  msdata <- prepare_data(
    data = dat, id = "ID", structure = ms,
    time_map = list(
      Responded = "time_Responded",
      Unresponded = "time_Unresponded",
      Stabilized = "time_Stabilized",
      Progressed = "time_Progressed",
      Death = "time_Death"
    ),
    censor_col = "time_censored",
    covariates = c("age", "sex", "BMI", "treatment")
  )

  fit <- suppressWarnings(rfmstate(
    msdata, covariates = c("age", "sex", "BMI", "treatment"),
    num.trees = 50, min_events = 3, seed = 42
  ))

  imp <- importance(fit)
  expect_s3_class(imp, "rfmstate_importance")
  expect_true(nrow(imp$importance) > 0)
  expect_true(ncol(imp$importance_matrix) > 0)
})

test_that("diagnose.rfmstate works", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 200, structure = ms)
  msdata <- prepare_data(
    data = dat, id = "ID", structure = ms,
    time_map = list(
      Responded = "time_Responded",
      Unresponded = "time_Unresponded",
      Stabilized = "time_Stabilized",
      Progressed = "time_Progressed",
      Death = "time_Death"
    ),
    censor_col = "time_censored",
    covariates = c("age", "sex", "BMI", "treatment")
  )

  fit <- suppressWarnings(rfmstate(
    msdata, covariates = c("age", "sex", "BMI", "treatment"),
    num.trees = 50, min_events = 3, seed = 42
  ))

  diag <- diagnose(fit)
  expect_s3_class(diag, "rfmstate_diag")
  expect_true(nrow(diag$oob_error) > 0)
  expect_true(nrow(diag$concordance) > 0)
  expect_false("bias_variance" %in% names(diag))
  expect_equal(diag$edge_oob$oob_concordance,
               1 - diag$edge_oob$prediction_error)
})
