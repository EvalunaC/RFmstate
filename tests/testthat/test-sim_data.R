test_that("sim_clinical_data generates valid data", {
  ms <- clinical_states()
  set.seed(123)
  dat <- sim_clinical_data(n = 50)

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 50)
  expect_true(all(c("ID", "age", "sex", "BMI", "treatment") %in% names(dat)))

  # Check that time columns are generated dynamically from structure
  expected_time_cols <- c(
    paste0("time_", ms$state_names[-1]),
    "time_censored"
  )
  expect_true(all(expected_time_cols %in% names(dat)))

  # Every patient should have at least one non-NA event time or censoring
  has_event_or_cens <- apply(
    dat[, expected_time_cols, drop = FALSE],
    1, function(x) any(!is.na(x))
  )
  expect_true(all(has_event_or_cens))

  # Event times should be positive
  for (col in expected_time_cols) {
    vals <- dat[[col]][!is.na(dat[[col]])]
    if (length(vals) > 0) {
      expect_true(all(vals >= 0), info = paste("Negative times in", col))
    }
  }
})

test_that("sim_clinical_data with custom structure", {
  ms <- define_multistate(
    state_names = c("S1", "S2", "S3"),
    absorbing = "S3",
    transitions = list(S1 = c("S2", "S3"), S2 = c("S3"))
  )
  set.seed(42)
  dat <- sim_clinical_data(n = 20, structure = ms)
  expect_equal(nrow(dat), 20)

  # Columns should match the custom structure, not the clinical defaults
  expect_true("time_S2" %in% names(dat))
  expect_true("time_S3" %in% names(dat))
  expect_true("time_censored" %in% names(dat))
  expect_false("time_Responded" %in% names(dat))

  # At least some patients should have non-NA event times
  expect_true(any(!is.na(dat$time_S2)) || any(!is.na(dat$time_S3)))
})

test_that("sim_clinical_data works with 4-state structure", {
  ms <- define_multistate(
    state_names = c("Healthy", "Sick", "Recovered", "Dead"),
    absorbing = "Dead",
    transitions = list(
      Healthy = c("Sick", "Dead"),
      Sick = c("Recovered", "Dead"),
      Recovered = c("Dead")
    )
  )
  set.seed(99)
  dat <- sim_clinical_data(n = 50, structure = ms)
  expect_equal(nrow(dat), 50)
  expect_true(all(c("time_Sick", "time_Recovered", "time_Dead",
                     "time_censored") %in% names(dat)))

  # Should be able to prepare data from this simulation
  msdata <- prepare_data(
    data = dat, id = "ID", structure = ms,
    time_map = list(
      Sick = "time_Sick",
      Recovered = "time_Recovered",
      Dead = "time_Dead"
    ),
    censor_col = "time_censored",
    covariates = c("age", "sex", "BMI", "treatment")
  )
  expect_s3_class(msdata, "msdata")
  expect_true(nrow(msdata) > 0)
})

test_that("sim_clinical_data is reproducible with seed", {
  dat1 <- sim_clinical_data(n = 20, seed = 42)
  dat2 <- sim_clinical_data(n = 20, seed = 42)
  expect_identical(dat1, dat2)
})

test_that("simulation retains full precision and truncates at censoring", {
  ms <- clinical_states()
  dat <- sim_clinical_data(n = 500, structure = ms, seed = 42)
  event_cols <- paste0("time_", setdiff(ms$state_names, ms$initial_state))
  for (i in seq_len(nrow(dat))) {
    values <- unlist(dat[i, event_cols, drop = FALSE], use.names = FALSE)
    values <- values[!is.na(values)]
    if (!is.na(dat$time_censored[i])) {
      expect_true(all(values < dat$time_censored[i]))
    }
    expect_equal(anyDuplicated(values), 0L)
  }
  all_events <- unlist(dat[event_cols], use.names = FALSE)
  all_events <- all_events[!is.na(all_events)]
  expect_true(any(abs(all_events * 10 - round(all_events * 10)) > 1e-8))
})
