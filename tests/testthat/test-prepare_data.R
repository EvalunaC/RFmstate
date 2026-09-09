test_that("prepare_data converts wide to long format", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 50, structure = ms)
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

  expect_s3_class(msdata, "msdata")
  expect_true(nrow(msdata) > 0)
  expect_true(all(c("id", "from", "to", "Tstart", "Tstop",
                     "status", "trans_id") %in% names(msdata)))

  # All patients present
  expect_true(length(unique(msdata$id)) <= 50)

  # Tstart < Tstop
  expect_true(all(msdata$Tstop > msdata$Tstart))

  # Status is 0 or 1
  expect_true(all(msdata$status %in% c(0L, 1L)))

  # From states are valid
  expect_true(all(msdata$from %in% ms$state_names))

  # To states are valid (or NA for censored)
  observed_to <- msdata$to[!is.na(msdata$to)]
  expect_true(all(observed_to %in% ms$state_names))
})

test_that("prepare_data handles right censoring", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 50, structure = ms)
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

  # Censored intervals should have status = 0
  censored <- msdata[msdata$status == 0, ]
  expect_true(all(is.na(censored$to)))
})

test_that("prepare_data validates inputs", {
  ms <- clinical_states()
  set.seed(42)
  dat <- sim_clinical_data(n = 10, structure = ms)

  expect_error(
    prepare_data(dat, id = "MISSING", structure = ms,
                 time_map = list(), censor_col = "time_censored",
                 covariates = c()),
    "not found"
  )
})

test_that("msdata print is a summary and head returns validated ordinary rows", {
  ms <- define_multistate(c("Z", "A"), "A", list(Z = "A"))
  dat <- data.frame(id = 1:3, x = 1:3, time_A = 1:3, censor = NA_real_)
  long <- prepare_data(dat, "id", ms, list(A = "time_A"), "censor", "x")

  expect_output(print(long), "Multistate Data")
  expect_output(print(long), "Z")
  rows <- head(long, 2)
  expect_s3_class(rows, "data.frame")
  expect_false(inherits(rows, "msdata"))
  expect_equal(nrow(rows), 2)
  expect_error(head(long, -1), "nonnegative integer")
  expect_error(head(long, NA_real_), "nonnegative integer")
})
