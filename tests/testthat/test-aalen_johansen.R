test_that("aalen_johansen computes valid estimates", {
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

  aj <- aalen_johansen(msdata)

  expect_s3_class(aj, "aj_estimate")
  expect_true(length(aj$time) > 0)

  # State occupation probabilities should sum to ~1
  final_occ <- aj$state_occ[nrow(aj$state_occ), ]
  expect_true(abs(sum(final_occ) - 1) < 0.01)

  # All probabilities non-negative
  expect_true(all(aj$state_occ >= -0.01))

  # Transition probability matrices should have rows summing to ~1
  last_P <- aj$trans_prob[[length(aj$trans_prob)]]
  row_sums <- rowSums(last_P)
  expect_true(all(abs(row_sums - 1) < 0.05))
  expect_false("variance" %in% names(aj))
  expect_equal(aj$initial_state, "Baseline")
  expect_error(aalen_johansen(msdata, s = 1), "Only s = 0")
})

test_that("aalen_johansen handles empty data gracefully", {
  ms <- define_multistate(
    state_names = c("A", "B"),
    absorbing = "B",
    transitions = list(A = "B")
  )

  msdata <- data.frame(
    id = character(0), from = character(0), to = character(0),
    Tstart = numeric(0), Tstop = numeric(0),
    status = integer(0), trans_id = integer(0),
    stringsAsFactors = FALSE
  )
  attr(msdata, "structure") <- ms
  class(msdata) <- c("msdata", "data.frame")

  expect_error(aalen_johansen(msdata), "No events")
})

test_that("two-state AJ point estimates agree with Kaplan-Meier", {
  ms <- define_multistate(c("A", "B"), "B", list(A = "B"))
  dat <- data.frame(
    id = 1:8, x = 1:8,
    time_B = c(1, 2, 3, 4, 5, NA, NA, NA),
    censor = c(NA, NA, NA, NA, NA, 2.5, 4.5, 6)
  )
  long <- prepare_data(dat, "id", ms, list(B = "time_B"), "censor", "x")
  aj <- aalen_johansen(long)
  followup <- ifelse(is.na(dat$time_B), dat$censor, dat$time_B)
  status <- as.integer(!is.na(dat$time_B))
  km <- survival::survfit(survival::Surv(followup, status) ~ 1)
  km_at_events <- summary(km, times = aj$time, extend = TRUE)$surv
  expect_equal(aj$state_occ[, "A"], km_at_events, tolerance = 1e-12)
  expect_equal(aj$state_occ[, "B"], 1 - km_at_events, tolerance = 1e-12)
})
